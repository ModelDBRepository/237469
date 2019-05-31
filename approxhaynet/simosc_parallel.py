from mpi4py import MPI
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import time
import scipy.io
import pickle
import sys
import mutation_stuff
import approxhaynetstuff
import mytools


def simseedburst_func(Nmc=1, tstop=10200,mutID=0,rdSeed=1,Econ=0.0004,Icon=0.001,nseg=5,rateCoeff=1.0,gNoiseCoeff=1.0,gSynCoeff=1.0,oscfreq=1.0,phase=0.0,oscamp=0.25,Ncells2save=1,sparsedt=1.0,Nsyns2save=1,connM=[],gToBlock=[],blockEfficiency=0.0):
  myrandsavemat = int(100000*gSynCoeff+1000*Nmc+rdSeed)

  rank = MPI.COMM_WORLD.Get_rank()
  dt = 0.025

  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
  MT = mutation_stuff.getMT(False)
  defVals = mutation_stuff.getdefvals()
  keyList = defVals.keys()
  for idefval in range(0,len(keyList)):
    if type(defVals[keyList[idefval]]) is not list:
      defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
  updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
  whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
  unpicklefile = open('scalings_cs.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  theseCoeffsAllAll = unpickledlist[0]

  filename = 'pars_withmids_combfs_final'
  unpicklefile = open(filename+".sav", 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  par_names = unpickledlist[0]
  par_values = unpickledlist[1]
  paramdict = {}
  for i in range(0,len(par_names)):
    paramdict[par_names[i]] = par_values[i]

  unpicklefile = open('isidists_'+str(oscfreq)+'_'+str(oscamp)+'.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  probsE = numpy.array(unpickledlist[0])
  probsI = numpy.array(unpickledlist[1])
  tmax = unpickledlist[2]
  nphases = probsE.shape[0]
  nts = probsE.shape[1]
  cprobsE = numpy.cumsum(probsE,1)
  cprobsI = numpy.cumsum(probsI,1)

  h("""
{load_file("stdlib.hoc")}
{load_file("stdrun.hoc")}

initialization_tstart = startsw()

strdef fileName
objref fileObj

fileObj = new File()

rdSeed = """+str(rdSeed)+"""
Nmc = """+str(Nmc)+"""
connectivity = 1
nphases = """+str(nphases)+"""
nts = """+str(nts)+"""
tmax = """+str(tmax)+"""
oscfreq = """+str(oscfreq)+"""
pi = """+str(numpy.pi)+"""
phase = """+str(phase)+"""
pi2oscfreq001 = 2*pi*oscfreq*0.001
tmaxpernts = tmax/nts
nphasesper2pi = nphases/(2*pi)

tstop = """+str(tstop)+"""
rcpWeightFactor = 1.5 // the factor by which reciprocal weights are stronger than unidirectional weights
pT2Tr = 0.06 //probability of reciprocating an existing connection to another L5bPC
pT2T = 0.13 //probability of a L5bPC being connected to another L5bPC
Econ = """+str(Econ)+""" //excitatory synaptic conductance
Icon = """+str(Icon)+""" //inhibitory synaptic conductance
NcontE = 5 // number of excitatory synaptic contacts per connection
NsynE = 10000 // number of excitatory synapses
NsynI = 2500 // number of inhibitory synapses
gNoiseCoeff = """+str(gNoiseCoeff)+""" // scaling of background synaptic conductances
rateE = """+str(0.72*rateCoeff)+""" // average rate of presynaptic excitatory cells
rateI = """+str(7.0*rateCoeff)+""" // average rate of presynaptic inhibitory cells
mainBifurcation = 650

{Ncells2save = """+str(Ncells2save)+"""}
sparsedt = """+str(sparsedt)+""" // recordings of [Ca], I_SK and vApical are done with low temporal resolution to save memory
{gSynCoeff = """+str(gSynCoeff)+"""}

objref tempvec
tempvec = new Vector()
{tempvec.append(Nmc)}
{Ncells2save = """+str(Ncells2save)+"""}
""")
  print "Params OK!"

  h("""
{load_file("models/TTC_det.hoc")}

objref MC_TTC
objref sl //synaptic locations list

objref rds1,rds2
{rds1 = new Random(1000*rdSeed)}
{rds1.uniform(0,1)} //random for microcircuit connectivity and noisyst

objref conMat
conMat = new Matrix(Nmc,Nmc)

for(i=0;i<Nmc;i+=1){
        conMat.x[i][i]=0
}
""")
  if len(connM) == 0:
    h("""
for(i=0;i<(Nmc-2);i+=1){
        for(j=(i+1);j<Nmc;j+=1){
                if (connectivity){
                        pcon = rds1.repick()
                        if (pcon<pT2Tr){
                                conMat.x[i][j]=rcpWeightFactor*gSynCoeff
                                conMat.x[j][i]=rcpWeightFactor*gSynCoeff
                        } else {
                                if (pcon<(pT2Tr + 0.5*pT2T)){
                                        conMat.x[i][j]=gSynCoeff
                                        conMat.x[j][i]=0
                                } else {
                                        if (pcon<(pT2Tr + pT2T)){
                                                conMat.x[i][j]=0
                                                conMat.x[j][i]=gSynCoeff
                                        } else {
                                                conMat.x[i][j]=0
                                                conMat.x[j][i]=0
                                        }
                                }
                        }
                } else {
                        conMat.x[i][j]=0
                        conMat.x[j][i]=0
                }
        }
}
""")
  else:
    for i in range(0,Nmc):
      for j in range(0,Nmc):
        if connM[i][j]:
          h("conMat.x["+str(i)+"]["+str(j)+"]="+str(gSynCoeff*connM[i][j])) # Remember that conMat.x[i][j] is connection FROM j TO i
  print "Connectivity OK!"
  h("double cprobsE["+str(nphases)+"]["+str(nts)+"]")
  h("double cprobsI["+str(nphases)+"]["+str(nts)+"]")
  print "Copying probability tables..."
  for i in range(0,nphases):
    for j in range(0,nts):
      h("cprobsE["+str(i)+"]["+str(j)+"] = "+str(cprobsE[i,j]))
      h("cprobsI["+str(i)+"]["+str(j)+"] = "+str(cprobsI[i,j]))

  del(probsE)
  del(probsI)
  del(cprobsE)
  del(cprobsI)


  #Define HOC functions for randomly sampling from given PDFs (that are defined in arrays cprobsE and cprobsI)
  #$1: which phase to use
  h("""func picktE() {local imin,imax,inow,ii
  imin = 0
  imax = nts
  inow = int(nts/2)
  r = rds2.repick()
  for(ii=0;ii<int(log(nts)/log(2))+1;ii+=1) {
    if (r > cprobsE[$1][inow]) {
      imin = inow
      inow = int((imax+imin)/2)
    } else {
      imax = inow
      inow = int((imax+imin)/2)
    }
    if (imin==inow || imax==inow) {
      break
    }
  }
  if (r > cprobsE[$1][inow]) {
    inow = inow+1
  }
  return inow
}
""")

  #$1: which phase to use
  h("""func picktI() {local imin,imax,inow,ii
  imin = 0
  imax = nts
  inow = int(nts/2)
  r = rds2.repick()
  for(ii=0;ii<int(log(nts)/log(2))+1;ii+=1) {
    if (r > cprobsI[$1][inow]) {
      imin = inow
      inow = int((imax+imin)/2)
    } else {
      imax = inow
      inow = int((imax+imin)/2)
    }
    if (imin==inow || imax==inow) {
      break
    }
  }
  if (r > cprobsI[$1][inow]) {
    inow = inow+1
  }
  return inow
}
""")

  h("""
{load_file("netparmpi.hoc")}
objref epnm

epnm = new ParallelNetManager(Nmc)
{epnm.round_robin()}

strdef treename
objref NsynsE, NsynsI, preTrainList, preTrainIndexList
""")
  print "ParallelNetManager OK!"

  for i in range(0,Nmc):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
        print \"rank = """+str(rank)+""", gid \", i, \" exists\\n\"
        MC_TTC = new TTC()
        epnm.register_cell(i,MC_TTC)
        epnm.pc.gid2cell(i).initRand(1000*rdSeed+i)
        epnm.pc.gid2cell(i).setnetworkparameters(rcpWeightFactor,Econ,Icon,NsynE,NsynI,NcontE,1.0,1.0,1.0,gNoiseCoeff)
  }""")
    print "ParallelNetManager " + str(i)+ " OK!"

  approxhaynetstuff.setparams(paramdict,Nmc)

  for i in range(0,Nmc):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
    {rds2 = new Random(1000*rdSeed+i)}//random for presynaptic trains
    {rds2.uniform(0,1)} //random for microcircuit noisyst
    {preTrainList = new List()}

    for(i2=0;i2<(NsynE+NsynI);i2+=1){
      {preTrainList.append(new Vector())}
      pst=0 //presynaptic spike time
      while(pst < tstop){
        if (i2<NsynE) {
          pst+= tmaxpernts*(0.5+picktE(int(nphasesper2pi*(pst*pi2oscfreq001+phase))%nphases))
        } else {
          pst+= tmaxpernts*(0.5+picktI(int(nphasesper2pi*(pst*pi2oscfreq001+phase))%nphases))
        }
        //print \"rank = """+str(rank)+""", i = \", i, \", pst = \", pst, \"\\n\"
        {preTrainList.o[preTrainList.count()-1].append(pst)}
      }
    }
    {epnm.pc.gid2cell(i).distributeSyn()}                            
    {epnm.pc.gid2cell(i).distributeSyn2()}                           
    {epnm.pc.gid2cell(i).setpretrains(preTrainList)}                 
    {epnm.pc.gid2cell(i).queuePreTrains()}                           
    //print \"distributeSyn(), setpretrains(), queuePreTrains() OK! i=\", i, \".\"
    for(i2=0;i2<preTrainList.count();i2+=1){
      preTrainList.o[i2].resize(0)
    }
    preTrainList.remove_all()
  }
""")
    if i in [5,10,20,50,100,200,500,1000,2000,5000,10000]:
      print "Synapses set for "+str(i)+" neurons"

  print "Spike trains OK!"

  #LOAD MUTATIONS
  counter = -1
  found = 0
  for igene in range(0, len(MT)):
    for imut in range(0, len(MT[igene])):
      theseCoeffs = theseCoeffsAllAll[0][igene][imut]

      nVals = len(MT[igene][imut])*[0]
      thesemutvars = []
      for imutvar in range(0,len(MT[igene][imut])):
        thesemutvars.append(MT[igene][imut][imutvar][0])
        if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
          MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
        nVals[imutvar] = len(MT[igene][imut][imutvar][1])
      cumprodnVals = cumprod(nVals)
      allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
      allmutvals = []
      for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
        allmutvals.append([0]*len(thesemutvars))
      for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
        for imutvar in range(0,len(MT[igene][imut])):
          if imutvar==0:
            allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
          else:
            allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]

      for iallmutval in range(0, len(theseCoeffs)):
        if igene == 0 and imut == 0 and iallmutval == 0:
          iters = [-1,0,2,6,8]
        else:
          iters = [0,2,6,8]
        for iiter in range(0,len(iters)):
          iter = iters[iiter]
          counter = counter + 1
          if counter == mutID:
            found = 1
            break
        if found:
          break
      if found:
        break
    if found:
      break
  if not found:
    print "Not found corresponding mutID, mutID = " + str(mutID)
    sys.exit()

  if iter >= 0:
    thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
  else:
    thisCoeff = 0
  print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)

  mutText = ""
  for imutvar in range(0,len(MT[igene][imut])):
    if imutvar > 0 and imutvar%2==0:
      mutText = mutText+"\n"
    mutvars = allmutvars[iallmutval][imutvar]
    mutvals = allmutvals[iallmutval][imutvar]
    if type(mutvars) is str:
      mutvars = [mutvars]
    mutText = mutText + str(mutvars) + ": "
    for kmutvar in range(0,len(mutvars)):
      mutvar = mutvars[kmutvar]
      if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
        newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
        if mutvals >= 0 and kmutvar==0:
          mutText = mutText + "+" + str(mutvals) +" mV"
        elif kmutvar==0:
          mutText = mutText  + str(mutvals) +" mV"
      else:
        newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
        if kmutvar==0:
          mutText = mutText + "*" + str(mutvals)
      if kmutvar < len(mutvars)-1:
        mutText = mutText + ", "
      if mutvar.find('_Ih') > -1:
        updateThese = [1,1,1]
      elif mutvar.find('_Ca_HVA') > -1 or mutvar.find('_Ca_LVAst') > -1 or mutvar.find('_SKv3.1') > -1 or mutvar.find('_Ca_HVA') > -1 or mutvar.find('_SK_E2') > -1 or mutvar.find('_NaTa_t') > -1 or mutvar.find('_CaDynamics_E2') > -1:
        updateThese = [1,1,0]
      elif mutvar.find('_K_Pst') > -1 or mutvar.find('_K_Tst') > -1 or mutvar.find('_Nap_Et2') > -1:
        updateThese = [1,0,0]
      elif mutvar.find('_Im') > -1:
        updateThese = [0,1,0]
      else:
        print "Error: str=" + str(mutvar)
        updatedThese = [0,0,0]
      for iupdated in range(0,3):
        if updateThese[iupdated]:
          print """forsec epnm.pc.gid2cell(i)."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""" 
          for i in range(0,Nmc):
            h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  forsec epnm.pc.gid2cell(i)."""+str(updatedVars[iupdated])+""" {
    """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }
}""")

  print mutText
  h("""
thisCa = 0.0001
for(i=0;i<Nmc;i+=1){
  if (epnm.gid_exists(i)) {
    thisCa = epnm.pc.gid2cell(i).soma.minCai_CaDynamics_E2
  }
}
""")
  thisCa = h.thisCa

  myMechs = ['Ca_HVA','Ca_LVAst','Ih','Im','K_Pst','K_Tst','NaTa_t','Nap_Et2','SK_E2','SKv3_1','']
  myMechToAdd = ""
  for iblock in range(0,len(gToBlock)):
    for iMech in range(0,len(myMechs)):
      if gToBlock[iblock] in myMechs[iMech]:
        break
    if iMech <= 9:
      if type(blockEfficiency) is list:
        for i in range(0,Nmc):
          h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  epnm.pc.gid2cell(i).soma g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[0])+"""
  forsec epnm.pc.gid2cell(i).apical { g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[1])+""" }
}""")
        print("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  epnm.pc.gid2cell(i).soma g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[0])+"""
  forsec epnm.pc.gid2cell(i).apical { g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" = g"""+str(myMechs[iMech])+"""bar_"""+str(myMechs[iMech])+""" * """+str(blockEfficiency[1])+""" }
}""")
        myMechToAdd = myMechToAdd+myMechs[iMech]+'x'+str(blockEfficiency[0])+'-'+str(blockEfficiency[1])+'_'        
      else:
        print("forall if(ismembrane(\""+str(myMechs[iMech])+"\")) { g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" = g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" * "+str(blockEfficiency)+" }")
        h("forall if(ismembrane(\""+str(myMechs[iMech])+"\")) { g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" = g"+str(myMechs[iMech])+"bar_"+str(myMechs[iMech])+" * "+str(blockEfficiency)+" }")
        myMechToAdd = myMechToAdd+myMechs[iMech]+'x'+str(blockEfficiency)+'_'
    else:
      print "Error: No mechanism recognized"


  h("""
v_init = -80
cai0_ca_ion = thisCa
dt = """+str(dt)+"""
objref syninds, conMatRows

for(i=0;i<Nmc;i+=1){
  if (epnm.gid_exists(i)) {
    //epnm.pc.gid2cell(i).geom_nseg()
    epnm.pc.gid2cell(i).insertMCcons(conMat.getcol(i))
  }
}

{syninds = new Vector()}
{conMatRows = new List()}
for(i=0;i<Nmc;i+=1){
        syninds.append(2*3*"""+str(nseg)+""")
}

// appending the microcircuit connections
for(i=0;i<Nmc;i+=1){
        conMatRows.append(new Vector())
        for(j=0;j<Nmc;j+=1){
                conMatRows.o[i].insrt(j,conMat.x[j][i])
                if (conMat.x[j][i] != 0){
                        for(jj=0;jj<NcontE;jj+=1){
                                epnm.nc_append(j,i,syninds.x[i],1,0.5)
                                syninds.x[i] +=1
                        }
                }
        }
}
""")
  h("forall nseg="+str(nseg))
  print "Syninds OK!"

  h("""
objref vSomaList, tvecList, caSomaList, skSomaList, cahvaSomaList, calvaSomaList
objref natSomaList, napSomaList, ihSomaList, kv31SomaList, ktSomaList, kpSomaList, IList
objref apcvecList, apcList, netcon, nil, spikes, spikedCells
{spikes = new Vector()}
{spikedCells = new Vector()}

{apcvecList = new List()}
{apcList = new List()}
{vSomaList = new List()}
{caSomaList = new List()}
{skSomaList = new List()}
{cahvaSomaList = new List()}
{calvaSomaList = new List()}
{natSomaList = new List()}
{napSomaList = new List()}
{ihSomaList = new List()}
{kv31SomaList = new List()}
{ktSomaList = new List()}
{kpSomaList = new List()}
{IList = new List()}
{tvecList = new List()}
""")

  Nsyns = numpy.array(h.syninds)-2*3*nseg
  cumpNsyns = numpy.cumsum(Nsyns)/numpy.sum(Nsyns)
  randVec = [rand() for x in range(0,Nsyns2save)]
  if sum(Nsyns) > 0:
    cellsSynRecorded = [next(i for i,x in enumerate(cumpNsyns) if x > randVec[j]) for j in range(0,Nsyns2save)]
    synIndsRecorded = [int(2*3*nseg+rand()*Nsyns[i]) for i in cellsSynRecorded]
  else:
    cellsSynRecorded = []
    synIndsRecorded = []

  for i in range(0,int(h.Ncells2save)):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
    {vSomaList.append(new Vector())}
    {caSomaList.append(new Vector())}
    {skSomaList.append(new Vector())}
    {cahvaSomaList.append(new Vector())}
    {calvaSomaList.append(new Vector())}
    {natSomaList.append(new Vector())}
    {napSomaList.append(new Vector())}
    {ihSomaList.append(new Vector())}
    {kv31SomaList.append(new Vector())}
    {ktSomaList.append(new Vector())}
    {kpSomaList.append(new Vector())}
    {tvecList.append(new Vector())}

    access epnm.pc.gid2cell(i).soma
    {vSomaList.o[vSomaList.count()-1].record(&v(0.5),dt)}
    {caSomaList.o[caSomaList.count()-1].record(&cai(0.5),sparsedt)}
    {skSomaList.o[skSomaList.count()-1].record(&ik_SK_E2(0.5),sparsedt)}
    {cahvaSomaList.o[skSomaList.count()-1].record(&ica_Ca_HVA(0.5),sparsedt)}
    {calvaSomaList.o[skSomaList.count()-1].record(&ica_Ca_LVAst(0.5),sparsedt)}
    {natSomaList.o[skSomaList.count()-1].record(&ina_NaTa_t(0.5),sparsedt)}
    {napSomaList.o[skSomaList.count()-1].record(&ina_Nap_Et2(0.5),sparsedt)}
    {ihSomaList.o[skSomaList.count()-1].record(&ihcn_Ih(0.5),sparsedt)}
    {kv31SomaList.o[skSomaList.count()-1].record(&ik_SKv3_1(0.5),sparsedt)}
    {ktSomaList.o[skSomaList.count()-1].record(&ik_K_Tst(0.5),sparsedt)}
    {kpSomaList.o[skSomaList.count()-1].record(&ik_K_Pst(0.5),sparsedt)}
  }
""")
    indSynIndsRecorded = [ix for ix,x in enumerate(cellsSynRecorded) if x==i]
    for isyn in range(0,len(indSynIndsRecorded)):
      h("""
  if (epnm.gid_exists(i)) {
    {IList.append(new Vector())}
    {IList.o[IList.count()-1].record(&epnm.pc.gid2cell(i).synlist.o["""+str(synIndsRecorded[indSynIndsRecorded[isyn]])+"""].i, sparsedt)}
  }
""")
    if i < Ncells2save:
      h("""
  if (epnm.gid_exists(i)) {
    {vSomaList.append(new Vector())}
    {vSomaList.o[vSomaList.count()-1].record(&v(0.5),dt)}
  }
""")

  for i in range(0,Nmc):
    h("""
  i = """+str(i)+"""
  if (epnm.gid_exists(i)) {
          access epnm.pc.gid2cell(i).soma
          {apcList.append(new APCount(0.5))}
          {apcvecList.append(new Vector())}
          apcList.o[apcList.count()-1].thresh= -40
          {apcList.o[apcList.count()-1].record(apcvecList.o[apcList.count()-1])}
          {netcon = new NetCon(&v(0.5), nil)}
          netcon.threshold = -20
          {netcon.record(spikes, spikedCells)  }
}""")
    print "epnm.gid " + str(i)+ " OK!"

  h("""
{epnm.set_maxstep(100)}

stdinit()

if (epnm.gid_exists(0)) {
        print \"\\n\"
        sim_tstart = startsw()
        initializationtime = (sim_tstart-initialization_tstart)/3600
        print \"Initialization completed. Initialization took \", initializationtime, \" hours\\n\"
        print \"Starting simulation\\n\"
        print \"\\n\"
}
""")
  h("""
{epnm.psolve(tstop)}

if (epnm.gid_exists(0)) {
        simruntime = (startsw() - sim_tstart)/3600
        print \"Simulation took \", simruntime, \" hours\\n\"
}

""")
  print "Simulation OK!"

  spikes = numpy.array(h.spikes)
  spikedCells = numpy.array(h.spikedCells)
  vSoma = numpy.array(h.vSomaList)
  picklelist = [spikes,spikedCells,vSoma]
  file = open('spikes_parallel_osc'+str(oscfreq)+'_'+myMechToAdd+str(nseg)+'_'+str(tstop)+'_'+str(mutID)+'_'+str(Econ)+'_'+str(Icon)+'_'+str(rateCoeff)+'_'+str(gNoiseCoeff)+'_'+str(gSynCoeff)+'_'+str(rdSeed)+'_'+str(rank)+'_of_'+str(Nmc)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()
  print 'Saved to spikes_parallel_osc'+str(oscfreq)+'_'+myMechToAdd+str(nseg)+'_'+str(tstop)+'_'+str(mutID)+'_'+str(Econ)+'_'+str(Icon)+'_'+str(rateCoeff)+'_'+str(gNoiseCoeff)+'_'+str(gSynCoeff)+'_'+str(rdSeed)+'_'+str(rank)+'_of_'+str(Nmc)+'.sav'

  h("""
{epnm.pc.runworker()}
{epnm.pc.done()}
""")
  print "runworkers OK!"

  return [spikes,spikedCells,vSoma]


