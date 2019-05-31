from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random

random.seed(1)

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
BACdt = 5.0
fs = 8
xs = range(700,1150,50);
tstop = 11000.0
currCoeff = 0.9 # use the threshold current I*1.1 for inducing the first spike, and I*1.1*2 for the second spike 
ITERS = 20#
PPIdts = [40,50,60,70,80,100,500]
currgs = [0.000025, 0.00005, 0.000075, 0.0001, 0.000125, 0.00015]
maxLens = [1300,1185]

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
defValsOrig = mutation_stuff.getdefvals()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
    defValsOrig[keyList[idefval]] = [defValsOrig[keyList[idefval]], defValsOrig[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments
unpicklefile = open('scalings_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
theseMutValsAllAll = unpickledlist[2]

unpicklefile = open('thresholddistalamp300_control.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
gs_control = unpickledlist[0]
Nsyns = unpickledlist[1]
maxSynsPerSeg = unpickledlist[2]

unpicklefile = open('synlocs300.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
synlocsAll = unpickledlist[3]

gCoeffsAllAll = []

def nanint(x):
  if isnan(x):
    return nan
  return int(x)
def nanintlist(x,mylist):
  if isnan(x):
    return nan
  return mylist[int(x)]

maxIDtab = array([[126, 190, 294, 314, nan, 334, nan, nan, 370, nan, nan, nan, 422, 442, nan]])
IDtab = r_[(maxIDtab-1)/4] # map itercounters to counters, i.e., 1-4 -> 0, 5-8 -> 1, ..., 457-460 -> 114

unpicklefile = open('mutindexlist.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
mutinds = unpickledlist[:]
IDtab = [[nanintlist(x,mutinds) for x in y] for y in IDtab]
mutIDs = IDtab[0]
print "mutIDs="+str(mutIDs)

for icell in range(0,1):
  synlocs = synlocsAll[icell]
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]

  gCoeffsAll = []

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1
st1 = new IClamp(0.5)
L5PC.soma st1
forsec L5PC.somatic {
}
forsec L5PC.apical {
}
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
objref sl,st2,ns,syn1,syn2,con1,isyn, tvec, syns["""+str(2*Nsyns)+"""]
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
st2 = new IClamp(siteVec[1])
st2.amp = 0
L5PC.apic[siteVec[0]] {
  st2
  syn1 = new epsp(siteVec[1])
  syn1.tau0 = 0.5
  syn1.tau1 = 5
  syn1.onset = 145 + """+str(BACdt)+""" 
  syn1.imax = 0
  syn2 = new epsp(siteVec[1])
  syn2.tau0 = 0.5
  syn2.tau1 = 5
  syn2.onset = 145 + """+str(BACdt)+""" 
  syn2.imax = 0
  cvode.record(&syn1.i,isyn,tvec)
}
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
""")
  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocs[istim][0])+"""
siteVec[1] = """+str(synlocs[istim][1])+"""
access L5PC.apic[siteVec[0]]
L5PC.apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000 + """+str(BACdt)+"""
  syns["""+str(istim+Nsyns)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim+Nsyns)+"""].e = 0
  syns["""+str(istim+Nsyns)+"""].tau = 5
  syns["""+str(istim+Nsyns)+"""].onset = 10000 + """+str(BACdt)+"""
}
""")
  #extra"""

  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  #cols = ['#aaffaa','#aaffaa','#66ff66','#66aaaa','#00aaaa','#00aaaa']
  #cols = ['#00aaaa','#00bb77','#11cc44','#11dd11','#55ee00','#99dd00']
  cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
  col_control = '#2222ff'
  #yplus = [1, 2, 3, 4, 5, 6, -1, -2, -3]
  #yplus = [x+3 for x in yplus]
  yplus = [0, 0, 0, 0, 0, 0, 0, 0, 0]
  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

  close("all")
  f, axarr = plt.subplots(2, 2)
  iters = [0, 2, 5, 6, 8, -1]

  unpicklefile = open('thresholddistalamp300_cs'+str(icell)+'_comb.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  gsThisComb = unpickledlist[1]
  gCoeffsThisComb = []
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    if iter==5:
      continue
    defVals = mutation_stuff.getdefvals()
    keyList = defVals.keys()
    for idefval in range(0,len(keyList)):
      if type(defVals[keyList[idefval]]) is not list:
        defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
    mutText = ""
    gCoeffsThisIter = []
    gsThisIter = gsThisComb[iiter]
    for imutID in range(0,len(mutIDs)):
       if type(mutIDs[imutID]) is not list:
         continue
       igene = mutIDs[imutID][0]
       imut = mutIDs[imutID][1]
       iallmutval = mutIDs[imutID][2]

       nVals = len(MT[igene][imut])*[0]
       thesemutvars = []
       theseCoeffs = theseCoeffsAll[igene][imut]
       for imutvar in range(0,len(MT[igene][imut])):
         thesemutvars.append(MT[igene][imut][imutvar][0])
         if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
           MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
         nVals[imutvar] = len(MT[igene][imut][imutvar][1])
       cumprodnVals = cumprod(nVals)
       allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
       allmutvals = []
       for iallmutvaltmp in range(0,cumprodnVals[len(MT[igene][imut])-1]):
         allmutvals.append([0]*len(thesemutvars))
       for iallmutvaltmp in range(0,cumprodnVals[len(MT[igene][imut])-1]):
         for imutvar in range(0,len(MT[igene][imut])):
           if imutvar==0:
             allmutvals[iallmutvaltmp][imutvar] = MT[igene][imut][imutvar][1][iallmutvaltmp%nVals[imutvar]]
           else:
             allmutvals[iallmutvaltmp][imutvar] = MT[igene][imut][imutvar][1][(iallmutvaltmp/cumprodnVals[imutvar-1])%nVals[imutvar]]
    
       if iter >= 0:
         thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
       else:
         thisCoeff = 0
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
           defVals[mutvar] = newVal[:]

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
               print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
               h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
    print mutText
    thisCa = h.L5PC.soma[0].minCai_CaDynamics_E2
    times_all = []
    Vsoma_all = []

    f,axarr = subplots(len(currgs),len(PPIdts))
    for iPPI in range(0,len(PPIdts)):
         times_thisPPI = []
         Vsoma_thisPPI = []
         PPIdt = PPIdts[iPPI]

         for iterI in range(0,len(currgs)):
           for istim in range(0,Nsyns):
             h("syns["+str(istim)+"].gmax = "+str(gsThisIter*currCoeff))
             h("syns["+str(istim+Nsyns)+"].gmax = "+str(currgs[iterI]))
             h("syns["+str(istim+Nsyns)+"].onset = "+str(10000+BACdt+PPIdt))
           h("""
tstop = """+str(tstop)+"""
cai0_ca_ion = """+str(thisCa)+"""
v_init = """+str(v0)+"""
st1.amp = 0
st1.del = 200
st1.dur = 10
""")
           h.init()
           try:
             h.run()
           except RuntimeError:
             print "Too large I!"
             if iterI > 1:
               nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*(nextCoeffs[0]+nextCoeffs[2])]
             continue
           times=np.array(h.tvec)
           Casoma=np.array(h.casoma)
           Cadend=np.array(h.cadend)
           Vsoma=np.array(h.vsoma)
           Vdend=np.array(h.vdend)
           nSpikes_total = len(mytools.spike_times(times,Vsoma,-35,-37.5))
           axarr[iterI,iPPI].text(10050,-70,"nSp="+str(nSpikes_total),fontsize=4)
           mycol = col_control
           if iter >= 0:
             mycol = cols[iter]
           axarr[iterI,iPPI].plot(times,Vsoma,color=mycol)
           axarr[iterI,iPPI].set_xlim([10000,10000+PPIdt+100])
           if iPPI==0:
             axarr[iterI,0].set_ylabel("currg="+str(currgs[iterI]),fontsize=5)
           for tick in axarr[iterI,iPPI].xaxis.get_major_ticks()+axarr[iterI,iPPI].yaxis.get_major_ticks():
             tick.label.set_fontsize(3)
           times_thisPPI.append(times[:])
           Vsoma_thisPPI.append(Vsoma[:])
         times_all.append(times_thisPPI[:])
         Vsoma_all.append(Vsoma_thisPPI[:])
         axarr[0,iPPI].set_title('PPI='+str(PPIdts),fontsize=5)
         f.savefig("testppi300_fixed_comb_iiter"+str(iiter)+".eps")

    picklelist = [times_all,Vsoma_all]
    file = open('subthppitest_cs'+str(icell)+'_comb_iiter'+str(iiter)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()

    for idefval in range(0,len(defVals.keys())):
         thisdefval = defVals.keys()[idefval]
         if thisdefval.find('_Im') > -1:
           h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))

         else:
           h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))

    #Restore default values:
    for imutID in range(0,len(mutIDs)):
      if type(mutIDs[imutID]) is not list:
        continue
      igene = mutIDs[imutID][0]
      imut = mutIDs[imutID][1]
      iallmutval = mutIDs[imutID][2]

      nVals = len(MT[igene][imut])*[0]
      thesemutvars = []
      theseCoeffs = theseCoeffsAll[igene][imut]
      for imutvar in range(0,len(MT[igene][imut])):
        thesemutvars.append(MT[igene][imut][imutvar][0])
        if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
          MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
        nVals[imutvar] = len(MT[igene][imut][imutvar][1])
      cumprodnVals = cumprod(nVals)
      allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
      for imutvar in range(0,len(MT[igene][imut])):
        mutvars = allmutvars[iallmutval][imutvar]
        if type(mutvars) is str:
          mutvars = [mutvars]
        for kmutvar in range(0,len(mutvars)):
          mutvar = mutvars[kmutvar]
          newVal = defVals[mutvar]
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
              print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
              h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
