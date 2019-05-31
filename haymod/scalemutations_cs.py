from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
import os.path

import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments

unpicklefile = open('control_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spikfreqs_control_All = unpickledlist[0]
timesc_control_All = unpickledlist[1]
Vsomac_control_All = unpickledlist[2]
VDerivc_control_All = unpickledlist[3]
VDcoeff_control_All = unpickledlist[4]
Is_control = unpickledlist[19]
Is = [0.2,0.4,0.6,0.8,1.0,1.2,1.4]

theseCoeffsAllAll = []

squareAmps = [[0.696,0,1.137],[0.872,0,0.993]]
epsp_gmaxs = [[0,0.0612,0.100],[0,0.455,0.518]]

for icell in range(0,2):
  spikfreqs_control = mytools.interpolate(Is_control,spikfreqs_control_All[icell],Is)
  Vsomac_control = Vsomac_control_All[icell]
  VDerivc_control = VDerivc_control_All[icell]
  VDcoeff_control = VDcoeff_control_All[icell]
  timesc_control = timesc_control_All[icell]
  print spikfreqs_control

  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  v0 = -80
  ca0 = 0.0001
  proximalpoint = 400
  distalpoint = 620
  BACdt = 5.0

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
objref sl,st2,ns,syn1,con1,isyn, tvec
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
st2 = new IClamp(siteVec[1])
st2.amp = 0
L5PC.apic[siteVec[0]] {
  st2
  syn1 = new AlphaSynapse(siteVec[1])
  syn1.e = 0
  syn1.tau = 5
  syn1.onset = 200 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
cadend2 = new Vector()
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
sl = new List()
sl = L5PC.locateSites("apic","""+str(proximalpoint)+""")
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
access L5PC.apic[siteVec[0]]
recSite = new IClamp(siteVec[1])
recSite.amp = 0
L5PC.apic[siteVec[0]] {
        recSite
}
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend2,tvec)
cvode.record(&cai(siteVec[1]),cadend2,tvec)
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  ITERS = 20
  theseCoeffsAll = []
  theseMutValsAll = []
  theseMutVarsAll = []

  counter = -1
  for igene in range(0,len(MT)):
   theseCoeffsGene = []
   for imut in range(0,len(MT[igene])):
    theseCoeffsMut = []
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    for imutvar in range(0,len(MT[igene][imut])):
      thesemutvars.append(MT[igene][imut][imutvar][0])
      if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
        MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
      nVals[imutvar] = len(MT[igene][imut][imutvar][1])
    cumprodnVals = cumprod(nVals)
    allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars[:]]
    allmutvals = []
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      allmutvals.append([0]*len(thesemutvars))
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar==0:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
        else:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
    theseMutValsAll.append(allmutvals[:])  
    theseMutVarsAll.append(allmutvars[:])  
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      counter = counter + 1
      if len(sys.argv) > 1 and int(float(sys.argv[1])) != counter:
        continue
      nextCoeffs = [0.0,2.0,1.0]
      for iter in range(0,ITERS+2+3):
        thisCoeff = nextCoeffs[min(iter,2)]
   
        mutText = ""
        for imutvar in range(0,len(MT[igene][imut])):
          if imutvar > 0 and imutvar%2==0:
            mutText = mutText+"\n"
          mutvars = allmutvars[iallmutval][imutvar]
          if type(mutvars) is str:
            mutvars = [mutvars]
          mutText = mutText + str(mutvars) + ": "
          mutvals = allmutvals[iallmutval][imutvar]
          for kmutvar in range(0,len(mutvars)):
            if mutvars[kmutvar].find('offm') > -1 or mutvars[kmutvar].find('offh') > -1 or mutvars[kmutvar].find('ehcn') > -1:
              newVal = [x+thisCoeff*mutvals for x in defVals[mutvars[kmutvar]]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals*thisCoeff) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals*thisCoeff) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvars[kmutvar]]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals**thisCoeff)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
            if mutvars[kmutvar].find('_Ih') > -1:
              updateThese = [1,1,1]
            elif mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_Ca_LVAst') > -1 or mutvars[kmutvar].find('_SKv3.1') > -1 or mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_SK_E2') > -1 or mutvars[kmutvar].find('_NaTa_t') > -1 or mutvars[kmutvar].find('_CaDynamics_E2') > -1:
              updateThese = [1,1,0]
            elif mutvars[kmutvar].find('_K_Pst') > -1 or mutvars[kmutvar].find('_K_Tst') > -1 or mutvars[kmutvar].find('_Nap_Et2') > -1: 
              updateThese = [1,0,0]
            elif mutvars[kmutvar].find('_Im') > -1:
              updateThese = [0,1,0]
            else:
              print "Error: str=" + str(mutvars[kmutvar])
              updatedVars = [0,0,0]
            for iupdated in range(0,3):
              if updateThese[iupdated]:
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
        print mutText
    
        ############################################# Condition 1: Short burst #############################################
        tstop = 500.0
        squareAmp = squareAmps[icell][0]
        squareDur = 150.0
        epsp_gmax = 0.0
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
        h.init()
        h.run()

        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes1 = len(spikes)

        close("all")

        f, axarr = plt.subplots(2, 3)
        axarr[0,0].plot(times, Vsoma)
        axarr[0,0].set_title("Perisomatic firing, nspikes="+str(nSpikes1))
        axarr[0,0].set_ylim([-100,40])
        for ix in range(0,3):
          for iy in range(0,2):
            axarr[iy,ix].set_position([0.05+0.3*ix, 0.05+0.4*(1-iy), 0.23, 0.3])

        ############################################# Condition 2: Distal EPSC #############################################
        squareAmp = 0.0
        epsp_gmax = epsp_gmaxs[icell][1]
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
        h.init()
        h.run()
        times=np.array(h.tvec)
        Vdend=np.array(h.vdend)
        Vdend2=np.array(h.vdend2)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes2 = len(spikes)
        axarr[1,0].plot(times, Vsoma)
        axarr[1,0].plot(times, Vdend2)
        axarr[1,0].plot(times, Vdend)
        axarr[1,0].set_ylim([-100,40])
        axarr[1,0].set_title("Strong distal EPSC, nspikes="+str(nSpikes2))

        ############################################# Condition 3: Somatic stim + EPSC #############################################
        squareAmp = squareAmps[icell][2]
        squareDur = 10
        epsp_gmax = epsp_gmaxs[icell][2]
        h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
        h.init()
        h.run()
        times=np.array(h.tvec)
        Vdend=np.array(h.vdend)
        Vdend2=np.array(h.vdend2)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        nSpikes3 = len(spikes)
        axarr[1,1].plot(times, Vsoma)
        axarr[1,1].plot(times, Vdend2)
        axarr[1,1].plot(times, Vdend)
        axarr[1,1].set_ylim([-100,40])
        axarr[1,1].set_title("Short stimulus at soma + weak EPSC, nspikes="+str(nSpikes3))

        ############################################# Condition 4: IF curve #############################################
        spikfreqs = len(Is)*[0]
        for iI in range(0,len(Is)):
          tstop = 4000.0
          squareAmp = Is[iI]
          squareDur = 3800.0
          epsp_gmax = 0.0
          h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = """+str(squareAmp)+"""
st1.del = 200
st1.dur = """+str(squareDur)+"""
syn1.gmax = """+str(epsp_gmax)+"""
syn1.onset = 200 + """+str(BACdt)+""" 
""")
          h.init()
          h.run()

          times=np.array(h.tvec)
          Vsoma=np.array(h.vsoma)
          spikes = mytools.spike_times(times,Vsoma,-35,100)
          spikfreqs[iI] = sum([1 for x in spikes if x >= 500.0])/3.5
          if iI==4: # use the memb. pot. time course of 1.0nA for the limit cycle
            times_lc = times[:]
            Vsoma_lc = Vsoma[:]
            spikes_lc = spikes[:]

        axarr[0,2].plot(Is, spikfreqs)
        axarr[0,2].set_xlim([0,1.25])
        axarr[0,2].set_ylim([0,20])
        spikfreqdiffsum = sum([abs(x-y) for x,y in zip(spikfreqs,spikfreqs_control)])
        spikfreqdiffrel = spikfreqdiffsum/sum(spikfreqs_control)
        axarr[0,2].set_title("IF curve, diff="+str(spikfreqdiffrel))
        print "IF curve, diff="+str(spikfreqdiffrel)

        ############################################# Condition 5: Limit cycle #############################################
        if len(spikes_lc) < 3:
          lcdiff = 1e6
        else:
          spts = spikes_lc[len(spikes_lc)-3:len(spikes_lc)]
          istart = next((i for i,x in enumerate(times_lc) if x > spts[0]))
          iend = next((i for i,x in enumerate(times_lc) if x > spts[1]))+4
          nsteps = iend-istart-1
          Vsomac = Vsoma_lc[istart:iend]
          timesc = times_lc[istart:iend]
          VDerivc = mytools.membpotderivs(timesc,Vsomac)
          VDcoeff =  mytools.limitcyclescaledv(Vsomac,VDerivc,Vsomac,VDerivc)
          lcdiff1 = mytools.limitcyclediff(Vsomac[1:nsteps-1],VDerivc,Vsomac_control,VDerivc_control,VDcoeff_control)
          lcdiff2 = mytools.limitcyclediff(Vsomac_control,VDerivc_control,Vsomac[1:nsteps-1],VDerivc,VDcoeff_control)
          lcdiff = 0.5*(lcdiff1+lcdiff2)
          axarr[1,2].plot(Vsomac[1:nsteps-1],VDerivc)
          axarr[1,2].set_title("Limit cycle")
          axarr[1,2].set_xlim([-70,30])
          axarr[1,2].set_ylim([-200,600])

        f.suptitle(mutText)
        if iter < ITERS+2:
          f.savefig("vrecs_cs"+str(icell)+"_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_ITER"+str(iter)+".png")
        else:
          f.savefig("vrecs_cs"+str(icell)+"_mut"+str(igene)+"_"+str(imut)+"_"+str(iallmutval)+"_TEST"+str(iter-ITERS-2)+".png")

        #Print the parameters and their default values:
        for idefval in range(0,len(defVals.keys())):
          thisdefval = defVals.keys()[idefval]
          if thisdefval.find('_Im') > -1:
            h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
            #) #+" (def="+str(defVals[thisdefval])+")"
          else:
            h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))
            #h('print L5PC.soma[0]."+thisdefval) #+" (def="+str(defVals[thisdefval])+")"

        isChanged = nSpikes1 != 4 or nSpikes2 != 1 or nSpikes3 != 2 or spikfreqdiffrel > 0.15 or lcdiff > 600.0
        print isChanged
        if iter==0 and isChanged:
          print "Even null mutation causes different spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          continue
        if iter==1 and not isChanged:
          print "This mutation effect does not alter spiking even when doubled!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
          continue
        if iter>=2 and iter < ITERS+2:
          if isChanged:
            nextCoeffs = [nextCoeffs[0],nextCoeffs[2],0.5*nextCoeffs[0]+0.5*nextCoeffs[2]]
          else:
            nextCoeffs = [nextCoeffs[2],nextCoeffs[1],0.5*nextCoeffs[1]+0.5*nextCoeffs[2]]
        if iter == ITERS+1:
          nextCoeffs = [nextCoeffs[2],nextCoeffs[2],nextCoeffs[2]*0.99]
        if iter == ITERS+2:
          nextCoeffs = [nextCoeffs[0],nextCoeffs[0],nextCoeffs[0]*1.0]
        if iter == ITERS+3:
          nextCoeffs = [nextCoeffs[0],nextCoeffs[0],nextCoeffs[0]*1.01]
      

      #Restore default values:
      for imutvar in range(0,len(MT[igene][imut])):
        mutvars = allmutvars[iallmutval][imutvar]
        if type(mutvars) is str:
          mutvars = [mutvars]
        mutvals = allmutvals[iallmutval][imutvar]
        for kmutvar in range(0,len(mutvars)):
          newVal = defVals[mutvars[kmutvar]]
          if mutvars[kmutvar].find('_Ih') > -1:
            updateThese = [1,1,1]
          elif mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_Ca_LVAst') > -1 or mutvars[kmutvar].find('_SKv3.1') > -1 or mutvars[kmutvar].find('_Ca_HVA') > -1 or mutvars[kmutvar].find('_SK_E2') > -1 or mutvars[kmutvar].find('_NaTa_t') > -1 or mutvars[kmutvar].find('_CaDynamics_E2') > -1:
            updateThese = [1,1,0]
          elif mutvars[kmutvar].find('_K_Pst') > -1 or mutvars[kmutvar].find('_K_Tst') > -1 or mutvars[kmutvar].find('_Nap_Et2') > -1: 
            updateThese = [1,0,0]
          elif mutvars[kmutvar].find('_Im') > -1:
            updateThese = [0,1,0]
          else:
            print "Error: str=" + str(mutvars[kmutvar])
            updatedVars = [0,0,0]
          for iupdated in range(0,3):
            if updateThese[iupdated]:
              h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
"""+mutvars[kmutvar]+""" = """+str(newVal[whichDefVal[iupdated]])+"""
}""")
      theseCoeffsMut.append(nextCoeffs[0]+0.0)
      picklelist = [nextCoeffs[0]+0.0,igene,imut,iallmutval,counter,MT]
      file = open('scalings_cs'+str(icell)+'_'+str(counter)+'.sav', 'w')
      pickle.dump(picklelist,file)
      file.close()

    theseCoeffsGene.append(theseCoeffsMut[:])
   theseCoeffsAll.append(theseCoeffsGene[:])
  theseCoeffsAllAll.append(theseCoeffsAll[:])

