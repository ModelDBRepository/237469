#Copied from findthresholddistalamps300.py 22.11.2017
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
ITERS = 20
tstop = 11000.0

unpicklefile = open('thresholddistalamp300_control.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Nsyns = unpickledlist[1]
maxSynsPerSeg = unpickledlist[2]
maxLens = [1300,1185]

unpicklefile = open('synlocs300.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
synlocsAll = unpickledlist[3]

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

gsAllAll = []

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
  gsAll = []
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  theseCoeffsAll = theseCoeffsAllAll[icell]
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
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
objref sl,ns, tvec, syns["""+str(Nsyns)+"""]
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
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
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
}
""")
  styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
  coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

  gsThisComb = []
  for iter in [0, 2, 5, 6, 8, -1]:
    defVals = mutation_stuff.getdefvals()
    keyList = defVals.keys()
    for idefval in range(0,len(keyList)):
      if type(defVals[keyList[idefval]]) is not list:
        defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
    mutText = ""
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
       if iter == -1 and (igene > 0 or imut > 0 or iallmutval > 0):
         continue # do the control only once!
       print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)
          
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
    if icell==0:
      nextgs = [0.00,0.003,0.0015]
    if icell==1:
      nextgs = [0.00,0.06,0.03]
    hasSpiked = 0
    hasErred = 0
    for iterg in range(0,ITERS+2):
            thisg = nextgs[min(iterg,2)]
            h("st1.amp = "+str(thisg))
            for istim in range(0,Nsyns):
              h("syns["+str(istim)+"].gmax = "+str(thisg))
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
              hasErred = 1
              print "Too large g!"
              if iterg == 1:
                nextgs = [0.0,4.0,3.0]
                continue
              else:
                nextgs = [nextgs[0],nextgs[2],nextgs[0]+nextgs[2]]
              continue
  
            times=np.array(h.tvec)
            Vsoma=np.array(h.vsoma)
            spikes = mytools.spike_times(times,Vsoma,-35,-45)
            nSpikes1 = len(spikes)
            hasSpiked = hasSpiked or (nSpikes1 > 0)
  
            print "iterg="+str(iterg)+" done, g="+str(thisg)+", "+str(nSpikes1)+" spikes, iter="+str(iter)
            if iterg==0 and nSpikes1 > 0:
              print "Even zero g causes spiking!! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)+", iter="+str(iter)
              nextgs = [0.0,0.0,0.0]
              break
            if iterg==1 and not hasSpiked:
              print "No spiking with iterg==1, adding 25% to the current! igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
              nextgs = [nextgs[0],2.0*nextgs[1],1.25*nextgs[min(iterg,2)]]
              continue

            if iterg>=2 and iterg < ITERS+2:
              if nSpikes1 > 0:
                nextgs = [nextgs[0],nextgs[2],0.5*nextgs[0]+0.5*nextgs[2]]
              else:
                nextgs = [nextgs[2],nextgs[1],0.5*nextgs[1]+0.5*nextgs[2]]

    #Print the parameters and their default values:
    for idefval in range(0,len(defValsOrig.keys())):
         thisdefval = defValsOrig.keys()[idefval]
         if thisdefval.find('_Im') > -1:
           h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defValsOrig[thisdefval][1]))
         else:
           h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defValsOrig[thisdefval][0]))
  
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
          newVal = defValsOrig[mutvar]
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
    gsThisComb.append(nextgs[2])
    picklelist = [theseCoeffsAll,gsThisComb,MT]
    file = open('thresholddistalamp300_cs'+str(icell)+'_comb.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()
  
  gsAll.append(gsThisComb[:])
  picklelist = [theseCoeffsAll,gsThisComb,MT]
  file = open('thresholddistalamp300_cs'+str(icell)+'_comb.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()
