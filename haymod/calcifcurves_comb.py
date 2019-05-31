#TODO: change counterID to triple (gene,mut,mutval) ID

from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import scipy.io

morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
Is = [0.35+0.05*x for x in range(0,22)]
fs = 8
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

import mutation_stuff
MT = mutation_stuff.getMT()
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
theseMutValsAll = unpickledlist[2]

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

spTimesAllAll = []

for icell in range(0,1):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  spTimesAll = []
  spTimesAll2 = []
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

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
objref st1, tvec, vsoma, casoma
L5PC.soma st1 = new IClamp(0.5)
tvec = new Vector()
vsoma = new Vector()
casoma = new Vector()
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
""")

  counter = -1
  for imutcomb in range(0,len(IDtab)):
    muts = IDtab[imutcomb]

    spTimesThisVal = []
    spTimesThisVal2 = []
    iters = [0, 2, 6, 8]
    for iiter in range(0,len(iters)):
      iter = iters[iiter]
      defVals = mutation_stuff.getdefvals() # reload this and update after each application of variant
      keyList = defVals.keys()
      for idefval in range(0,len(keyList)):
        if type(defVals[keyList[idefval]]) is not list:
          defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]

      for jmut in range(0,len(muts)):
        if type(muts[jmut]) is not list:
          continue
        igene = muts[jmut][0]
        imut = muts[jmut][1]
        print "imutcomb="+str(imutcomb)+", iiter="+str(iiter)+". Adding variant igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(muts[jmut][2])+"..."
        if imut >= len(MT[igene]) or imut < 0:
          print "                   something's wrong"
          continue

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
        for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
          allmutvals.append([0]*len(thesemutvars))
        for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
          for imutvar in range(0,len(MT[igene][imut])):
            if imutvar==0:
              allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
            else:
              allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]
      
        iallmutval = muts[jmut][2]

        if iter >= 0:
          thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][1]*(1.0 - 0.5*theseCoeffs[iallmutval])
        else:
          thisCoeff = 0

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
              newVal = [x+mutvals*thisCoeff for x in defVals[mutvar]]
              defVals[mutvar] = [x+mutvals*thisCoeff for x in defVals[mutvar]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              defVals[mutvar] = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
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
                print """forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }"""
                h("""forsec L5PC."""+str(updatedVars[iupdated])+""" {
  """+mutvar+""" = """+str(newVal[whichDefVal[iupdated]])+"""
  }""")
      print mutText

      spTimesThisCoeff = []
      spTimesThisCoeff2 = []
      for iI in range(0,len(Is)):
        tstop = 8000.0
        squareAmp = Is[iI]
        squareDur = 7800.0
        epsp_Imax = 0.0
        h("""
  tstop = """+str(tstop)+"""
  v_init = """+str(v0)+"""
  cai0_ca_ion = """+str(ca0)+"""
  st1.amp = """+str(squareAmp)+"""
  st1.del = 200
  st1.dur = """+str(squareDur)+"""
  """)
        h.init()
        h.run()
  
        times=np.array(h.tvec)
        Vsoma=np.array(h.vsoma)
        spikes = mytools.spike_times(times,Vsoma,-35,-45)
        spikes2 = mytools.spike_times(times,Vsoma,-35,100)

        spTimesThisCoeff.append(spikes[:])
        spTimesThisCoeff2.append(spikes2[:])


      if iter==-1:
        picklelist = [theseCoeffsAllAll,spTimesThisCoeff,spTimesThisCoeff2,MT]
        file = open('ifcurve_comb_cs'+str(icell)+'_control.sav', 'w')
        pickle.dump(picklelist,file)
        file.close()
      spTimesThisVal.append(spTimesThisCoeff[:])
      spTimesThisVal2.append(spTimesThisCoeff2[:])
  

    picklelist = [theseCoeffsAllAll,spTimesThisVal,spTimesThisVal2,IDtab,MT]
    file = open('ifcurves_comb_cs'+str(icell)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()
   
    defVals = mutation_stuff.getdefvals() # restore the default values                                                                                            
    keyList = defVals.keys()
    for idefval in range(0,len(keyList)):
      if type(defVals[keyList[idefval]]) is not list:
        defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]       
    #Print the parameters and their default values:
    for idefval in range(0,len(defVals.keys())):
      thisdefval = defVals.keys()[idefval]
      if thisdefval.find('_Im') > -1:
        h('print "L5PC.apic[0].'+thisdefval+' = ", L5PC.apic[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][1]))
      else:
        h('print "L5PC.soma[0].'+thisdefval+' = ", L5PC.soma[0].'+thisdefval+', "Default = ", '+str(defVals[thisdefval][0]))

    #Restore default values:
    for imutvar in range(0,len(MT[igene][imut])):
        mutvars = allmutvars[iallmutval][imutvar]
        mutvals = allmutvals[iallmutval][imutvar]
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
  

  
