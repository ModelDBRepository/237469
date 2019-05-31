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

if len(sys.argv) > 1:
  upfromd = int(float(sys.argv[1]))
else:
  upfromd = 300

unpicklefile = open('thresholddistalamp'+str(upfromd)+'_control.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
gmaxes = unpickledlist[0]
Nsyns = unpickledlist[1]
maxSynsPerSeg = unpickledlist[2]
maxLens = [1300,1185]

unpicklefile = open('synlocs'+str(upfromd)+'.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
synlocsAll = unpickledlist[3]

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
theseMutValsAllAll = unpickledlist[2]

gsAllAll = []

for icell in range(0,1):
  synlocs = synlocsAll[icell]
  gsAll = []
  theseCoeffsAll = theseCoeffsAllAll[icell]

  counter = -1
  
  for igene in range(0,len(MT)):
   gsThisGene = []
   for imut in range(0,len(MT[igene])):
    gsThisMut = []
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
      
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      counter = counter + 1                                                                                                                                                               
      unpicklefile = open('thresholddistalamp'+str(upfromd)+'_cs'+str(icell)+'_'+str(counter)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      gsThisMut.append(unpickledlist[1])
    gsThisGene.append(gsThisMut[:])
   gsAll.append(gsThisGene[:])
  gsAllAll.append(gsAll[:])
  
picklelist = [theseCoeffsAllAll,gsAllAll,MT]
file = open('thresholddistalamp'+str(upfromd)+'_cs.sav', 'w')
pickle.dump(picklelist,file)
file.close()
  
