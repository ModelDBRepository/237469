import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
BACdt = 5.0
fs = 8
xs = range(700,1150,50);
tstop = 11000.0
currCoeff = 1.1 # use the threshold current I*1.1 for inducing the first spike, and I*1.1*2 for the second spike 
ITERS = 20
PPIdts = range(0,500,2)
maxLens = [1300,1185]

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

unpicklefile = open('thresholddistalamp300_control.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
gs_control = unpickledlist[0]
Nsyns = unpickledlist[1]
maxSynsPerSeg = unpickledlist[2]

unpicklefile = open('thresholddistalamp300_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
gsAllAll = unpickledlist[1]

DataAllAll = []
for icell in range(0,1):
  DataAll = []
  theseCoeffsAll = theseCoeffsAllAll[icell]

  counter = -1
  for igene in range(0,len(MT)):
   DataThisGene = []
   for imut in range(0,len(MT[igene])): 
    DataThisMut = []
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
      try:
        unpicklefile = open('ppispthrcoeff300_relthr_scaledonly_cs'+str(icell)+'_'+str(counter)+'.sav', 'r')
        unpickledlist = pickle.load(unpicklefile)
        unpicklefile.close()
        Data = unpickledlist[3]
      except:
        try:
          unpicklefile = open('ppispthrcoeff300_relthr_scaledonly_cs'+str(icell)+'_'+str(counter)+'_one.sav', 'r')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
          Data = unpickledlist[3]
        except:
          print 'ppispthrcoeff300_relthr_scaledonly_cs'+str(icell)+'_'+str(counter)+'_one.sav nor ppispthrcoeff300_relthr_scaledonly_cs'+str(icell)+'_'+str(counter)+'.sav found'
          Data = []
      DataThisMut.append(Data[:])
    DataThisGene.append(DataThisMut[:])
  
   DataAll.append(DataThisGene[:])

  DataAllAll.append(DataAll[:])
#          picklelist = 
#          file = open('ppispthrcoeff_relthr_cs'+str(icell)+'_control.sav', 'w')
#          pickle.dump(picklelist,file)
#          file.close()
picklelist = [theseCoeffsAllAll,gsAllAll,PPIdts,DataAllAll,MT]
file = open('ppispthrcoeff300_relthr_cs.sav', 'w')
pickle.dump(picklelist,file)
file.close()

