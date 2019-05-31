import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys

Is = [0.2,0.4,0.6,0.8,1.0,1.2,1.4]

import mutation_stuff
MT = mutation_stuff.getMT()

theseCoeffsAllAll = []
for icell in range(0,2):
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
     try:
       unpicklefile = open('scalings_cs'+str(icell)+'_'+str(counter)+'.sav', 'r')
       unpickledlist = pickle.load(unpicklefile)
       unpicklefile.close()
       theseCoeffsMut.append(unpickledlist[0])
     except:
       theseCoeffsMut.append([])
   theseCoeffsGene.append(theseCoeffsMut[:])
  theseCoeffsAll.append(theseCoeffsGene[:])
 theseCoeffsAllAll.append(theseCoeffsAll[:])
 print theseCoeffsAll
 
picklelist = [theseCoeffsAllAll,theseMutVarsAll,theseMutValsAll,MT]
file = open('scalings_cs.sav', 'w')
pickle.dump(picklelist,file)
file.close()


