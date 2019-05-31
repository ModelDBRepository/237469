import mytools
import simseedburst_func
from pylab import *
from neuron import h
import pickle
from os.path import exists
import time

Nmc = 150
rates = [0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]
seeds = range(1,1001)
cols = ['#0000FF','#400090','#900040','#FF0000','#904000','#409000','#00FF00','#30A030','#606060']

extensions_final = ['_final', '_extended_final', '_withisi_final', '_withisi_extended_final', '_final_alternative', '_extended_final_alternative']
igsyn = 0
mutID = int(sys.argv[4])
irate = int(sys.argv[5])
iseed = int(sys.argv[6])

rateCoeff = rates[irate]
gSynCoeff = 1.07
gNoiseCoeff = 1.07
myseed = seeds[iseed]
if exists('spikes_parallel'+str(Nmc)+'_mutID'+str(mutID)+'_'+str(rateCoeff)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav'):
  print 'spikes_parallel'+str(Nmc)+'_mutID'+str(mutID)+'_'+str(rateCoeff)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav exists'
else:
  Q=simseedburst_func.simseedburst_func(Nmc,11000,mutID,myseed,0.00039,0.0006,5,rateCoeff,gNoiseCoeff,gSynCoeff,0)



