#Copied from calcmutgains 18.10.2016
import mytools
import simosc_parallel
from pylab import *
from neuron import h
import pickle
from os.path import exists
import time

Nmc = 150
gNoiseCoeff = 1.07
gSynCoeff = 1.07
seeds = range(1,1001)
cols = ['#0000FF','#400090','#900040','#FF0000','#904000','#409000','#00FF00','#30A030','#606060']

oscamp = 0.25
oscfreqs = [0.5,0.625,0.75,0.875,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0,15.0]
oscphase = 0.0

igsyn = 0
mutID = int(sys.argv[4])
iosc = int(sys.argv[5])
iseed = int(sys.argv[6])

oscfreq = oscfreqs[iosc]
myrate = 1.0
myseed = seeds[iseed]
if exists('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_mutID'+str(mutID)+'_'+str(myrate)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav'):
  print 'spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_mutID'+str(mutID)+'_'+str(myrate)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav exists, continuing...'
else:
  Q=simosc_parallel.simseedburst_func(Nmc,11000,mutID,myseed,0.00039,0.0006,5,myrate,gNoiseCoeff,gSynCoeff,oscfreq,oscphase,oscamp,0,10.0,0)
