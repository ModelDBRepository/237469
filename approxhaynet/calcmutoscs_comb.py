import simosc_parallel_comb_varconn
from pylab import *
import pickle
from os.path import exists

Nmc = 150
gNoiseCoeff = 1.07
gSynCoeff = 1.07
seeds = range(1,1001)

oscamp = 0.25
oscfreqs = [0.5,0.625,0.75,0.875,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0,15.0]
oscphase = 0.0

mutcombID = int(sys.argv[4])
iosc = int(sys.argv[5])
iseed = int(sys.argv[6])

oscfreq = oscfreqs[iosc]
myrate = 1.0
myseed = seeds[iseed]

maxIDtab = array([[126, 190, 294, 314, nan, 334, nan, nan, 370, nan, nan, nan, 422, 442, nan]])
IDtab = r_[maxIDtab, maxIDtab-1, maxIDtab+1, maxIDtab+2] #epsilon=1/2, epsilon=1/4, epsilon=-1/4, epsilon=-1/2

if exists('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_combmutID'+str(mutcombID)+'_'+str(myrate)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav'):
  print 'spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_combmutID'+str(mutcombID)+'_'+str(myrate)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav exists, continuing'
else:
  Q=simosc_parallel_comb_varconn.simseedburst_func(Nmc,11000,IDtab[mutcombID],mutcombID,myseed,0.00039,0.0006,5,myrate,gNoiseCoeff,gSynCoeff,oscfreq,oscphase,oscamp,0,10.0,0)
