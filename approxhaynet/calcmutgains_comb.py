#Copied from calcspikes_conn_more_more.py 10.5.2016
import mytools
import simseedburst_func_comb_varconn
from pylab import *
from neuron import h
import pickle
from os.path import exists
import time

Nmc = 150
rates = [0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]
#gSynCoeffs = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]
seeds = range(1,1001)
cols = ['#0000FF','#400090','#900040','#FF0000','#904000','#409000','#00FF00','#30A030','#606060']

extensions_final = ['_final', '_extended_final', '_withisi_final', '_withisi_extended_final', '_final_alternative', '_extended_final_alternative']
mutcombID = int(sys.argv[4])
igsyn = 0
irate = int(sys.argv[5])
iseed = int(sys.argv[6])

gSynCoeff = 1.07
gNoiseCoeff = 1.07
rateCoeff = rates[irate]
myseed = seeds[iseed]

maxIDtab = array([[126, 190, 294, 314, nan, 334, nan, nan, 370, nan, nan, nan, 422, 442, nan]])
IDtab = r_[maxIDtab, maxIDtab-1, maxIDtab+1, maxIDtab+2] #epsilon=1/2, epsilon=1/4, epsilon=-1/4, epsilon=-1/2


if exists('spikes_parallel'+str(Nmc)+'_mutcombID'+str(mutcombID)+'_'+str(rateCoeff)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav'):
  print 'spikes_parallel'+str(Nmc)+'_mutcombID'+str(mutcombID)+'_'+str(rateCoeff)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav exists'
else:
  Q=simseedburst_func_comb_varconn.simseedburst_func(Nmc,11000,IDtab[mutcombID],mutcombID,myseed,0.00039,0.0006,5,rateCoeff,gNoiseCoeff,gSynCoeff,1)
