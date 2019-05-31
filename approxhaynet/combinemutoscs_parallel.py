#Combine spikes_*_of_* files (where population spike trains of each MPI are saved to different files) into one spikes_ file.
#run as python combinemutgains_parallel.py $imut $iosc $iseed $NUMP

#import simseedburst_func
#import mytools
#from pylab import *
#from neuron import h
import pickle
import sys

print "importings OK"

Nmc = 150
seeds = range(1,1001)
cols = ['#0000FF','#400090','#900040','#FF0000','#904000','#409000','#00FF00','#30A030','#606060']
seeds = range(1,1001)

oscamp = 0.25
oscfreqs = [0.5,0.625,0.75,0.875,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0,15.0]
oscphase = 0.0

gSynCoeff = 1.07
gNoiseCoeff = 1.07
mutID = int(sys.argv[1])
iosc = int(sys.argv[2])
iseed = int(sys.argv[3])

oscfreq = oscfreqs[iosc]
myseed = seeds[iseed]



if len(sys.argv) > 4:
  nCPUs = int(sys.argv[4])
else:
  nCPUs = Nmc

nseg = 5
tstop = 11000
Econ=0.00039
Icon=0.0006
rateCoeff = 1.0

spikes = []
spikedCells = []
for i in range(0,nCPUs):
  unpicklefile = open('spikes_parallel_osc'+str(oscfreq)+'_'+str(nseg)+'_'+str(tstop)+'_'+str(mutID)+'_'+str(Econ)+'_'+str(Icon)+'_'+str(rateCoeff)+'_'+str(gNoiseCoeff)+'_'+str(gSynCoeff)+'_'+str(myseed)+'_'+str(i)+'_of_'+str(Nmc)+'.sav','r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  spikes.append(unpickledlist[0])
  spikedCells.append(unpickledlist[1])
  print "load "+str(i)+ " OK"

picklelist = [spikes,spikedCells]
file = open('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_mutID'+str(mutID)+'_'+str(rateCoeff)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'.sav', 'w')
pickle.dump(picklelist,file)
file.close()


