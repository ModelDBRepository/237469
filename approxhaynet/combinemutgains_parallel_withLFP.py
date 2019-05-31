import pickle
import sys
import numpy

print "importings OK"

Nmc = 150
gNoiseCoeff = 1.07
seeds = range(1,1001)
cols = ['#0000FF','#400090','#900040','#FF0000','#904000','#409000','#00FF00','#30A030','#606060']
rates = [0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]
seeds = range(1,1001)

gSynCoeff = 1.07
gNoiseCoeff = 1.07
imut = int(sys.argv[1])
irate = int(sys.argv[2])
iseed = int(sys.argv[3])
myseed = seeds[iseed]

if len(sys.argv) > 4:
  nCPUs = int(sys.argv[4])
else:
  nCPUs = Nmc

nseg = 5
tstop = 11000
mutID = imut
Econ=0.00039
Icon=0.0006
rateCoeff = rates[irate]

spikes = []
spikedCells = []
dipolesE = numpy.zeros([440000,3])
for i in range(0,nCPUs):
  unpicklefile = open('spikes_parallel_'+str(nseg)+'_'+str(tstop)+'_'+str(mutID)+'_'+str(Econ)+'_'+str(Icon)+'_'+str(rateCoeff)+'_'+str(gNoiseCoeff)+'_'+str(gSynCoeff)+'_'+str(myseed)+'_'+str(i)+'_of_'+str(Nmc)+'_withLFP.sav','r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  spikes.append(unpickledlist[0])
  spikedCells.append(unpickledlist[1])
  dipolesE = dipolesE + unpickledlist[3]
  print "load "+str(i)+ " OK"

picklelist = [spikes,spikedCells,dipolesE]
file = open('spikes_parallel'+str(Nmc)+'_mutID'+str(mutID)+'_'+str(rateCoeff)+'_gNoise'+str(gNoiseCoeff)+'_gsyn'+str(gSynCoeff)+'_seed'+str(myseed)+'_withLFP.sav', 'w')
pickle.dump(picklelist,file)
file.close()


