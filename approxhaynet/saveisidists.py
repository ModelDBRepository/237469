import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import time
import scipy.io
import pickle
import sys
import mutation_stuff
import mytools

Nmc = 150
rateE = 0.72
rateI = 7.0
NsynE = 10000
NsynI = 2500

rdSeeds = range(1,51)
tstop = 12000 #ms
oscfreq = 1.0 #Hz

rateCoeffDiff = 0.25

for i in range(0,len(sys.argv)):
  print "sys.argv["+str(i)+"]: "+sys.argv[i]
if len(sys.argv) > 1:
  oscfreq = float(sys.argv[1])
if len(sys.argv) > 2:
  rateCoeffDiff = float(sys.argv[2])

rateCoeffs = [1-rateCoeffDiff, 1+rateCoeffDiff]

phases = [2.0*numpy.pi*x/50 for x in rdSeeds]


tmax = 10000
nts = 100000
nphases = 30

probsE = numpy.zeros([nphases,nts])
probsI = numpy.zeros([nphases,nts])

for irdSeed in range(0,len(rdSeeds)):
  try:
    rdSeed = rdSeeds[irdSeed]
    phase = phases[irdSeed]
    unpicklefile = open('backgroundspikes_'+str(oscfreq)+'_'+str(rateCoeffDiff)+'_'+str(tstop)+'_'+str(rdSeed)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    for icell in range(0,len(unpickledlist)):
      theseSpikes = unpickledlist[icell][0]
      iprev = -1
      tprev = 0
      for ispike in range(0,len(theseSpikes)):
        t = theseSpikes[ispike][1]
        if theseSpikes[ispike][0] != iprev:
          tprev = 0
        iprev = theseSpikes[ispike][0]
        dt = t - tprev        
        myphase = phase + 2*numpy.pi*oscfreq*0.001*tprev
        while myphase >= 2*numpy.pi:
          myphase = myphase-2*numpy.pi
        phaseslot = int(myphase/(2*numpy.pi)*nphases)
        dtslot = int(dt/tmax*nts)
        if dtslot >= nts:
          dtslot = nts-1
        tprev = t
        probsE[phaseslot,dtslot] = probsE[phaseslot,dtslot] + 1

      theseSpikes = unpickledlist[icell][1]
      iprev = 0
      tprev = 0
      for ispike in range(0,len(theseSpikes)):
        t = theseSpikes[ispike][1]
        if theseSpikes[ispike][0] != iprev:
          tprev = 0
        iprev = theseSpikes[ispike][0]
        dt = t - tprev        
        myphase = phase + 2*numpy.pi*oscfreq*0.001*tprev
        while myphase >= 2*numpy.pi:
          myphase = myphase-2*numpy.pi
        phaseslot = int(myphase/(2*numpy.pi)*nphases)
        dtslot = int(dt/tmax*nts)
        if dtslot >= nts:
          dtslot = nts-1
        tprev = t
        probsI[phaseslot,dtslot] = probsI[phaseslot,dtslot] + 1
  except:
    print 'Could not open backgroundspikes_'+str(oscfreq)+'_'+str(rateCoeffDiff)+'_'+str(tstop)+'_'+str(rdSeed)+'.sav'
for iphase in range(0,nphases):
  probsE[iphase,:] = probsE[iphase,:]/sum(probsE[iphase,:])
  probsI[iphase,:] = probsI[iphase,:]/sum(probsI[iphase,:])

file = open('isidists_'+str(oscfreq)+'_'+str(rateCoeffDiff)+'.sav', 'w')
pickle.dump([probsE,probsI,tmax],file)
file.close()

