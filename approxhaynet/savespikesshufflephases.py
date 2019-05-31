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

rdSeed = 1
tstop = 12000 #ms
oscfreq = 1.0 #Hz
phase = 0 #in [0,2pi]

rateCoeffDiff = 0.25

for i in range(0,len(sys.argv)):
  print "sys.argv["+str(i)+"]: "+sys.argv[i]

if len(sys.argv) > 1:
  rdSeed = int(sys.argv[1])
if len(sys.argv) > 2:
  tstop = int(sys.argv[2])
if len(sys.argv) > 3:
  oscfreq = float(sys.argv[3])
if len(sys.argv) > 4:
  rateCoeffDiff = float(sys.argv[4])

rateCoeffs = [1-rateCoeffDiff, 1+rateCoeffDiff]

phase = 2.0*numpy.pi*rdSeed/50

numpy.random.seed(rdSeed)
sps = []
for i in range(0,Nmc):
  spE=mytools.oscillatorypoissontimeseries(NsynE,rateE*rateCoeffs[0],rateE*rateCoeffs[1],oscfreq,phase,tstop,0)
  spI=mytools.oscillatorypoissontimeseries(NsynI,rateI*rateCoeffs[0],rateI*rateCoeffs[1],oscfreq,phase,tstop,0)
  print "cell "+str(i)+" done"
  sps.append([spE[:],spI[:]])

file = open('backgroundspikes_'+str(oscfreq)+'_'+str(rateCoeffDiff)+'_'+str(tstop)+'_'+str(rdSeed)+'.sav', 'w')
pickle.dump(sps,file)
file.close()

