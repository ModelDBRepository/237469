import matplotlib
matplotlib.use('Agg')
import mytools
from pylab import *
from neuron import h
import pickle
from os.path import exists
import numpy
import scipy.io
import time

Nmc = 150
rates = [1.0]
seeds = range(1,1000)
cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
col_control = '#2222ff'
oscamp = 0.25
oscfreqs = [0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5, 10.0, 15.0]
oscphase = 0.0

icell = 0
iext = 3

gsyn = 1.07
gNoise = 1.07
myrate = rates[0]

coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
dt = 0.025

dodraw = False
if len(sys.argv) > 1:
  dodraw = int(sys.argv[1]) > 0

import mutation_stuff
MT = mutation_stuff.getMT()
geneNames = mutation_stuff.getgenenames()
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
theseMutValsAll = unpickledlist[2]

Nspikes_samples = []
Nspikes_before1000_samples = []
Nspikes_after1000_samples = []
Nspikes500_11000_samples = []
gSynCoeffs_samples = []
spikets_samples = []
spikers_samples = []
rates_samples = []
oscfreqs_samples = []

ft_df = 0.00001      # kHz
ft_maxf = 0.05       # kHz
ft_fs = [ft_df*x for x in range(0,int(round(ft_maxf/ft_df))+1)]

counter = -1
icell = 0
theseCoeffsAll = theseCoeffsAllAll[icell]

f,axarr = subplots(3,1)
axnew = []
for i in range(0,3):
  axarr[i].set_position([0.12,0.7-0.3*i,0.43,0.28])
  axnew.append(f.add_axes([0.34, 0.7-0.3*i+0.18, 0.18, 0.09],axisbg='w'))

unpicklefile = open('spikes_parallel_5_11000_0_0.00039_0.0006_1.0_1.07_1.07_1_0_of_1_withLFP.sav', 'r') #from runsinglecellLFP.py
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spike_ts = unpickledlist[0]
Vs = unpickledlist[2][0]
dipoles = [x[1] for x in unpickledlist[3]]

unpicklefile = open('spikes_parallel_5_11000_0_0.00039_0.0006_1.0_1.07_1.07_1_0_of_1_withLFP_EEG.sav', 'r') #from runsinglecellLFP.py and calcEEG_uncombined.py
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
EEGs = unpickledlist[0][0]

axarr[0].plot([0.025*i for i in range(0,len(Vs))],Vs,'b-')
axnew[0].plot([0.025*i for i in range(0,len(Vs))],Vs,'b-')
axarr[1].plot([0.025*i for i in range(0,len(dipoles))],dipoles,'b-')
axnew[1].plot([0.025*i for i in range(0,len(dipoles))],dipoles,'b-')
axarr[2].plot([0.5*i for i in range(0,len(EEGs))],EEGs,'b-')
axnew[2].plot([0.5*i for i in range(0,len(EEGs))],EEGs,'b-')

t_zoom = 3000
dt_zoom = 50

axarr[0].plot([t_zoom-dt_zoom,t_zoom+dt_zoom,t_zoom+dt_zoom,t_zoom-dt_zoom,t_zoom-dt_zoom], [-80,-80,40,40,-80],'r-')
axarr[1].plot([t_zoom-dt_zoom,t_zoom+dt_zoom,t_zoom+dt_zoom,t_zoom-dt_zoom,t_zoom-dt_zoom], [-600,-600,500,500,-600],'r-')
axarr[2].plot([t_zoom-dt_zoom,t_zoom+dt_zoom,t_zoom+dt_zoom,t_zoom-dt_zoom,t_zoom-dt_zoom], [-600,-600,500,500,-600],'r-')

axarr[0].set_ylim([-80,50])
axarr[1].set_ylim([-1200,1200])
axarr[2].set_ylim([-1200,1200])
axnew[0].set_ylim([-80,40])
axnew[1].set_ylim([-600,500])
axnew[2].set_ylim([-600,500])
for i in range(0,3):
  axarr[i].set_xlim([0,7000])
  axarr[i].set_xticks([0,2000,4000,6000])
  #axnew[i].set_xlim([spike_ts[3]-50,spike_ts[3]+50])
  axnew[i].set_xlim([t_zoom-dt_zoom,t_zoom+dt_zoom])

axarr[0].set_yticks([-80,-40,0,40])
axnew[0].set_yticks([-80,-40,0,40])
axnew[1].set_yticks([-500,0])
axnew[2].set_yticks([-500,0])
for i in range(0,3):
  #axnew[i].set_xticks([t_zoom-dt_zoom,t_zoom+dt_zoom])
  axnew[i].set_xticks([])
  axnew[i].set_yticklabels([])
axarr[1].set_yticks([-1000,-500,0,500,-1000])
axarr[2].set_yticks([-1000,-500,0,500,-1000])
axarr[0].set_xticklabels([])
axarr[1].set_xticklabels([])



axarr[0].set_ylabel('mV')
axarr[1].set_ylabel('$\mu$m$\cdot$nA')
axarr[2].set_ylabel('pV')
axarr[2].set_xlabel('$t$ (ms)')

f.text(0.01, 0.93, 'A', fontsize=30)
f.text(0.01, 0.93-0.3*1, 'B', fontsize=30)
f.text(0.01, 0.93-0.3*2, 'C', fontsize=30)

f.savefig("stationary_EEG.eps")
