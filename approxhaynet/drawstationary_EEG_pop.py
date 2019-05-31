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
  axnew.append(f.add_axes([0.18, 0.7-0.3*i+0.18, 0.18, 0.09],axisbg='w'))

unpicklefile = open('spikes_parallel150_mutID0_1.0_gNoise1.07_gsyn1.07_seed1_withLFP.sav', 'r') #from runmanycellsLFP.py
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Nplaced = 0
spikedCells_all = []
for j in range(0,len(unpickledlist[1])):
  spikedCells = unpickledlist[1][j]
  spikedCellsUnique = unique(spikedCells)
  spikedCells2 = zeros(spikedCells.shape)
  for i in range(0,len(spikedCellsUnique)):
    spikedCells2[spikedCells == spikedCellsUnique[i]] = Nplaced + i
  Nplaced = Nplaced + len(spikedCellsUnique)
  spikedCells_all = hstack([spikedCells_all, spikedCells2])
spikes = [hstack(unpickledlist[0]),spikedCells_all]

dipoles = [x[1] for x in unpickledlist[2]]

unpicklefile = open('spikes_parallel150_mutID0_1.0_gNoise1.07_gsyn1.07_seed1_withLFP_EEG.sav', 'r') #from runsinglecellLFP.py and calcEEG.py
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
EEGs = unpickledlist[0][0]

axarr[0].plot(spikes[0],spikes[1],'b.',markersize=0.5,mew=0.5)
axnew[0].plot(spikes[0],spikes[1],'b.',markersize=1.0,mew=1.0)
axarr[1].plot([0.025*i for i in range(0,len(dipoles))],dipoles,'b-',linewidth=0.125)
axnew[1].plot([0.025*i for i in range(0,len(dipoles))],dipoles,'b-')
axarr[2].plot([0.5*i for i in range(0,len(EEGs))],EEGs,'b-',linewidth=0.125)
axnew[2].plot([0.5*i for i in range(0,len(EEGs))],EEGs,'b-')

t_zoom = 5550
dt_zoom = 300
axarr[0].plot([t_zoom-dt_zoom,t_zoom+dt_zoom,t_zoom+dt_zoom,t_zoom-dt_zoom,t_zoom-dt_zoom], [0,0,149.5,149.5,0],'r-')
axarr[1].plot([t_zoom-dt_zoom,t_zoom+dt_zoom,t_zoom+dt_zoom,t_zoom-dt_zoom,t_zoom-dt_zoom], [-8000,-8000,-2000,-2000,-8000],'r-')
axarr[2].plot([t_zoom-dt_zoom,t_zoom+dt_zoom,t_zoom+dt_zoom,t_zoom-dt_zoom,t_zoom-dt_zoom], [-8000,-8000,-2000,-2000,-8000],'r-')

axarr[0].set_ylim([-1,150])
axarr[1].set_ylim([-8500,1500])
axarr[2].set_ylim([-8500,1500])
axnew[0].set_ylim([-1,150])
axnew[1].set_ylim([-8000,-2000])
axnew[2].set_ylim([-8000,-2000])

axnew[0].set_yticks([])
axnew[1].set_yticks([-8000,-4000])
axnew[2].set_yticks([-8000,-4000])
for i in range(0,3):
  #axnew[i].set_xticks([t_zoom-dt_zoom,t_zoom+dt_zoom])
  axnew[i].set_xticks([])
  axnew[i].set_yticklabels([])
axarr[0].set_yticks([0,50,100,150])
axarr[1].set_yticks([-8000,-4000,0])
axarr[2].set_yticks([-8000,-4000,0])
axarr[0].set_xticklabels([])
axarr[1].set_xticklabels([])

for i in range(0,3):
  axarr[i].set_xlim([0,7000])
  axarr[i].set_xticks([0,2000,4000,6000])
  axnew[i].set_xlim([t_zoom-dt_zoom,t_zoom+dt_zoom])

axarr[0].set_ylabel('')
axarr[1].set_ylabel('$\mu$m$\cdot$nA')
axarr[2].set_ylabel('pV')
axarr[2].set_xlabel('$t$ (ms)')

f.text(0.01, 0.93, 'D', fontsize=30)
f.text(0.01, 0.93-0.3*1, 'E', fontsize=30)
f.text(0.01, 0.93-0.3*2, 'F', fontsize=30)

f.savefig("stationary_EEG_pop.eps")
