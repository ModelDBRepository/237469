import mytools
from pylab import *
from neuron import h
import pickle
from os.path import exists
import numpy
import scipy.io
import time

Nmc = 150
seeds = range(1,1000)

oscamp = 0.25

gsyn = 1.07
gNoise = 1.07
myrate = 1.0

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

ft_df = 0.00001      # kHz
ft_maxf = 0.02       # kHz
ft_fs = [ft_df*x for x in range(0,int(round(ft_maxf/ft_df))+1)]

counter = -1
combmutIDs = [1, 0, 2, 3]
for counter in range(0,4):
  oscfreq = float(sys.argv[1])
  combmutID = combmutIDs[counter]
  if not exists('spectrum_freq'+str(oscfreq)+'_comb'+str(counter)+'.sav'):
    FRft = []
    for iseed in range(0,8):
      myseed = seeds[iseed]
      if exists('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_combmutID'+str(combmutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'):
        unpicklefile = open('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_combmutID'+str(combmutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav','r')
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
        spikes = spikes[0]
        sps = array([spikes[i] for i in range(0,len(spikes)) if spikes[i] >= 2000 and spikes[i] < 11000])
        FRftthis = [0.0 for x in ft_fs]
        for ifreq in range(0,len(ft_fs)):
          FRftthis[ifreq] = sum(exp(-2*pi*1j*sps*ft_df*(ifreq-1)))
        FRft.append(FRftthis[:])
      else:
        print 'spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_combmutID'+str(combmutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav not found'
    if len(FRft) > 0:
      picklelist = [ft_fs, FRft]
      file = open('spectrum_freq'+str(oscfreq)+'_comb'+str(counter)+'.sav','w')
      pickle.dump(picklelist,file)
      file.close()
