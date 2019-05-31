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

rates = [0.1*x for x in range(4,17)]
gsyn = 1.07
gNoise = 1.07

irate_chosen = 6

cols = ['#666666','#012345','#cc00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
col_control = '#2222ff'  

seeds = [1,2,3,4,5]
gauss_std = 25 #ms

f,axarr = subplots(2,1)
for irate in range(0,len(rates)):
  myrate = rates[irate]
  foundOne = False
  for iseed in range(0,len(seeds)):
    myseed = seeds[iseed]
    mutIDs = [0]
    if irate == irate_chosen:
      mutIDs = [0,70]
    for imut in range(0,len(mutIDs)):
      mutID = mutIDs[imut]
      mycolor = col_control
      mystyle = 'b-'
      if imut == 1:
        mycolor = cols[2]
        mystyle = 'b--'
      if exists('spikes_parallel150_mutID'+str(mutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'):
        print 'loading spikes_parallel150_mutID'+str(mutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'
        unpicklefile = open('spikes_parallel150_mutID'+str(mutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav', 'r')
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
        if irate == irate_chosen:
          axarr[0].plot(spikes[0][::5],[x+irate*150 for x in spikes[1]][::5],'b.',color='#000000',markersize=0.5)
        else:
          axarr[0].plot(spikes[0][::5],[x+irate*150 for x in spikes[1]][::5],'b.',markersize=0.5)
        foundOne = True
      else:
        print 'spikes_parallel150_mutID'+str(mutID)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav not found'
      if irate == irate_chosen:
        nSpikesBins = zeros([11001,])
        nSpikesBins_conv = zeros([11001,])
        tBins = [0.5+x for x in range(0,11001)]
        print "Calculating nSpikes"
        for ispike in range(0,len(spikes[0])):
          nSpikesBins[int(spikes[0][ispike])] = nSpikesBins[int(spikes[0][ispike])] + 1
        for ibin in range(0,len(nSpikesBins)):
          nSpikesBins_conv[ibin] = sum(exp(-1/2*(tBins[ibin]-spikes[0])**2/gauss_std**2))/sqrt(pi*gauss_std**2)
        print "Done"

        axarr[1].plot(tBins, nSpikesBins_conv/150*1000, mystyle, color=mycolor)
    if foundOne:
      break
axarr[0].set_xlim([0,11000])
axarr[0].set_xticks([0,2000,4000,6000,8000,10000])
axarr[0].set_xlabel('$t$ (ms)')
axarr[0].set_yticks([])
axarr[0].set_ylabel('neurons')
for irate in range(0,len(rates)):
      axarr[0].plot([0,11000],[irate*150,irate*150],'b-',linewidth=0.25,color='#777777')
      axarr[0].text(-670+10017,irate*150+38,'r = '+str(rates[irate]),fontsize=8,color='#ffffff',fontweight='bold')
      axarr[0].text(-670+10017,irate*150+52,'r = '+str(rates[irate]),fontsize=8,color='#ffffff',fontweight='bold')
      axarr[0].text(-670+10000-17,irate*150+38,'r = '+str(rates[irate]),fontsize=8,color='#ffffff',fontweight='bold')
      axarr[0].text(-670+10000-17,irate*150+52,'r = '+str(rates[irate]),fontsize=8,color='#ffffff',fontweight='bold')
      axarr[0].text(-670+10000,irate*150+45,'r = '+str(rates[irate]),fontsize=8,color='#ff0000',fontweight='bold')
axarr[0].set_ylim([0,150*13])

axarr[1].set_xlim([0,11000])
axarr[1].set_ylim([0,16])
axarr[1].set_yticks([0,5,10,15])
axarr[1].set_xlabel('$t$ (ms)')
axarr[1].set_ylabel('$f$ (Hz)')
axarr[0].set_position([0.125, 0.536363636364, 0.375, 0.9-0.536363636364])
axarr[1].set_position([0.125, 0.1, 0.375, 0.463636363636-0.1])
          
print 'saving fig1ab.eps'
          
f.text(0.05, 0.86, 'A', fontsize=28)
f.text(0.05, 0.45, 'B', fontsize=28)
f.savefig('fig1ab.eps')
