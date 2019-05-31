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
from matplotlib.collections import PatchCollection

seeds = range(1,6)
oscfreqs = [0.5,0.625,0.75,0.875,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,7.5,10.0,15.0]
iosc_chosen = 6

cols = ['#0000FF','#400090','#900040','#FF0000','#904000','#409000','#00FF00','#30A030','#606060']
col_control = '#2222ff'  

f,axarr = subplots(2,1)
for iosc in range(0,len(oscfreqs)):
  oscfreq = oscfreqs[iosc]
  foundOne = False
  for iseed in range(0,len(seeds)):
    myseed = seeds[iseed]
    mutIDs = [0]
    if iosc == iosc_chosen:
      mutIDs = [0,70]
    for imut in range(0,len(mutIDs)):
      mutID = mutIDs[imut]
      mycolor = col_control
      mystyle = 'b-'
      if imut == 1:
        mycolor = cols[2]
        mystyle = 'b--'
      if exists('spikes_parallel_osc'+str(oscfreq)+'_150_mutID'+str(mutID)+'_1.0_gNoise1.07_gsyn1.07_seed'+str(myseed)+'.sav'):
        unpicklefile = open('spikes_parallel_osc'+str(oscfreq)+'_150_mutID'+str(mutID)+'_1.0_gNoise1.07_gsyn1.07_seed'+str(myseed)+'.sav', 'r')
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
        for ipol in range(0,int(11*oscfreq)+1):
          polygon = Polygon(array([[ipol*1000/oscfreq + x for x in [0,500/oscfreq,500/oscfreq,0,0]],[150*iosc + x for x in [0,0,150,150,0]]]).T, True)
          p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
          p.set_facecolor('#D0D0FF')
          p.set_edgecolor('none')
          axarr[0].add_collection(p)
        if iosc == iosc_chosen:
          axarr[0].plot(spikes[0][::10],[x+iosc*150 for x in spikes[1]][::10],'b.',color='#000000',markersize=0.5)
        else:
          axarr[0].plot(spikes[0][::10],[x+iosc*150 for x in spikes[1]][::10],'b.',markersize=0.5)
        foundOne = True
      else:
        print 'spikes_parallel_osc'+str(oscfreq)+'_150_mutID'+str(mutID)+'_1.0_gNoise1.07_gsyn1.07_seed'+str(myseed)+'.sav not found'
      if iosc == iosc_chosen:
        unpicklefile = open('spectrum_freq'+str(oscfreq)+'_'+str(mutID)+'.sav','r')
        unpickledlist = pickle.load(unpicklefile)
        unpicklefile.close()
        fs = unpickledlist[0]
        ft_fs = unpickledlist[1]
        fts = [mean([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs))]

        unpicklefile = open('spectrum_freq'+str(oscfreq)+'_0.sav','r')
        unpickledlist = pickle.load(unpicklefile)
        unpicklefile.close()
        fs_control = unpickledlist[0]
        ft_fs_control = unpickledlist[1]
        fts_control = [mean([abs(ft_fs_control[isamp][ifreq])**2 for isamp in range(0,len(ft_fs_control))]) for ifreq in range(0,len(fs_control))]
        axarr[1].plot([1000*x for x in fs_control], fts_control, color=col_control)
        axarr[1].plot([1000*x for x in fs], fts, 'b--',color=cols[2])#,linewidth=0.5)
        axnew = f.add_axes([0.27, 0.17, 0.2, 0.25],axisbg='w')
        axnew.plot([1000*x for x in fs_control], fts_control, color=col_control)
        axnew.plot([1000*x for x in fs], fts, 'b--',color=cols[2])#,linewidth=0.5)
    if foundOne:
      break
            
axarr[0].set_xlim([0,11000])
axarr[0].set_xticks([0,2000,4000,6000,8000,10000])
axarr[0].set_xlabel('$t$ (ms)')
axarr[0].set_yticks([])
axarr[0].set_ylabel('Neurons')
for iosc in range(0,len(oscfreqs)):
  axarr[0].plot([0,11000],[iosc*150,iosc*150],'b-',linewidth=0.25,color='#777777')
  axarr[0].text(-1000+10017,iosc*150+23,str(oscfreqs[iosc])+' Hz',fontsize=8,color='#ffffff',fontweight='bold')
  axarr[0].text(-1000+10017,iosc*150+37,str(oscfreqs[iosc])+' Hz',fontsize=8,color='#ffffff',fontweight='bold')
  axarr[0].text(-1000+10000-17,iosc*150+23,str(oscfreqs[iosc])+' Hz',fontsize=8,color='#ffffff',fontweight='bold')
  axarr[0].text(-1000+10000-17,iosc*150+37,str(oscfreqs[iosc])+' Hz',fontsize=8,color='#ffffff',fontweight='bold')
  axarr[0].text(-1000+10000,iosc*150+30,str(oscfreqs[iosc])+' Hz',fontsize=8,color='#ff0000',fontweight='bold')
axarr[0].set_ylim([0,150*17])

axarr[1].set_xlim([0,20])
axarr[1].set_ylim([0,8e7])
axarr[1].plot([1,2,2,1,1],[0,0,5e7,5e7,0],'k-',linewidth=0.25)
axarr[1].set_xlabel('$f$ (Hz)')
axarr[1].set_ylabel('Power')
axarr[1].set_yticks([0,2e7,4e7,6e7,8e7])
          
axnew.set_xlim([1,2])
axnew.set_ylim([0,5e7])
axnew.set_xlabel('$f$ (Hz)',fontsize = 8)
axnew.set_ylabel('Power',fontsize = 8)
for tick in axnew.xaxis.get_major_ticks()+axnew.yaxis.get_major_ticks():
  tick.label.set_fontsize(6)
t = axnew.yaxis.get_offset_text()
t.set_size(6)

axarr[0].set_position([0.125, 0.536363636364, 0.375, 0.9-0.536363636364])
axarr[1].set_position([0.125, 0.1, 0.375, 0.463636363636-0.1])

f.text(0.07, 0.86, 'A', fontsize=28)
f.text(0.07, 0.43, 'B', fontsize=28)
f.savefig('fig2ab.eps')
