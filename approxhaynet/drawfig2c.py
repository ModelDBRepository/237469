import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import sys
import time
from os.path import exists

coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
import mutation_stuff
MT = mutation_stuff.getMT()
defVals = mutation_stuff.getdefvals()
geneNames = mutation_stuff.getgenenames()
keyList = defVals.keys()
for idefval in range(0,len(keyList)):
  if type(defVals[keyList[idefval]]) is not list:
    defVals[keyList[idefval]] = [defVals[keyList[idefval]], defVals[keyList[idefval]]] #make the dictionary values [somatic, apical]
updatedVars = ['somatic','apical','basal'] # the possible classes of segments that defVals may apply to
whichDefVal = [0,1,0]                      # use the defVal[0] for somatic and basal segments and defVal[1] for apical segments

variants = [[0,5,1],[2,4,7],[1,2,13],[3,1,0],[5,0,0],[8,3,0],[12,2,0],[13,4,0]] #from drawallmeangains (maxCountersAll)

oscamp = 0.25
## These are the frequencies used in the article figure 2:
#oscfreqs = [0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5, 10.0, 15.0]
# This is a smaller set of frequencies that captures most of the shape of the response curve:
oscfreqs = [0.5, 0.625, 0.75, 0.875, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.5, 10.0, 15.0]

gsyn = 1.07
gNoise = 1.07

fts_control = []
ft_df = 0.00001

if exists('spectrum_freq1.0_0.sav'): #Do this just to get the fs_control
  unpicklefile = open('spectrum_freq1.0_0.sav','r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  fs_control = unpickledlist[0] 

print "Loading control..."
amps_control = []
amps_control_std = []
if exists('spectra_fig2_0.sav'):
  unpicklefile = open('spectra_fig2_0.sav','r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  fts = unpickledlist[0]
  fts_std = unpickledlist[1]
  amps_control = unpickledlist[2]
  amps_control_std = unpickledlist[3]
else:
  for iosc in range(0,len(oscfreqs)):
    oscfreq = oscfreqs[iosc]
    if exists('spectrum_freq'+str(oscfreq)+'_0.sav'):
      unpicklefile = open('spectrum_freq'+str(oscfreq)+'_0.sav','r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      fs_control = unpickledlist[0]
      ft_fs = unpickledlist[1]
      fts = [mean([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs_control))]
      fts_std = [std([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs_control))]
      iosc_fs = next(i for i,x in enumerate(fs_control) if x>oscfreq/1000)
      amps_control.append(fts[iosc_fs])
      amps_control_std.append(fts_std[iosc_fs])
    else:
      print 'spectrum_freq'+str(oscfreq)+'_0.sav not found!'
  print "Loading done"
  picklelist = [fts, fts_std, amps_control, amps_control_std]
  file = open('spectra_fig2_0.sav','w')
  pickle.dump(picklelist,file)
  file.close()

styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
cols = ['#666666','#012345','#cc00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
col_control = '#2222ff'  

counter = -1

close("all")                               
f, axarr = plt.subplots(1, len(variants)+1)
lenvarper2 = 4
for ix in range(0,lenvarper2):
  for iy in range(0,2):
    axarr[ix+lenvarper2*iy].set_position([0.08+0.176*ix, 0.1+0.44*(1-iy), 0.176, 0.37])
axarr[8].set_position([0.08+0.176*4, 0.1+0.44*(1-1), 0.176, 0.37])

timesAll = []
VsomaAll = []
spikeFreqsAll = []

for igene in range(0,len(MT)):
  for imut in range(0,len(MT[igene])):
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    for imutvar in range(0,len(MT[igene][imut])):
      thesemutvars.append(MT[igene][imut][imutvar][0])
      if type(MT[igene][imut][imutvar][1]) is int or type(MT[igene][imut][imutvar][1]) is float:
        MT[igene][imut][imutvar][1] = [MT[igene][imut][imutvar][1]]
      nVals[imutvar] = len(MT[igene][imut][imutvar][1])
    cumprodnVals = cumprod(nVals)
    allmutvars = cumprodnVals[len(MT[igene][imut])-1]*[thesemutvars]
    allmutvals = []
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      allmutvals.append([0]*len(thesemutvars))
    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar==0:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][iallmutval%nVals[imutvar]]
        else:
          allmutvals[iallmutval][imutvar] = MT[igene][imut][imutvar][1][(iallmutval/cumprodnVals[imutvar-1])%nVals[imutvar]]

    for iallmutval in range(0,cumprodnVals[len(MT[igene][imut])-1]):
      mutval = allmutvals[iallmutval]

      thisCoeff = 1
      mutText = ""
      for imutvar in range(0,len(MT[igene][imut])):
        if imutvar > 0 and imutvar%2==0:
          mutText = mutText+"\n"
        mutvars = allmutvars[iallmutval][imutvar]
        mutvals = allmutvals[iallmutval][imutvar]
        if type(mutvars) is str:
          mutvars = [mutvars]
        mutText = mutText + str(mutvars) + ": "
        for kmutvar in range(0,len(mutvars)):
          mutvar = mutvars[kmutvar]
          if mutvar.find('offm') > -1 or mutvar.find('offh') > -1 or mutvar.find('ehcn') > -1:
            newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
            if mutvals >= 0 and kmutvar==0:
              mutText = mutText + "+" + str(mutvals) +" mV"
            elif kmutvar==0:
              mutText = mutText  + str(mutvals) +" mV"
          else:
            newVal =  [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
            if kmutvar==0:
              mutText = mutText + "*" + str(mutvals)
          if kmutvar < len(mutvars)-1:
            mutText = mutText + ", "

      if igene == 0 and imut == 0 and iallmutval == 0:
        iters = [-1, 0, 2, 6, 8]
      else:
        iters = [0, 2, 6, 8]
      doSkip = True
      ivar = -1
      for iiter in range(0,len(iters)):
        iter = iters[iiter]
        counter = counter+1
        for ivar2 in range(0,len(variants)):
          if igene == variants[ivar2][0] and imut == variants[ivar2][1] and iallmutval == variants[ivar2][2]:
            ivar = ivar2
            doSkip = False
            break
        if doSkip:
          continue
        
        if exists('spectra_fig2_'+str(counter)+'.sav'):
          unpicklefile = open('spectra_fig2_'+str(counter)+'.sav','r')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
          fts = unpickledlist[0]
          fts_std = unpickledlist[1]
          amps = unpickledlist[2]
          amps_std = unpickledlist[3]
        else:
          amps = []
          amps_std = []
          for iosc in range(0,len(oscfreqs)):
            oscfreq = oscfreqs[iosc]
            iosc_fs = next(i for i,x in enumerate(fs_control) if x>oscfreq/1000)
            if exists('spectrum_freq'+str(oscfreq)+'_'+str(counter)+'.sav'):
              unpicklefile = open('spectrum_freq'+str(oscfreq)+'_'+str(counter)+'.sav','r')
              unpickledlist = pickle.load(unpicklefile)
              unpicklefile.close()
              fs = unpickledlist[0]
              ft_fs = unpickledlist[1]
              fts = [mean([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs))]
              fts_std = [std([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs))]
              amps.append(fts[iosc_fs])
              amps_std.append(fts_std[iosc_fs])
            else:
              print 'spectrum_freq'+str(oscfreq)+'_'+str(counter)+'.sav does not exist' 
              amps.append([])
              amps_std.append([])
          picklelist = [fts, fts_std, amps, amps_std]
          file = open('spectra_fig2_'+str(counter)+'.sav','w')
          pickle.dump(picklelist,file)
          file.close()
        if len(amps) > 0:
          if not any([type(amps[iosc]) is list for iosc in range(0,len(amps))]):
            axarr[ivar].semilogx(oscfreqs, [amps[iosc] for iosc in range(0,len(amps))], 'b-', color=cols[iter])
          else:
            myvec = zeros([len(amps),])
            for iosc in range(0,len(amps)):
              if type(amps[iosc]) is not list:
                myvec[iosc] = amps[iosc]
              else: 
                myvec[iosc] = nan
            axarr[ivar].semilogx(oscfreqs, myvec, 'b.-', color=cols[iter])
            print "Some frequencies not found in counter="+str(counter)+", iter="+str(iter)+"!"
        else:
          print "No frequencies not found in counter="+str(counter)+", iter="+str(iter)+"!"
        print "igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)+", ivar="+str(ivar)+", counter="+str(counter)+", iter="+str(iter)
      if doSkip:
        continue
      print mutText

      axarr[ivar].semilogx(oscfreqs, amps_control, 'b-', color=col_control)
      axarr[ivar].set_title(geneNames[igene])
      axarr[ivar].set_xlim([0.4,15])
      axarr[ivar].set_ylim([0,7.5e7])
      axarr[ivar].set_yticks([0, 3e7, 6e7])
      if ivar < lenvarper2:
        axarr[ivar].set_xticklabels(['']*len(axarr[ivar].get_xticks()))
      elif ivar == 6:
        axarr[ivar].set_xlabel('Input frequency $f$ (Hz)')
      if ivar % lenvarper2 > 0:
        axarr[ivar].set_yticklabels(['', '', ''])        
      elif ivar == 0:
        axarr[ivar].set_ylabel('Power of the frequency component corresponding $f$                                            ')
  f.savefig("fig2c.eps")

iters = [0, 2, 6, 8]
combmutIDnums = [0, 1, 2, 3]
ivar = 8

amps_thismut = []
for iiter in range(0,4):
  iter = iters[iiter]
  if exists('spectra_fig2_comb'+str(combmutIDnums[iiter])+'.sav'):
    unpicklefile = open('spectra_fig2_comb'+str(combmutIDnums[iiter])+'.sav','r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    fts = unpickledlist[0]
    fts_std = unpickledlist[1]
    amps = unpickledlist[2]
    amps_std = unpickledlist[3]
  else:
    amps = []
    amps_std = []
    for iosc in range(0,len(oscfreqs)):
      oscfreq = oscfreqs[iosc]
      iosc_fs = next(i for i,x in enumerate(fs_control) if x>oscfreq/1000)
      if exists('spectrum_freq'+str(oscfreq)+'_comb'+str(combmutIDnums[iiter])+'.sav'):
        unpicklefile = open('spectrum_freq'+str(oscfreq)+'_comb'+str(combmutIDnums[iiter])+'.sav','r')
        unpickledlist = pickle.load(unpicklefile)
        unpicklefile.close()
        fs = unpickledlist[0]
        ft_fs = unpickledlist[1]
        fts = [mean([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs))]
        fts_std = [std([abs(ft_fs[isamp][ifreq])**2 for isamp in range(0,len(ft_fs))]) for ifreq in range(0,len(fs))]
        amps.append(fts[iosc_fs])
        amps_std.append(fts_std[iosc_fs])
      else:
        print 'spectrum_freq'+str(oscfreq)+'_comb'+str(combmutIDnums[iiter])+'.sav does not exist'
        amps.append(nan)
        amps_std.append(nan)
    picklelist = [fts, fts_std, amps, amps_std]
    file = open('spectra_fig2_comb'+str(combmutIDnums[iiter])+'.sav','w')
    pickle.dump(picklelist,file)
    file.close()

  if len(amps) > 0:
    if not any([type(amps[iosc]) is list for iosc in range(0,len(amps))]):
      axarr[ivar].semilogx(oscfreqs, [amps[iosc] for iosc in range(0,len(amps))], 'b-', color=cols[iter])
    else:
      myvec = zeros([len(amps),])
      for iosc in range(0,len(amps)):
        if type(amps[iosc]) is not list:
          myvec[iosc] = amps[iosc]
        else: 
          myvec[iosc] = nan
      axarr[ivar].semilogx(oscfreqs, myvec, 'b.-', color=cols[iter])
      print "Some frequencies not found in counter="+str(counter)+", iter="+str(iter)+"!"
  else:
    print "No frequencies not found in counter="+str(counter)+", iter="+str(iter)+"!"

  amps_thismut.append(amps[:])
axarr[ivar].semilogx(oscfreqs, amps_control, 'b-', color=col_control)
axarr[ivar].set_ylim([0,8e7])

axarr[ivar].set_title("Combination")
axarr[ivar].set_xlim([0.4,15])
axarr[ivar].set_ylim([0,7.5e7])
axarr[ivar].set_yticks([0, 3e7, 6e7])
axarr[ivar].set_yticklabels(['', '', ''])        

for ivar in range(0,9):
  t = axarr[ivar].yaxis.get_offset_text()
  if type(t) is matplotlib.text.Text:
    if t.get_text()[0:2] == '1e':
      t.set_text('$\\times 10^{'+t.get_text()[2:]+'}$')
    t.set_position((t.get_position()[0]-0.15,t.get_position()[1]))
f.text(0.01, 0.9, 'C', fontsize=31)
f.savefig("fig2c.eps")
