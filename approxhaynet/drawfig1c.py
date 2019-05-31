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

variants = [[0,5,1],[2,4,7],[1,2,13],[3,1,0],[5,0,0],[8,3,0],[12,2,0],[13,4,0]] #from drawallmeangains (maxCountersAll) #Removed KCNN3!

rates = [0.1*x for x in range(4,17)]


icell = 0
gsyn = 1.07
gNoise = 1.07

cols = ['#666666','#012345','#cc00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
col_control = '#2222ff'

if exists('nSpikes_control.sav'):
  unpicklefile = open('nSpikes_control.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  nSpikes_control = unpickledlist[0]
else:
  nSpikes_control = []
  for irate in range(0,len(rates)):
    myrate = rates[irate]
    nSpikesThisRate = []
    print "loading control rate="+str(myrate)
    for myseed in range(1,10):
      if exists('spikes_parallel150_mutID0_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'):
        unpicklefile = open('spikes_parallel150_mutID0_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav', 'r')
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
        nSpikesThisRate.append(len(spikes[0]))
      else:
        print 'spikes_parallel150_mutID0_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav not found'
    nSpikes_control.append(nSpikesThisRate[:])
  file = open('nSpikes_control.sav', 'w')
  pickle.dump([nSpikes_control],file)
  file.close()

for counter in range(0,461):
  if not exists('nSpikes'+str(counter)+'.sav'):
    nSpikesThisIter = []
    print "loading mutID="+str(counter)
    for irate in range(0,len(rates)):
      myrate = rates[irate]
      nSpikesThisRate = []
      print "loading rate="+str(myrate)
      for myseed in range(1,10):
        if exists('spikes_parallel150_mutID'+str(counter)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'):
          unpicklefile = open('spikes_parallel150_mutID'+str(counter)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav', 'r')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
        else:
          continue
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

        nSpikesThisRate.append(len(spikes[0]))
        if len(spikes[0]) < 20:
          print str(spikes[0])
      nSpikesThisIter.append(nSpikesThisRate[:])
      file = open('nSpikes'+str(counter)+'.sav', 'w')
      pickle.dump([nSpikesThisIter],file)
      file.close()

print "Loading control f-I curves..."
ispDef = 0
Is = [0.35+0.05*x for x in range(0,22)]
unpicklefile = open('../haymod/ifcurvesmut_cs'+str(icell)+'_0_0_0.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spTimesThisMutVal = unpickledlist[1+ispDef]
nSpikes_control_if = [sum([1 for x in spTimesThisMutVal[5][j] if x >= 500]) for j in range(0,len(Is))]

counter = 0

close("all")                               
f, axarr = plt.subplots(2, len(variants)+1)              
lenvarper2 = 4 # len(variants)/2
for ix in range(0,lenvarper2):
  for iy in range(0,2):
    axarr[0,ix+lenvarper2*iy].set_position([0.1+0.176*ix, 0.1+0.44*(1-iy), 0.176, 0.37])
    axarr[1,ix+lenvarper2*iy].set_position([0.121+0.176*ix, 0.33+0.44*(1-iy), 0.075, 0.11])

axarr[0,8].set_position([0.1+0.176*4, 0.1+0.44*(1-1), 0.176, 0.37])
axarr[1,8].set_position([0.121+0.176*4, 0.33+0.44*(1-1), 0.075, 0.11])

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
            newVal =  [x+mutvals for x in defVals[mutvar]]
            if mutvals >= 0 and kmutvar==0:
              mutText = mutText + "+" + str(mutvals) +" mV"
            elif kmutvar==0:
              mutText = mutText  + str(mutvals) +" mV"
          else:
            newVal =  [x*mutvals for x in defVals[mutvar]]
            if kmutvar==0:
              mutText = mutText + "*" + str(mutvals)
          if kmutvar < len(mutvars)-1:
            mutText = mutText + ", "

      iters = [0, 2, 6, 8]
      iters_if = [0, 2, 5, 6, 8, -1]
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
        if exists('nSpikes'+str(counter)+'.sav'):
          unpicklefile = open('nSpikes'+str(counter)+'.sav', 'r')
          unpickledlist = pickle.load(unpicklefile)
          unpicklefile.close()
          nSpikesThisIter = unpickledlist[0]
        else:
          print 'nSpikes'+str(counter)+'.sav not found!'
          continue
        axarr[0,ivar].plot(rates, [mean(x)/150./11.0 for x in nSpikesThisIter], 'b-',color=cols[iter])
        print "igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)+", ivar="+str(ivar)+", counter="+str(counter)+", iter="+str(iter)
      if doSkip:
        continue
      print mutText

      print "Loading f-I curves..."
      unpicklefile = open('../haymod/ifcurvesmut_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      spTimesThisMutVal = unpickledlist[1+ispDef]
      for iiter in range(0,len(iters_if)):
        iter = iters_if[iiter]
        if iter == 5 or iter == -1:
          continue
        nSpikes_if = [sum([1 for x in spTimesThisMutVal[iiter][j] if x >= 500]) for j in range(0,len(Is))]
        axarr[1,ivar].plot(Is, [x/7.5 for x in nSpikes_if], 'b-', color=cols[iter])
      axarr[1,ivar].plot(Is, [x/7.5 for x in nSpikes_control_if], 'b-', color=col_control)
      
      axarr[0,ivar].plot(rates, [mean(x)/150.0/11. for x in nSpikes_control], 'b-',color=col_control)
      axarr[0,ivar].set_title(geneNames[igene])

      axarr[0,ivar].set_xlim([0.4,1.6])
      axarr[0,ivar].set_ylim([0,12])
      axarr[0,ivar].set_xticks([0.6, 1.0, 1.4])
      axarr[0,ivar].set_yticks([0, 3, 6, 9, 12])
      if ivar < lenvarper2:
        axarr[0,ivar].set_xticklabels(['', '', ''])
      elif ivar==6:
        #axarr[0,ivar].set_xlabel('rate factor r                                ')
        axarr[0,ivar].set_xlabel('rate factor $r$')
        #axarr[0,ivar].set_xlabel('rate factor $c_{\mathrm{rate}}$                            ')
      if ivar % lenvarper2 > 0:
        axarr[0,ivar].set_yticklabels(['', '', '', '', ''])        
      else:
        axarr[0,ivar].set_ylabel('$f$ (Hz)')

      axarr[1,ivar].set_xlim([0.35,1.4])
      axarr[1,ivar].set_ylim([0,20])
      axarr[1,ivar].set_xticks([0.5,1.0])
      axarr[1,ivar].set_yticks([0,10,20])
      axarr[1,ivar].set_xlabel('$I$ (nA)',fontsize=8)
      for tick in axarr[1,ivar].xaxis.get_major_ticks()+axarr[1,ivar].yaxis.get_major_ticks():
        tick.label.set_fontsize(6)
      axarr[1,ivar].yaxis.set_tick_params(pad=3)      

      f.savefig("fig1c.eps")

iters = [0, 2, 6, 8]
iters_if = [0, 2, 6, 8]
combmutIDnums = [1, 0, 2, 3]
ivar = 8

for iiter in range(0,len(iters)):
  iter = iters[iiter]
  combmutIDnum = combmutIDnums[iiter]
  if not exists('nSpikes_comb_'+str(iiter)+'.sav'):
    nSpikesThisIter = []
    print "loading mutID="+str(combmutIDnum)
    for irate in range(0,len(rates)):
      myrate = rates[irate]
      nSpikesThisRate = []
      print "loading rate="+str(myrate)
      for myseed in range(1,12):
        if exists('spikes_parallel150_mutcombID'+str(combmutIDnum)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'):
          unpicklefile = open('spikes_parallel150_mutcombID'+str(combmutIDnum)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav', 'r')
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
          nSpikesThisRate.append(len(spikes[0]))
          if len(spikes[0]) < 20:
            print str(spikes[0])
        else:
          print 'spikes_parallel150_mutcombID'+str(combmutIDnum)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav not found'
      nSpikesThisIter.append(nSpikesThisRate[:])
      file = open('nSpikes_comb_'+str(iiter)+'.sav', 'w')
      pickle.dump([nSpikesThisIter],file)
      file.close()

  if exists('nSpikes_comb_'+str(iiter)+'.sav'):
    unpicklefile = open('nSpikes_comb_'+str(iiter)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    nSpikesThisIter = unpickledlist[0]
  else:
    print 'nSpikes_comb_'+str(iiter)+'.sav not found!'
    continue
  axarr[0,ivar].plot(rates, [mean(x) for x in nSpikesThisIter], 'b-',color=cols[iter])

unpicklefile = open('../haymod/ifcurves_comb_cs'+str(icell)+'.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
spTimesThisMutVal = unpickledlist[1+ispDef]

for iiter in range(0,len(iters_if)):
  iter = iters_if[iiter]
  if iter == 5 or iter == -1:
    continue
  nSpikes_if = [sum([1 for x in spTimesThisMutVal[iiter][j] if x >= 500]) for j in range(0,len(Is))]
  axarr[1,ivar].plot(Is, [x/7.5 for x in nSpikes_if], 'b-', color=cols[iter])
  print "Comb nSpikes="+str(nSpikes_if)

axarr[1,ivar].plot(Is, [x/7.5 for x in nSpikes_control_if], 'b-', color=col_control)

axarr[0,ivar].plot(rates, [mean(x) for x in nSpikes_control], 'b-',color=col_control)
axarr[0,ivar].set_title("Combination")
axarr[0,ivar].set_xlim([0.4,1.6])
axarr[0,ivar].set_ylim([0,20000])
axarr[0,ivar].set_xticks([0.6, 1.0, 1.4])
axarr[0,ivar].set_yticks([0, 5000, 10000, 15000, 20000])
axarr[0,ivar].set_yticklabels(['', '', '', '', ''])        

axarr[1,ivar].set_xlim([0,1.4])
axarr[1,ivar].set_ylim([0,20])
axarr[1,ivar].set_xticks([0,0.5,1.0])
axarr[1,ivar].set_yticks([0,10,20])
axarr[1,ivar].set_xlabel('$I$ (nA)',fontsize=8)
for tick in axarr[1,ivar].xaxis.get_major_ticks()+axarr[1,ivar].yaxis.get_major_ticks():
  tick.label.set_fontsize(6)
      
f.text(0.025, 0.88, 'C', fontsize=36)
f.savefig("fig1c.eps")
