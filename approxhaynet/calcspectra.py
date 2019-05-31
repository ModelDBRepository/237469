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

    print "igene="+str(igene)+", imut="+str(imut)+", iallmutval="+str(iallmutval)
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
      print mutText


      iters = [0, 2, 6, 8]
      for iiter in range(0,len(iters)):
        iter = iters[iiter]
        counter = counter+1
        if int(sys.argv[1]) != counter:
          continue
        print str(counter)
        oscfreq = float(sys.argv[2])
        if not exists('spectrum_freq'+str(oscfreq)+'_'+str(counter)+'.sav'):
          FRft = []
          for iseed in range(0,15):              
            myseed = seeds[iseed]
            if exists('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_mutID'+str(counter)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav'):
              unpicklefile = open('spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_mutID'+str(counter)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav','r')
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
              print 'spikes_parallel_osc'+str(oscfreq)+'_'+str(Nmc)+'_mutID'+str(counter)+'_'+str(myrate)+'_gNoise'+str(gNoise)+'_gsyn'+str(gsyn)+'_seed'+str(myseed)+'.sav not found'
          if len(FRft) > 0:
            picklelist = [ft_fs, FRft]
            file = open('spectrum_freq'+str(oscfreq)+'_'+str(counter)+'.sav','w')
            pickle.dump(picklelist,file)
            file.close()
        print 'spectrum_freq'+str(oscfreq)+'_'+str(counter)+'.sav exists'
