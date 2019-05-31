#cp ../haymod3c/drawppiranges.py drawppiranges.py
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
#import matplotlib.pyplot
#import matplotlib.lines


morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
fs = 8
DI = 70 #distance between images

import mutation_stuff
MT = mutation_stuff.getMT()
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

spTimesAll = []
styles = ['g-','g-','g-','g-','g-','g-','g-','g-','g-']
#cols = ['#aaffaa','#aaffaa','#66ff66','#66aaaa','#00aaaa','#00aaaa']
#cols = ['#00aaaa','#00bb77','#11cc44','#11dd11','#55ee00','#99dd00']
#cols = ['#00aaaa','#11cc44','#55ee00','#bbaa00','#ee6600','#ff0000', '#aa00aa','#772277','#333333']
#cols = ['#666666','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#00aaaa','#772277','#00cc00']
cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
yplus = [1, 2, 3, 4, 5, 6, -1, -2, -3]
yplus = [x+3 for x in yplus]
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]

params = {'text.latex.preamble': [r"\usepackage{upgreek}"],
          'text.usetex': True}
plt.rcParams.update(params)

genelabels = ['CACNA1C (Kudrnac et al. 2009)','CACNA1C (Kudrnac et al. 2009)','CACNB2 (Cordeiro et al. 2009)','CACNB2 (Massa et al. 1995)','CACNB2 (Link et al. 2009)',
          'CACNA1D (Tan et al. 2011; \n Bock et al. 2011)','CACNA1D (Tan et al. 2011; \n Bock et al. 2011)','CACNA1D (Zhang et al. 2011; \n Perez-Alvarez et al. 2011)',
          'CACNA1I (Murbartian et al. 2004)','CACNA1S (Pirone et al. 2010)','CACNA1S (Tuluc et al. 2009)','ATP2A2 (Ji et al. 2000)','ATP2B2 (Fakira et al. 2012)',
          'ATP2B2 (Empson et al. 2010)','ATP2B2 (Ficarella et al. 2007)','SCN1A (Cestele et al. 2008)','SCN1A (Vanmolkot et al. 2007)','SCN9A (Estacion et al. 2011)',
          'SCN9A (Estacion et al. 2008)','SCN9A (Han et al. 2006)','SCN9A (Dib-Hajj et al. 2005)','KCNS3 (Shepard and Rae 1999)','KCNB1 (Bocksteins et al. 2011)',
          'KCNB1 (Bocksteins et al. 2011)','KCNB1 (Bocksteins et al. 2011)','KCNB1 (Bocksteins et al. 2011)','KCNB1 (Bocksteins et al. 2011)',
          'KCNB1 (Bocksteins et al. 2011)','KCNN3 (Wittekindt et al. 2004)','HCN1 (Ishii et al. 2007)','','','','','','','','','','','','','','','','','','','','','','','','']
labelxplus = [0,0,0,0,0,-0.6,-0.6,-0.6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
labelyplus = [0,0,0,1,0,-6,-6,-6,-3,1,0,0,0,1,0,3,0,1,1,1,0,1,0,0,0,0,0,0,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]


unpicklefile = open('thresholddistalamp300_control.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
gmaxes = unpickledlist[0]
Nsyns = unpickledlist[1]
maxSynsPerSeg = unpickledlist[2]

unpicklefile = open('ppiextraspthrcoeff300_relthr_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
gsAllAll = unpickledlist[1]
PPIdts = unpickledlist[2]
DataAllAll = unpickledlist[3]



for icell in range(0,1):
  theseCoeffsAll = theseCoeffsAllAll[icell]
  close("all")
  f, axarr = plt.subplots(1, 1)
  counter = -1
  labelcounter = -1

  thisy = [x*1.1 for x in DataAllAll[icell][0][0][0][5]]
  firstabovethree = -1
  lastabovethree = -1
  for iy in range(0,len(thisy)):
    if thisy[iy] > 3 and firstabovethree == -1:
      firstabovethree = iy
    if thisy[iy] > 3:
      lastabovethree = iy
  xplusthis1 = 2.0*(3.0-thisy[firstabovethree-1])/(thisy[firstabovethree]-thisy[firstabovethree-1])
  xplusthis2 = 2.0*(3.0-thisy[lastabovethree])/(thisy[lastabovethree+1]-thisy[lastabovethree])
  range_control = [PPIdts[firstabovethree-1]+xplusthis1,PPIdts[lastabovethree]+xplusthis2]
  
  for igene in [0,1,2,3,4,5,6,8,9,10,13,11,12]:
   for imut in range(0,len(MT[igene])):
    labelcounter = labelcounter+1
    nVals = len(MT[igene][imut])*[0]
    thesemutvars = []
    theseCoeffs = theseCoeffsAll[igene][imut]
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
      counter = counter+1
      #print str(DataAllAll[icell][igene][imut][iallmutval])
      #if igene != 13:
      #  unpicklefile = open('../haymod3/steadystate2_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      #else:
      #  unpicklefile = open('steadystate2_cs'+str(icell)+'_'+str(igene)+'_'+str(imut)+'_'+str(iallmutval)+'.sav', 'r')
      #unpickledlist = pickle.load(unpicklefile)
      #unpicklefile.close()

      iters = [0, 2, 5, 6, 8]
      yplus = [0, 1, 2, 0, 1]
      xplus = [0, 0, 0, 0.3, 0.3]
      for iiter in range(0,len(iters)):
        iter = iters[iiter]
        thisy = [x*1.1 for x in DataAllAll[icell][igene][imut][iallmutval][iiter]]
        firstabovethree = -1
        lastabovethree = -1
        for iy in range(0,len(thisy)):
          if thisy[iy] > 3 and firstabovethree == -1:
            firstabovethree = iy
          if thisy[iy] > 3:
            lastabovethree = iy
        if firstabovethree == -1:
          thisRange = [nan,nan]
        else:
          xplusthis1 = 2.0*(3.0-thisy[firstabovethree-1])/(thisy[firstabovethree]-thisy[firstabovethree-1])
          xplusthis2 = 2.0*(3.0-thisy[lastabovethree])/(thisy[lastabovethree+1]-thisy[lastabovethree])
          thisRange = [PPIdts[firstabovethree-1]+xplusthis1,PPIdts[lastabovethree]+xplusthis2]
        
        if igene==0 and imut==0 and iallmutval==0 and iiter < 3:
          axarr.plot([-2,-2], [x+yplus[iiter]*DI for x in range_control], 'b-',color='#0000FF')
          axarr.plot([-1,98], [x+yplus[iiter]*DI for x in [range_control[0], range_control[0]]], 'k-',color='#AADDFF')
          axarr.plot([-1,98], [x+yplus[iiter]*DI for x in [range_control[1], range_control[1]]], 'k-',color='#AADDFF')
        if iter >= 0:
          thisCoeff = coeffCoeffs[iter][0]*theseCoeffs[iallmutval] + coeffCoeffs[iter][0]*(1.0 - 0.5*theseCoeffs[iallmutval])
        else:
          thisCoeff = 0
        print "iter="+str(iter)+", thisCoeff="+str(thisCoeff)
          
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
            if mutvar.find('off') > -1 or mutvar.find('ehcn') > -1:
              newVal =  [x+mutvals*thisCoeff for x in defVals[mutvar]]
              if mutvals >= 0 and kmutvar==0:
                mutText = mutText + "+" + str(mutvals) +" mV"
              elif kmutvar==0:
                mutText = mutText  + str(mutvals) +" mV"
            else:
              newVal = [x*(mutvals**thisCoeff) for x in defVals[mutvar]]
              if kmutvar==0:
                mutText = mutText + "*" + str(mutvals)
            if kmutvar < len(mutvars)-1:
              mutText = mutText + ", "
        print mutText

        axarr.plot([xplus[iiter]+counter,xplus[iiter]+counter], [x+yplus[iiter]*DI for x in thisRange], styles[iter],color=cols[iter])
        #if igene==0 and imut==0 and iallmutval==0 and iiter == 0:
        axarr.set_xlim([-4,150])
        axarr.set_xticks([])
        #axarr.set_ylim([0.000185+0.000014,0.00039+0.000014])
        #axarr.set_ylim([0.000105+0.000014,0.00039+0.000014])
        #axarr.set_yticks([min(Casoma_control),max(Casoma_control),min(Casoma_control)+DI,max(Casoma_control)+DI,min(Casoma_control)+2*DI,max(Casoma_control)+2*DI,0.000275+2*DI,0.0003+2*DI,0.000325+2*DI])
        for tick in axarr.yaxis.get_major_ticks():
          tick.label.set_fontsize(fs)
        labels = ['','','','','','','%0.3f' % 0.275,'%0.3f' % 0.3,'%0.3f' % 0.325]
        #for itick in range(0,3):
        #  labels[itick*2] = '%0.3f' % (1000*min(Casoma_control))
        #  labels[itick*2+1] = '%0.3f' % (1000*max(Casoma_control))
        #axarr.set_yticklabels(labels)
    axarr.text(counter-0.5*cumprodnVals[len(MT[igene][imut])-1]+labelxplus[labelcounter], 0.000225+0.000001*labelyplus[labelcounter], genelabels[labelcounter], {},rotation=90,fontsize=fs)

    counter = counter+1
   counter = counter+1
  axarr.text(-2.5, 0.000229, 'control', {},rotation=90,fontsize=fs)
  y1=0.00021+0.000014
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000002,y1+0.000004],'k-');
  thisline[0].set_clip_on(False)
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000000,y1+0.000002],'k-');
  thisline[0].set_clip_on(False)
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000001,y1+0.000003],'k-',color='#FFFFFF',zorder=100,linewidth=2.0);
  thisline[0].set_clip_on(False)
  y1=0.000245+0.000014
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000002,y1+0.000004],'k-');
  thisline[0].set_clip_on(False)
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000000,y1+0.000002],'k-');
  thisline[0].set_clip_on(False)
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000001,y1+0.000003],'k-',color='#FFFFFF',zorder=100,linewidth=2.0);
  thisline[0].set_clip_on(False)
  y1=0.000275+0.000014
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000002,y1+0.000004],'k-');
  thisline[0].set_clip_on(False)
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000000,y1+0.000002],'k-');
  thisline[0].set_clip_on(False)
  thisline = axarr.plot([-4.6,-3.4],[y1+0.000001,y1+0.000003],'k-',color='#FFFFFF',zorder=100,linewidth=2.0);
  thisline[0].set_clip_on(False)
  axarr.set_ylabel("PPI window span")
  
  f.savefig("ppiranges_cs"+str(icell)+".eps")
  
  
