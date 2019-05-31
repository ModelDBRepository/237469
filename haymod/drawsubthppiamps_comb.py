#cp ../haymod3d/drawppiamps_relthr.py drawppiamps_relthr.py
from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *
import mytools
import pickle
import time
import sys
import random

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
theseMutValsAllAll = unpickledlist[2]

barxs = [1,2,3,-1,-2,0]
styles = ['b-','b-','b-','b-','b-','b-','b-','b-','b-','b-']
downstyles = ['b--','b--','b--','b--','b--','b--','b--','b--','b--','b--']
cols = ['#444444','#012345','#aa00aa','#bbaa00','#ee6600','#ff0000', '#009999','#772277','#00cc00']
col_control = '#2222ff'
coeffCoeffs = [[0.25,0],[0.125,0],[0.5,0],[0.5,1.0/3],[0.5,2.0/3],[0.5,1.0],[-0.25,0],[-0.125,0],[-0.5,0]]
lw = 0.5

CasomasAllAll = []
CadendsAllAll = []
VsomasAllAll = []
VdendsAllAll = []
timesAllAll = []

CasomaAmpsAllAll = []
CadendAmpsAllAll = []
VsomaAmpsAllAll = []
VdendAmpsAllAll = []

unpicklefile = open('ppispthrcoeff300_relthr_cs.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
theseCoeffsAllAll = unpickledlist[0]
gsAllAll = unpickledlist[1]
PPIdts = unpickledlist[2]
gCoeffsAllAll = unpickledlist[3]
gCoeffs_control = gCoeffsAllAll[0][0][0][0][5]
gs_control = gsAllAll[0][0][0][0][5]

for icell in range(0,1):
  gsAll = gsAllAll[icell]
  counter = -1

  unpicklefile = open('thresholddistalamp300_cs'+str(icell)+'_comb.sav', 'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  gsThisComb = unpickledlist[1]

  close("all")
  f, axarr = plt.subplots(1, 1)

  iters = [0, 2, 5, 6, 8]
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    if iter == 5:
      continue
    gs = gsThisComb[iiter]
    axarr.plot([0,500],[gs, gs], 'b--', color=cols[iter], linewidth=2.0, dashes=(1,2))
  axarr.plot([0,500],[gs_control, gs_control], 'b--', color=col_control, linewidth=2.0, dashes=(1,2))
  for iiter in range(0,len(iters)):
    iter = iters[iiter]
    if iter == 5:
      continue
    gs = gsThisComb[iiter]

    unpicklefile = open('ppispthrcoeff300_relthr_scaledonly_cs0_comb_iiter'+str(iiter)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    gCoeffsThisComb = unpickledlist[1]
    PPIdts = unpickledlist[2]
    gCoeffs = gCoeffsThisComb[0]

    axarr.plot(PPIdts,[1.1*gCoeffs[i]*gs for i in range(0,len(PPIdts))], color=cols[iter], linewidth=2.0)

  axarr.plot(PPIdts,[1.1*gCoeffs_control[i]*gs_control for i in range(0,len(PPIdts))], color=col_control, linewidth=2.0)
  axarr.set_xlabel('ISI (ms)',fontsize=15)
  axarr.set_ylabel('Threshold $g$ (nS)',fontsize=15)
  axarr.set_xticks([0, 100, 200, 300, 400])
  axarr.set_yticks([0, 0.00005, 0.0001, 0.00015, 0.0002, 0.00025])
  axarr.set_yticklabels(['0', '0.05', '0.1', '0.15', '0.2', '0.25'])

  f.savefig('subthppiamps_relthr_comb.eps')

  PPIdts = [40,50,60,70,80,100,500]
  currgs = [0.000025, 0.00005, 0.000075, 0.0001, 0.000125, 0.00015]

  iPPIs = [0,2,5]
  icurrs = [1,3,5]
      
  for tick in axarr.xaxis.get_major_ticks()+axarr.yaxis.get_major_ticks():
    tick.label.set_fontsize(15)
  axnew = f.add_axes([0.5, 0.45, 0.35, 0.4],axisbg='w')
  iys = [1,0,nan,3,4,2]
  for iiter in range(0,len(iters)+1):
    iter = -1
    if iiter < len(iters):
      iter = iters[iiter]
    if iter == 5:
      continue
    unpicklefile = open('subthppitest_cs0_comb_iiter'+str(iiter)+'.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    mycol = col_control
    if iter >= 0:
      mycol = cols[iter]
    times_all = unpickledlist[0]
    Vsoma_all = unpickledlist[1]
    for iiPPI in range(0,len(iPPIs)):
      for iicurr in range(0,len(icurrs)):
        times = times_all[iPPIs[iiPPI]][icurrs[iicurr]]
        Vsoma = Vsoma_all[iPPIs[iiPPI]][icurrs[iicurr]]
        istart = next((i for i,x in enumerate(times) if x > 10005))
        iend = next((i for i,x in enumerate(times) if x > 10155))
        axnew.plot([200*iiPPI+x for x in times[istart:iend]],[600*iicurr+95*iys[iiter]+x for x in Vsoma[istart:iend]],color=mycol,linewidth=2)
        if iiPPI == 0:
          axnew.text(9740,600*iicurr+95*3,str(1000*currgs[icurrs[iicurr]])+' nS:',fontsize=15)
      axnew.text(10000+iiPPI*200,-260,str(PPIdts[iPPIs[iiPPI]])+' ms',fontsize=15)
      axnew.plot([200*iiPPI+x for x in [10005,10005+PPIdts[iPPIs[iiPPI]]]],[-100]*2,'k.',markersize=2,mew=2)
  axnew.text(9850,-260,'ISI =',fontsize=15)
  axnew.set_xlim([9900,10600])
  axnew.set_ylim([-170,1640])
  axnew.plot([9910,9910+50],[50,50],'k-',linewidth=2)
  axnew.plot([9910,9910],[50,150],'k-',linewidth=2)

  axnew.set_axis_off()
  f.savefig('subthppiamps_relthr_comb_withV.eps')

