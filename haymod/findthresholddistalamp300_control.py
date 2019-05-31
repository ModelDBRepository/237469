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

random.seed(1)

v0 = -80
ca0 = 0.0001
proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
fs = 8
maxSynsPerSeg = 100
Nsyns = 3000
maxLens = [1300,1185]

gmaxes = []

unpicklefile = open('synlocs300.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Nsyns = unpickledlist[0]
maxSynsPerSeg = unpickledlist[1]
synlocsAll = unpickledlist[3]

for icell in range(0,2):
  synlocs = synlocsAll[icell]
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"

  h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
objref cvode
cvode = new CVode()
cvode.active(1)
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
objref st1, st2
st1 = new IClamp(0.5)
st2 = new IClamp(0.5)
L5PC.soma st1
L5PC.soma st2
forsec L5PC.somatic {
}
forsec L5PC.apical {
}
L5PC.distribute_channels("apic","gIhbar_Ih",2,-0.8696,3.6161,0.0,1.0*2.0870,0.0002)
L5PC.distribute_channels("apic","gCa_HVAbar_Ca_HVA",3,1.0,0.1,685.0,885.0,1.0*0.000555)
L5PC.distribute_channels("apic","gCa_LVAstbar_Ca_LVAst",3,1.0,0.01,685.0,885.0,1.0*0.0187)
objref sl,ns,syn1,con1,isyn, tvec, syns["""+str(Nsyns)+"""]
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "distalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "distalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
access L5PC.apic[siteVec[0]]
L5PC.apic[siteVec[0]] {
  syn1 = new AlphaSynapse(siteVec[1])
  syn1.e = 0
  syn1.tau = 5
  syn1.onset = 10000 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
""")
  for istim in range(0,Nsyns):
    h("""
siteVec[0] = """+str(synlocs[istim][0])+"""
siteVec[1] = """+str(synlocs[istim][1])+"""
access L5PC.apic[siteVec[0]]
L5PC.apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000 + """+str(BACdt)+""" 
}
""")
  
  h("""
objref vsoma, vdend, recSite, vdend2, isoma, cadend, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
access L5PC.soma
cvode.record(&v(0.5),vsoma,tvec)
cvode.record(&cai(0.5),casoma,tvec)
sl = new List()
sl = L5PC.locateSites("apic","""+str(proximalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
print "proximalpoint gCa_HVA: ", L5PC.apic[siteVec[0]].gCa_HVAbar_Ca_HVA
print "proximalpoint gCa_LVA: ", L5PC.apic[siteVec[0]].gCa_LVAstbar_Ca_LVAst
access L5PC.apic[siteVec[0]]
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
recSite = new IClamp(siteVec[1])
recSite.amp = 0
L5PC.apic[siteVec[0]] {
        recSite
}
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  ITERS = 27
  if icell==0:
    nextgs = [0.000,0.001,0.0005]
  if icell==1:
    nextgs = [0.000,0.01,0.005]
  for iter in range(0,ITERS):
    tstop = 11000.0
    squareAmp = 0.0
    squareDur = 10.0
    epsp_gmax = nextgs[min(iter,2)]
    print epsp_gmax
    h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
""")
    for istim in range(0,Nsyns):
      h("syns["+str(istim)+"].gmax = "+str(epsp_gmax))
    h.init()
    h.run()

    times=np.array(h.tvec)
    Casoma=np.array(h.casoma)
    Cadend=np.array(h.cadend)
    Vsoma=np.array(h.vsoma)
    Vdend=np.array(h.vdend)
    nSpikes = len(mytools.spike_times(times,Vsoma,-35,100))
    if iter > 2 and nSpikes > 0:
      nextgs = [nextgs[0],nextgs[2],0.5*(nextgs[0]+nextgs[2])]
    if iter > 2 and nSpikes == 0:
      nextgs = [nextgs[2],nextgs[1],0.5*(nextgs[2]+nextgs[1])]
    print str(nSpikes)+", nextgs="+str(nextgs)

  gmaxes.append(nextgs[2])
picklelist = [gmaxes,Nsyns,maxSynsPerSeg]
file = open('thresholddistalamp300_control.sav', 'w')
pickle.dump(picklelist,file)
file.close()
