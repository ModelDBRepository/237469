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

proximalpoint = 400
distalpoint = 620
BACdt = 5.0
fs = 8
maxLens = [1300,1185]

maxSynsPerSeg = 50
Nsyns = 3000

synlocsAll = []

for icell in range(0,2):
  morphology_file = "morphologies/cell"+str(icell+1)+".asc"
  biophys_file = "models/L5PCbiophys3.hoc"
  template_file = "models/L5PCtemplate.hoc"
  synlocs = []

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
objref vsoma, vdend, recSite, vdend2, isoma, cadend, casoma
vsoma = new Vector()
casoma = new Vector()
vdend = new Vector()
cadend = new Vector()
vdend2 = new Vector()
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
cvode.record(&v(siteVec[1]),vdend,tvec)
cvode.record(&cai(siteVec[1]),cadend,tvec)
recSite = new IClamp(siteVec[1])
recSite.amp = 0
L5PC.apic[siteVec[0]] {
        recSite
}
L5PC.apic[siteVec[0]] {
  syn1 = new AlphaSynapse(siteVec[1])
  syn1.e = 0
  syn1.tau = 5
  syn1.onset = 10000 + """+str(BACdt)+""" 
  cvode.record(&syn1.i,isyn,tvec)
}
""")
  synsInSegs = [0]*len(h.L5PC.apic)
  for istim in range(0,Nsyns):
    myiseg = -1
    while myiseg == -1:
      x = 300.0+(maxLens[icell]-300)*random.random()
      h("""sl = L5PC.locateSites("apic","""+str(x)+""")
Nsegs_x = sl.count()
""")
      iseg = random.randint(0,h.Nsegs_x-1)
      h("iseg = sl.o["+str(iseg)+"].x[0]")
      if synsInSegs[int(h.iseg)] < maxSynsPerSeg:
        myiseg = int(h.iseg)
        break
      print "istim = "+str(istim)+", x = "+str(x)+", continue searching for iseg..."
    synsInSegs[myiseg] = synsInSegs[myiseg] + 1
    h("""
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
access L5PC.apic[siteVec[0]]
L5PC.apic[siteVec[0]] {
  syns["""+str(istim)+"""] = new AlphaSynapse(siteVec[1])
  syns["""+str(istim)+"""].e = 0
  syns["""+str(istim)+"""].tau = 5
  syns["""+str(istim)+"""].onset = 10000 + """+str(BACdt)+""" 
}
""")
    synlocs.append([h.siteVec[0],h.siteVec[1]])
  
  h("""
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
access L5PC.soma
isoma = new Vector()
cvode.record(&st1.i,isoma,tvec)
""")

  synlocsAll.append(synlocs[:])
picklelist = [Nsyns,maxSynsPerSeg,maxLens,synlocsAll]
file = open('synlocs300.sav', 'w')
pickle.dump(picklelist,file)
file.close()
