from neuron import h
import matplotlib
matplotlib.use('Agg')
import numpy
from pylab import *

def setparams(params,Nmc):
  global dists_apical, dists_basal

  keys = params.keys()
  lengthChanged = False
  for ikey in range(0,len(keys)):
    key = keys[ikey]
    if key[0:2] == "L_":
      lengthChanged = True
    underscoreind = key.rfind('_')
    section = key[underscoreind+1:len(key)]
    if section == "*":
      for i in range(0,Nmc):
        h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  forsec epnm.pc.gid2cell(i).all """+key[0:underscoreind]+""" = """+str(params[key])+"""
}
""")
      h("forall "+key[0:underscoreind]+" = "+str(params[key]))
    else:
      for i in range(0,Nmc):
        h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
  epnm.pc.gid2cell(i)."""+section+""" """+key[0:underscoreind]+""" = """+str(params[key])+"""
}
""")

  if lengthChanged:
    for i in range(0,Nmc):
      h("""
i = """+str(i)+"""
if (epnm.gid_exists(i)) {
    epnm.pc.gid2cell(i).soma diam = """+str(360.132/params['L_soma'])+"""
    epnm.pc.gid2cell(i).dend diam = """+str(2821.168/params['L_dend'])+"""
    epnm.pc.gid2cell(i).apic[0] diam = """+str(4244.628/params['L_apic[0]'])+"""
    epnm.pc.gid2cell(i).apic[1] diam = """+str(2442.848/params['L_apic[1]'])+"""

    lengthA = 0
    lengthB = 0
    forsec epnm.pc.gid2cell(i).apical {
      lengthA = lengthA + L
    }
    forsec epnm.pc.gid2cell(i).basal {
      lengthB = lengthB + L
    }
    epnm.pc.gid2cell(i).pA = lengthA/(lengthA + lengthB)
}
""")

    #No need for dists or dendritic recordings/stimuli?
    #dists_apical = []
    #dists_basal = []
    #for j in range(0,nrecsperseg):
    #  dists_apical.append(h.distance(xsperseg[j],sec=h.L5PC.apic[0]))
    #for j in range(0,nrecsperseg):
    #  dists_apical.append(h.distance(xsperseg[j],sec=h.L5PC.apic[1]))
    #for j in range(0,nrecsperseg):
    #  dists_basal.append(h.distance(xsperseg[j],sec=h.L5PC.dend))
    #if params['L_apic[0]'] > 200:
    #  for i in range(0,Nmc):
    #    h("""
    #i = """+str(i)+"""
    #if (epnm.gid_exists(i)) {
    #  epnm.pc.gid2cell(i)."""+section+""" """+key[0:underscoreind]+""" = """+str(params[key])+"""
    #}
    #""")
    #  h("L5PC.apic[0] st2.loc("+str(200.0/params['L_apic[0]'])+")")
    #elif params['L_apic[0]'] + params['L_apic[1]'] > 200:
    #  h("L5PC.apic[1] st2.loc("+str((200.0-params['L_apic[0]'])/params['L_apic[1]'])+")")
    #else:
    #  h("L5PC.apic[1] st2.loc(1.0)")
    #if params['L_apic[0]'] > 620:
    #  h("L5PC.apic[0] syn1.loc("+str(620.0/params['L_apic[0]'])+")")
    #  h("L5PC.apic[0] syni.loc("+str(620.0/params['L_apic[0]'])+")")
    #  h("L5PC.apic[0] cvode.record(&v("+str(620.0/params['L_apic[0]'])+"),vdend,tvec)")
    #  h("L5PC.apic[0] cvode.record(&cai("+str(620.0/params['L_apic[0]'])+"),cadend,tvec)")
    #elif params['L_apic[0]'] + params['L_apic[1]'] > 620:
    #  h("L5PC.apic[1] syn1.loc("+str((620.0-params['L_apic[0]'])/params['L_apic[1]'])+")")
    #  h("L5PC.apic[1] syni.loc("+str((620.0-params['L_apic[0]'])/params['L_apic[1]'])+")")
    #  h("L5PC.apic[1] cvode.record(&v("+str((620.0-params['L_apic[0]'])/params['L_apic[1]'])+"),vdend,tvec)")
    #  h("L5PC.apic[1] cvode.record(&cai("+str((620.0-params['L_apic[0]'])/params['L_apic[1]'])+"),cadend,tvec)")
    #else:
    #  h("L5PC.apic[1] syn1.loc(1.0)")
    #  h("L5PC.apic[1] syni.loc(1.0)")
    #  h("L5PC.apic[1] cvode.record(&v(1.0),vdend,tvec)")
    #  h("L5PC.apic[1] cvode.record(&cai(1.0),cadend,tvec)")
    #if params['L_apic[0]'] > 800:
    #  h("L5PC.apic[0] cvode.record(&v("+str(800.0/params['L_apic[0]'])+"),vdend2,tvec)")
    #  h("L5PC.apic[0] cvode.record(&cai("+str(800.0/params['L_apic[0]'])+"),cadend2,tvec)")
    #elif params['L_apic[0]'] + params['L_apic[1]'] > 800:
    #  h("L5PC.apic[1] cvode.record(&v("+str((800.0-params['L_apic[0]'])/params['L_apic[1]'])+"),vdend2,tvec)")
    #  h("L5PC.apic[1] cvode.record(&cai("+str((800.0-params['L_apic[0]'])/params['L_apic[1]'])+"),cadend2,tvec)")
    #else:
    #  h("L5PC.apic[1] cvode.record(&v(1.0),vdend2,tvec)")
    #  h("L5PC.apic[1] cvode.record(&cai(1.0),cadend2,tvec)")
