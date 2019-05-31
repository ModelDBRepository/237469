import pickle
from pylab import *
from os.path import exists
import time
import sys
import LFPy

filename = sys.argv[1]
if exists(filename):
  unpicklefile = open(filename,'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  dipolesE = unpickledlist[2]
  dt_data = 0.025 #ms
  dt_EEG = 0.5    #ms
  dtratio = int(dt_EEG/dt_data)

  dipolesE_lowres = zeros([int(11000/dt_EEG),3])

  for itime in range(0,dipolesE_lowres.shape[0]):
    dipolesE_lowres[itime,:] = mean(dipolesE[itime*dtratio:(itime+1)*dtratio,:],axis=0)

  radii = [79000., 80000., 85000., 90000.]
  sigmas = [0.3, 1.5, 0.015, 0.3]
  d_contact = 10.0
    
  somapos = array([0., 77500., 0.])
  eeg_coords_top = array([[0., radii[3]-d_contact, 0.]])
  four_sphere_top = LFPy.FourSphereVolumeConductor(radii, sigmas, eeg_coords_top)

  picklelist = []
  pot_db_4s_top = four_sphere_top.calc_potential(dipolesE_lowres, somapos)
  picklelist.append(array(pot_db_4s_top) * 1e9)

  file = open(filename[:-4]+'_EEG.sav','w')
  pickle.dump(picklelist,file)
  file.close()
  
