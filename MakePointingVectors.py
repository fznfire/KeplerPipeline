#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017
from __future__ import division
import matplotlib
matplotlib.use('Qt4Agg')

from glob import glob
import re
import numpy as np
from astropy.io import fits

import warnings
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.

import centroidfit
from ExtractFlux import PredeterminedAperture
import periodfinder_mag
import matplotlib.pyplot as pl
from GaussianFit import Gaussian2DFit


inputpath = '/home/prajwal/Downloads/PointingStars/*.fits'
SubFolder = 'PointingStars'
outputpath = 'PointingVectors'

filepaths = glob(inputpath)

i = 0
exc_list = []


while i<len(filepaths):
  starname = re.search('[0-9]{9}',filepaths[i]).group(0)
  print starname
  Campaign = re.search('c[0-9]{2}',filepaths[i]).group(0)[1:]
  CampaignInt = int(Campaign)
  if CampaignInt>8 and Campaign<20:
      print starname
      Campaign = re.search('c[0-9]{3}',filepaths[i]).group(0)[1:]

  SaveName = Campaign+"_Pointing_"+starname+".txt"
  t,f_t,Xc,Yc = PredeterminedAperture(filepaths[i],outputpath=outputpath, SubFolder=SubFolder)
  np.savetxt((outputpath+"/"+SaveName), np.transpose([t,f_t,Xc, Yc]), delimiter=",", header="Time,Flux,Xc,Yc")

  SaveName = Campaign+"_GaussPointing_"+starname+".txt"
  t,f_t,Xc,Yc = Gaussian2DFit(filepath=filepaths[i],outputpath=outputpath, SubFolder=SubFolder)
  np.savetxt((outputpath+"/"+SaveName), np.transpose([t,f_t,Xc, Yc]), delimiter=",", header="Time,Flux,Xc,Yc")

  i=i+1
