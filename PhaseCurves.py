#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017

from __future__ import division
from run_pipeline import run
from glob import glob
import re

import warnings
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.

outputpath = 'PCurves/'




inputpath = "/home/pniraula/Downloads/PhaseCurveFitsFiles/*.fits"
inputpath = "/home/pniraula/Downloads/*.fits"

filepaths = glob(inputpath)


i = 0
exc_list = []

#filepaths = filepaths[0:1]

while i < len(filepaths):

  EPIC_ID = re.search('[0-9]{9}',filepaths[i]).group(0)
  print "Currently running EPIC ID::", EPIC_ID


  try:
    #if "spd" in filepaths[i]: #only for the short cadence data
        run(filepath=filepaths[i],outputpath=outputpath,makelightcurve=True,find_transits=True, method='SFF')
  except Exception as inst:
    print inst
    exc_list.append(inst)

  i = i + 1
if exc_list:
  print 'Module failed for some stars:'
  print exc_list

print 'Done...'



#Notes for the future
#starnames= ["201637175"]
#starnames = [ f[4:13] for f in listdir(inputpath) if isfile(join(inputpath,f)) ] # get all files in the folder, just grab the EPIC number of the filename

#print starnames
#starnames = [starnames[10]] #from campaign1
#starnames = [starnames[1]] #campaign 8
#starnames = ["201563164"] #WD 1145 from Campaign 1
#starnames =["220171396"] #good star from campaign 8
#validated target 201365699 from campaign 1 good for developing the framework
#starnames=["220171396"] #interesting star from campaign 8

#starnames = 212676650 #13 magnitude star from campaign 6
#starnames=['229228362']#17  many negative values in the frame,
#sum of the values comes out to be negative
#starnames=['220177283'] #18 magnitude star from Campaign 8 with easy aperture applying potential
#starnames = ['229228362'] #negative values in the grid c7 white dwarf
