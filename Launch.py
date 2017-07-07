#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017
from __future__ import division
import matplotlib
matplotlib.use('Qt4Agg')


from run_pipeline import run
from glob import glob
import re

import warnings
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.

#Campaign = ...
#inputpath = '/Volumes/westep/prajwal/Campaign'+str(Campaign)+'/*.fits'
#SubFolder = 'Campaign'+str(Campaign)
#outputpath = 'Campaign'+str(Campaign)


inputpath = '/Volumes/westep/prajwal/ActiveStars/*.fits' #Campaign'+str(10)+'/*.fits'
SubFolder = 'ActiveStars'
outputpath = 'ActiveStars'

#inputpath = '/Volumes/westep/prajwal/PhaseCurves/*.fits' #Campaign'+str(10)+'/*.fits'
#SubFolder = 'PhaseCurves'
#outputpath = 'PhaseCurves'

filepaths = glob(inputpath)


i = 0
exc_list = []

while i < len(filepaths):
 EPIC_ID = re.search('[0-9]{9}',filepaths[i]).group(0)
 print "Currently running EPIC ID::", EPIC_ID
 try:
    #Don't run the spd for now
    if "spd" in filepaths[i]:
        pass
    else:
        chunksize = 100
        run(filepath=filepaths[i],outputpath=outputpath,chunksize=100,method ='SFF',SubFolder=SubFolder)
 except Exception as inst:
    print inst
    exc_list.append(inst)

 i = i + 1
 print str(i/len(filepaths)*100),"% completed"
print 'Done...'
