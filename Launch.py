#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017

import matplotlib
matplotlib.use('Qt4Agg')


from run_pipeline import run
from glob import glob
import re

import warnings
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.

outputpath = 'Temporary/'
CmpNum = 1
inputpath = '/Volumes/westep/prajwal/Campaign'+str(CmpNum)+'/*.fits'
SubFolder = 'Campaign'++str(CmpNum)


filepaths = glob(inputpath)


i = 0
exc_list = []

while i < len(filepaths):
  EPIC_ID = re.search('[0-9]{9}',filepaths[i]).group(0)
  print "Currently running EPIC ID::", EPIC_ID
  try:
    if "spd" in filepaths[i]:
        chunksize = 3000
    else:
        chunksize = 100
    run(filepath=filepaths[i],outputpath=outputpath,CNum=CmpNum,chunksize=100,method ='SFF',SubFolder=SubFolder)
  except Exception as inst:
    print inst
    exc_list.append(inst)

  i = i + 1
  print str(i/len(filepaths)*100),"% completed"
print 'Done...'
