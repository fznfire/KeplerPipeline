#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017

from __future__ import division
from run_pipeline import run
from glob import glob
from os import system
from os.path import exists
import time
import re

import warnings
warnings.filterwarnings("ignore")

CampaignNumber = [4,5,6,7]


RUNID = '_'.join(time.asctime().split(" ")[i] for i in [1,2,4,3])

#making directory if the folder does not already exist
if not(exists('Output')):
    system("mkdir Output")

#Making directory for a particular run
system("mkdir Output/%s" %(RUNID))
outputpath = "Output/"+RUNID
RecordFile = open(outputpath+"/RunSummary.txt","w")

#starting to run the campaign
for Campaign in CampaignNumber:

  inputpath = '/Volumes/westep/prajwal/Campaign'

  if Campaign == 9 or Campaign == 10:
      inputpath = inputpath+str(Campaign)+"/*1_*.fits"
  else:
      inputpath = inputpath+str(Campaign)+"/*.fits"

  filepaths = glob(inputpath)
  i = 0
  exc_list = []

  print "Starting Campaign::", str(Campaign)
  RecordFile.write(str(Campaign)+'\n')
  print "_"*75
  while i < len(filepaths):
    EPIC_ID = re.search('[0-9]{9}',filepaths[i]).group(0)
    print "Currently running EPIC ID::", filepaths[i], "  ",str(i+1), " out of ", str(len(filepaths))
    RecordFile.write(str(Campaign)+'::')
    RunSuccess = False
    inst = ''
    try:
      run(filepath=filepaths[i],outputpath=outputpath,CampaignNumber=Campaign,makelightcurve=True,find_transits=True, method='Spitzer')
      RunSuccess = True
    except Exception as inst:
      print str(inst) + "\n"
      exc_list.append(inst)
    RecordFile.write(str(EPIC_ID)+'::'+str(RunSuccess)+"  "+str(inst)+"\n")

    i = i + 1

if exc_list:
    print 'Module failed for some stars:'
    print exc_list

RecordFile.close()
print "Run completed"
