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
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.

#RUNID starts from 1000000 (1 million) and is increased by one every run
RUNID_File = open("RUNID.txt","r")
RUNID = RUNID_File.readlines()[-1][0:7]

if len(RUNID)<7:
    print "Error reading the RUNID. Update it manually"
RUNID_File.close()

Date = '_'.join(time.asctime().split(" ")[i] for i in [1,2,4,3])

#making directory if the folder does not already exist
if not(exists('Output')):
    system("mkdir Output")

#Making directory for a particular run
system("mkdir Output/%s" %(RUNID))

MethodName = "Spitzer"
PeriodFinder = "nf = 1e5 df = 1e-4"
SpecialNote = "Test Run"
ApertureSigma = "2.5"

CampaignNumber = [9]
CampaignStr = ','.join(str(i) for i in CampaignNumber)

outputpath = "Output/"+RUNID
RecordFile = open(outputpath+"/RunSummary.csv","w")
RecordFile.write("Date, %s \n" %(Date))
RecordFile.write("Campaigns, %s \n" %(CampaignStr))
RecordFile.write("Method, %s \n" %(MethodName))

#Manually update the variable
RecordFile.write("ApertureSigma, %s \n" %(ApertureSigma))
RecordFile.write("Period Finder values, %s \n" %(PeriodFinder))
RecordFile.write("Special Note, %s \n" %(SpecialNote)) #Have to manually change the value in extract flux
RecordFile.write("\n\n")
RecordFile.close()


#starting to run the campaig
for Campaign in CampaignNumber:

  #Header Data for each campaign
  RecordFile = open(outputpath+"/RunSummary.csv","a")
  RecordFile.write('CAMPAIGN,EPIC_ID,RA,DEC,MAG,SNR,PERIOD(days),RUNTIME(s),REDUCED,ERROR \n')
  RecordFile.close()

  inputpath = '/Volumes/westep/prajwal/Campaign'

  if Campaign == 9 or Campaign == 10:
      inputpath = inputpath+str(Campaign)+"/*1_*.fits"
  else:
      inputpath = inputpath+str(Campaign)+"/*.fits"

  filepaths = glob(inputpath)
  i = 0 #initiating counter for the files

  print "Starting Campaign::", str(Campaign)
  print "_"*75


  while i < len(filepaths)-5: #Need to modify this. Here here see this!

    EPIC_ID = re.search('[0-9]{9}',filepaths[i]).group(0) #extracting epic ID number
    print "Currently running EPIC ID::", filepaths[i], "  ",str(i+1), " out of ", str(len(filepaths))

    RecordFile = open(outputpath+"/RunSummary.csv","a")
    RecordFile.write(str(Campaign)+',')
    RecordFile.close()
    RunSuccess = "Failed"
    inst = ''
    try:
      run(filepath=filepaths[i],outputpath=outputpath,CNum=Campaign,makelightcurve=True,find_transits=True, method=MethodName)
      RecordFile = open(outputpath+"/RunSummary.csv","a")
      RecordFile.write('Success'+','+' ,\n')
      RecordFile.close()
    except Exception as inst:
      print str(inst) + "\n"
      RecordFile = open(outputpath+"/RunSummary.csv","a")
      RecordFile.write(',,,,,,'+'Failed'+','+str(inst)+'\n')
      RecordFile.close()

    i = i + 1


#Increase RUNID by one
RUNID_File = open("RUNID.txt","r+")
RUNID = int(RUNID_File.readlines()[-1])+1
RUNID_File.write(str(RUNID)+'\n')
RUNID_File.close()

print "Run completed"
