#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
import matplotlib
matplotlib.use('Qt4Agg')

from glob import glob
import re
import numpy as np
import os
from FindAperture import FindAperture
import time


import warnings
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.

inputpath = '/home/prajwal/Downloads/PhaseCurves/*.fits'
SubFolder = 'PhaseCurves'

outputpath = "Apertures/"
filepaths = glob(inputpath)
outputfolder = outputpath+SubFolder

TestPaths = [outputpath,outputfolder]

for path in TestPaths:
    if not os.path.exists(path):
        os.system("mkdir %s" %(path))

if not os.path.exists(outputpath+SubFolder+".csv"):
    FirstRun = True
    SummaryFile = open(outputpath+SubFolder+".csv",'w')
    SummaryFile.write("EPIC_ID,Run,Verified \n")
    SummaryFile.close()
else:
    FirstRun = False
    Status = np.loadtxt(outputpath+SubFolder+".csv",delimiter=',',skiprows=1,dtype=float)
    os.system("cp %s %s" %(outputpath+SubFolder+".csv",outputpath+SubFolder+"_"+str(time.time())+".csv"))
    SummaryFile = open(outputpath+SubFolder+".csv",'w')
    SummaryFile.write("EPIC_ID,Run,Verified \n")
    SummaryFile.close()



for filepath in filepaths:
  starname = str(re.search('[0-9]{9}',filepath).group(0))
  CampaignNum = re.search('c[0-9]{2}',filepath).group(0)
  CampaignNum = int(CampaignNum[1:])
  if FirstRun:
      print "Case 1 (First Time)::", starname
      if not("spd" in filepath):
        try:
            FindAperture(filepath=filepath,outputpath=outputpath,SubFolder=SubFolder, CampaignNum=CampaignNum)
        except Exception as inst:
            print inst
            SummaryFile = open(outputpath+SubFolder+".csv",'a')
            SummaryFile.write(starname+",0,0\n")
            SummaryFile.close()
      else:
          print "Skipping spd"
  else:
      '''See if the star image is validated properly'''
      Location =  np.where(Status[:,0]==int(starname))[0]
      if len(Location)==0:
          '''if the file is recently added'''
          print "New Addition"
          try:
              FindAperture(filepath=filepath,outputpath=outputpath,SubFolder=SubFolder, CampaignNum=CampaignNum)
          except Exception as inst:
              print inst
              SummaryFile = open(outputpath+SubFolder+".csv",'a')
              SummaryFile.write(starname+",0,0\n")
              SummaryFile.close()
      else:
          Location = Location[0]
          Verified = Status[Location][2]
          if not(Verified) and not("spd" in starname):
              print "Case 2 (Trying Again)::", starname
              try:
                  FindAperture(filepath=filepath,outputpath=outputpath,SubFolder=SubFolder, CampaignNum=CampaignNum)
              except Exception as inst:
                  print inst
                  SummaryFile = open(outputpath+SubFolder+".csv",'a')
                  SummaryFile.write(starname+",0,0\n")
                  SummaryFile.close()
          else:
                SummaryFile = open(outputpath+SubFolder+".csv",'a')
                SummaryFile.write(starname+",1,1\n")
                SummaryFile.close()



try:
    Status.close()
except:
    pass
