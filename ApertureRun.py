#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University

from glob import glob
import re
from numpy import loadtxt
import os
from FindAperture import FindAperture
import time

outputpath = 'Apertures/'
SubFolder = "PhaseCurves"


inputpath = "/home/pniraula/Downloads/PhaseCurveFitsFiles/*.fits" #For testing

filepaths = glob(inputpath)

outputpath = 'Apertures/'
SubFolder = "PhaseCurves"
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
    Status = loadtxt(outputpath+SubFolder+".csv",delimiter=',',skiprows=1,dtype=int)
    os.system("cp %s %s" %(outputpath+SubFolder+".csv",outputpath+SubFolder+"_"+str(time.time())+".csv"))
    SummaryFile = open(outputpath+SubFolder+".csv",'w')
    SummaryFile.write("EPIC_ID,Run,Verified \n")
    SummaryFile.close()


i = 0
for filepath in filepaths:
  starname = str(re.search('[0-9]{9}',filepath).group(0))
  print "Currently running EPIC ID::", starname

  if FirstRun:
      print "Case 1 (First Time)::", starname
      try:
          FindAperture(filepath=filepath,outputpath=outputpath,SubFolder=SubFolder)
      except Exception as inst:
          print inst
          SummaryFile = open(outputpath+SubFolder+".csv",'a')
          SummaryFile.write(starname+",0,0\n")
          SummaryFile.close()
  elif(Status[i][2]==0):
      print "Case 2 (Trying Again)::", starname
      try:
          FindAperture(filepath=filepath,outputpath=outputpath,SubFolder=SubFolder)
      except Exception as inst:
          print inst
          SummaryFile = open(outputpath+SubFolder+".csv",'a')
          SummaryFile.write(starname+",0,0\n")
          SummaryFile.close()
  else:
      print "Case 3 (Already Verified)::", starname
      SummaryFile = open(outputpath+SubFolder+".csv",'a')
      SummaryFile.write(starname+",1,1\n")
      SummaryFile.close()
  i+=1



try:
    Status.close()
except:
    pass
