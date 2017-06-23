#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017

from __future__ import division
from run_pipeline import run
from glob import glob
import re
from numpy import loadtxt
import warnings
warnings.filterwarnings("ignore") #To suppress the warning. Comment this to see the range of warning.
from
outputpath = 'Apertures/'
SubFolder = "PhaseCurves"


inputpath = "/home/pniraula/Downloads/TestBackUp"+"/*.fits" #For testing

filepaths = glob(inputpath)

outputpath = 'Apertures/'
SubFolder = "PhaseCurves"
outputfolder = outputpath+SubFolder

TestPaths = [outputpath,outputfolder]


FirstRun = True
for path in TestPaths:
    F
    if not os.path.exists(path):
        os.system("mkdir %s" %(path))

if not os.path.exists(outputpath+SubFolder+".csv",'w'):
    SummaryFile = open(outputpath+SubFolder+".csv",'w')
    SummaryFile.write("EPIC_ID,Run,Verified")
    SummaryFile.close()
else:
    Status = loadtxt(outputpath+SubFolder+".csv",delimiter=',',skiprows=1,dtype=int)


i = 0


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
  print str(i/len(filepaths)*100),"% completed"
if exc_list:
  print 'Module failed for some stars:'
  print exc_list

try:
    Status = loadtxt(outputpath+SubFolder+".csv",delimiter=',',skiprows=1,dtype=int)
except:
    pass
