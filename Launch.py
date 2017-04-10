#Author: Van Eylen code modified by Prajwal Niraula
#Institution: Wesleyan University
#Last Updated: March 1, 2017

from __future__ import division
from run_pipeline import run
from glob import glob
import re

outputpath = 'Output/'
inputpath = '/Volumes/westep/prajwal/Campaign'


CampaignNumber = 1
inputpath = inputpath+str(CampaignNumber)+"/*.fits"

inputpath = "/home/pniraula/Downloads/Test"+"/*.fits" #For testing
filepaths = glob(inputpath)


'''
SelectedStars = ["201563164"] #limit the stars by EPIC ID

def LimitByEpic(Filepaths, EPICLists): #select by epic name
    TempFilePaths = []
    for epic in EPICLists:
        for filepath in Filepaths:
            if epic in filepath:
                TempFilePaths.append(filepath)
    return TempFilePaths

if len(SelectedStars)>0:
    filepaths = LimitByEpic(filepaths, SelectedStars)

#filepaths = ["/home/pniraula/Downloads/ktwo201577035-c01_lpd-targ.fits"]
#filepaths = ["/Volumes/westep/prajwal/Campaign5/ktwo211719918-c05_spd-targ.fits"]
#filepaths = ["/Volumes/westep/prajwal/Campaign8/ktwo220717512-c08_lpd-targ.fits"]
#filepaths = ["/Volumes/westep/prajwal/Campaign4/ktwo210659779-c04_lpd-targ.fits"] #change to spd
'''


filepaths = [filepaths[0]]


i = 0
exc_list = []

while i < len(filepaths):
  EPIC_ID = re.search('[0-9]{9}',filepaths[i]).group(0)
  print "Currently running EPIC ID::", EPIC_ID

  try:
    run(filepath=filepaths[i],outputpath=outputpath,makelightcurve=True,find_transits=True, method='SplineDetrend')
  except Exception as inst:
    print inst
    exc_list.append(inst)

  i = i + 1
  print str(i/len(filepaths)*100),"% completed"
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
