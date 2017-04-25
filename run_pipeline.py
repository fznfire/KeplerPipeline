'''
% Main pipeline file to generate K2 photometry starting from .fits pixel files downloaded from K2 MAST
% Author Vincent Van Eylen
% Contact vincent@phys.au.dk
% See Van Eylen et al. 2015 (ApJ) for details. Please reference this work if you found this code helpful!
'''

# general python files
import os
import matplotlib.pyplot as pl
import numpy as np
import time

# pipeline files
import pixeltoflux
import centroidfit
import periodfinder
import re
from ExtractFlux import AdaptiveAperture, StandardAperture

def run(filepath,outputpath='',CNum=1,makelightcurve=True, find_transits=True,chunksize=300,cutoff_limit=2., method ='Spitzer'):

  #Initializing the time
  InitialTime = time.time()

  # Takes strings with the EPIC number of the star and input/outputpath. Campaign number is used to complete the correct filename as downloaded from MAST
  starname = re.search('[0-9]{9}',filepath).group(0)


  #handling case for spd vs lpd
  if "spd" in filepath:
    starname = starname+"_spd"


  outputfolder = os.path.join(outputpath,str(starname))
  outputlightcurvefolder = os.path.join(outputfolder,'lcs/') # one may want to put the LCs in a different folder e.g. to keep all together
  if not os.path.exists(outputlightcurvefolder):
    os.makedirs(outputlightcurvefolder)

  if makelightcurve:

    # makes raw light curve from pixel file
    if CNum==9 or CNum==10:
        AdditionalFilepath = filepath.replace("1_","2_")
        t1,f_t1,Xc1,Yc1 = StandardAperture(filepath,outputpath=outputpath,plot=False,Campaign=CNum)
        t1,f_t1,Xc1,Yc1 = centroidfit.find_thruster_events(t1,f_t1,Xc1,Yc1,starname=starname,outputpath=outputfolder)

        t2,f_t2,Xc2,Yc2 = StandardAperture(AdditionalFilepath,outputpath=outputpath,plot=False,Campaign=CNum)
        t2,f_t2,Xc2,Yc2 = centroidfit.find_thruster_events(t2,f_t2,Xc2,Yc2,starname=starname,outputpath=outputfolder)

    else:
        t,f_t,Xc,Yc = StandardAperture(filepath,outputpath=outputpath,plot=False,Campaign=CNum)
        #Remove the thruster events
        t,f_t,Xc,Yc = centroidfit.find_thruster_events(t,f_t,Xc,Yc,starname=starname,outputpath=outputfolder)

    #save the raw Curve
    if CNum==9 or CNum==10:
        temp_t = np.append(t1,t2)
        temp_ft = np.append(t1,t2)
        np.savetxt(os.path.join(outputlightcurvefolder, 'RawFlux.txt'),np.transpose([temp_t,temp_ft]),header='Time, Flux')
        del temp_t, temp_ft
    else:
        np.savetxt(os.path.join(outputlightcurvefolder, 'RawFlux.txt'),np.transpose([t,f_t]),header='Time, Flux')
    # now fit a polynomial to the data (inspired by Spitzer data reduction), ignore first data points which are not usually very high-quality
    if method == 'Spitzer':
        if CNum==9 or CNum==10:
            [t1,f_t1] = centroidfit.spitzer_fit(t1,f_t1,Xc1,Yc1,starname=starname,outputpath=outputpath,chunksize=chunksize)
            [t2,f_t2] = centroidfit.spitzer_fit(t2,f_t2,Xc2,Yc2,starname=starname,outputpath=outputpath,chunksize=chunksize)
            t = np.append(t1,t2)
            f_t = np.append(f_t1,f_t2)
            del t1, t2, f_t1, f_t2
        elif CampaignNumber==1:
            [t,f_t] = centroidfit.spitzer_fit(t[90:],f_t[90:],Xc[90:],Yc[90:],starname=starname,outputpath=outputpath,chunksize=chunksize)
        else:
            [t,f_t] = centroidfit.spitzer_fit(t,f_t,Xc,Yc,starname=starname,outputpath=outputpath,chunksize=chunksize)

    elif method == 'SFF':
        if CNum==9 or CNum==10:
            [t1,f_t1] = centroidfit.sff_fit(t1,f_t1,Xc1,Yc1,starname=starname,outputpath=outputpath,chunksize=chunksize)
            [t2,f_t2] = centroidfit.sff_fit(t2,f_t2,Xc2,Yc2,starname=starname,outputpath=outputpath,chunksize=chunksize)
            t = np.append(t1,t2)
            f_t = np.append(f_t1,f_t2)
        elif CNum==1:
            [t,f_t] = centroidfit.sff_fit(t[90:],f_t[90:],Xc[90:],Yc[90:],starname=starname,outputpath=outputpath,chunksize=chunksize)
        else:
            [t,f_t] = centroidfit.sff_fit(t,f_t,Xc,Yc,starname=starname,outputpath=outputpath,chunksize=chunksize)
    else:
        raise Exception('No valid method given.')

    np.savetxt(os.path.join(outputlightcurvefolder, 'ReducedData.txt'),np.transpose([t,f_t]),header='Time, Flux')
    [t,f_t] = centroidfit.clean_data(t,f_t) # do a bit of cleaning
    np.savetxt(os.path.join(outputlightcurvefolder, 'Cleaned.txt'),np.transpose([t,f_t]),header='Time, Flux')

    np.savetxt(os.path.join(outputlightcurvefolder, 'CentroidDetrended.txt'),np.transpose([t,f_t]),header='Time, Flux')
  else:
    outputlightcurvefolder = os.path.join(os.path.join(outputpath,str(starname)),'lcs')
    [t,f_t] = np.loadtxt(os.path.join(outputlightcurvefolder, 'CentroidDetrended.txt'),unpack=True,usecols=(0,1))

  if find_transits:
    folded,f_t_folded,period,freqlist,powers = periodfinder.get_period(t,f_t,outputpath=outputpath,starname=starname,get_mandelagolmodel=False)
    np.savetxt(os.path.join(outputfolder, 'powerspectrum_' + str(starname) + '.txt'),np.transpose([freqlist,powers]),header='Frequencies, Powers')
    periodfinder.make_combo_figure(filepath,t,f_t,period,freqlist,powers,starname=starname,outputpath=outputpath)


  TimeTaken = time.time() - InitialTime
  RecordFile = open(outputpath+"/RunSummary.csv","a")
  TimeTakenStr = "%.2f" %TimeTaken
  print "Time Taken::",TimeTakenStr
  RecordFile.write(TimeTakenStr+',')
  RecordFile.close()
  pl.close('all')
