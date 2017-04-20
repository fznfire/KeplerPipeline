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
#from Detrender import SplineDetrend


def run(filepath,outputpath='',CampaignNumber=1,makelightcurve=True, find_transits=True,chunksize=300,cutoff_limit=2., method = 'Spitzer'):
  # Takes strings with the EPIC number of the star and input/outputpath. Campaign number is used to complete the correct filename as downloaded from MAST
  starname = re.search('[0-9]{9}',filepath).group(0)
  InitialTime = time.time()
  #handling case for spd vs lpd
  if "spd" in filepath:
    starname = starname+"_spd"

  outputfolder = os.path.join(outputpath,str(starname))
  if makelightcurve:

    # makes raw light curve from pixel file
    print time.time()-InitialTime

    if CampaignNumber==9 or CampaignNumber==10:
        AdditionalFilepath = filepath.replace("1_","2_")
        t1,f_t1,Xc1,Yc1 = StandardAperture(filepath,outputpath=outputpath,plot=False)
        t2,f_t2,Xc2,Yc2 = StandardAperture(filepath,outputpath=outputpath,plot=False)
        t = np.append(t1,t2)
        f_t = np.append(f_t1,f_t2)
        Xc = np.append(Xc1, Xc2)
        Yc = np.append(Yc1, Yc2)
    else:
        print "Alternative"
        t,f_t,Xc,Yc = StandardAperture(filepath,outputpath=outputpath,plot=False)

    print time.time()-InitialTime
    t1, f_t1 = t[:], f_t[:]/np.median(f_t)
    # removes outlying data points where thrusters are fired
    t,f_t,Xc,Yc = centroidfit.find_thruster_events(t,f_t,Xc,Yc,starname=starname,outputpath=outputfolder)

    #Detrend(t,f_t,Xc,Yc)
    t2, f_t2 = t[:], f_t[:]/np.median(f_t)

    print time.time()-InitialTime
    # now fit a polynomial to the data (inspired by Spitzer data reduction), ignore first data points which are not usually very high-quality
    if method == 'Spitzer':
        [t,f_t] = centroidfit.spitzer_fit(t[90:],f_t[90:],Xc[90:],Yc[90:],starname=starname,outputpath=outputpath,chunksize=chunksize)
    elif method == 'SFF':
        [t,f_t] = centroidfit.sff_fit(t[90:],f_t[90:],Xc[90:],Yc[90:],starname=starname,outputpath=outputpath,chunksize=chunksize)
    elif method == 'SplineDetrend':
        [t,f_t] = SplineDetrend(t[200:],f_t[200:],Xc[200:],Yc[200:])
    else:
        print 'No valid method given. Using Spitzer Polynomial Correction'
        [t,f_t] = centroidfit.spitzer_fit(t[90:],f_t[90:],Xc[90:],Yc[90:],starname=starname,outputpath=outputpath,chunksize=chunksize)

    print time.time()-InitialTime
    t3, f_t3 = t[:], f_t[:]

    [t,f_t] = centroidfit.clean_data(t,f_t) # do a bit of cleaning
    t4, f_t4 = t[:], f_t[:]

    pl.figure(figsize=(20,10))

    pl.subplot(221)
    std1 = np.std(f_t1)
    pl.plot(t1,f_t1, "ko", MarkerSize=2)
    pl.title("Raw Light Curve::" +str(len(t1)))
    pl.ylim(ymin=(1-3*std1),ymax=(1+3*std1))

    pl.subplot(222)
    std2 = np.std(f_t2)
    pl.plot(t2,f_t2,"ko", MarkerSize=2 )
    pl.title("Removed thruster events:: " +str(len(t2)))
    pl.ylim(ymin=(1-3*std2),ymax=(1+3*std2))

    pl.subplot(223)
    std3 = np.std(f_t3)
    pl.plot(t3,f_t3,"ko", MarkerSize=2)
    pl.title(method + " Correction:: " +str(len(t3)))
    pl.ylim(ymin=(-3*std3),ymax=(3*std3))

    pl.subplot(224)
    std4 = np.std(f_t4)
    pl.plot(t4,f_t4, "ko", MarkerSize=2)
    pl.title("Further data cleaning::" +str(len(t4)))
    pl.ylim(ymin=(-3*std4),ymax=(3*std4))

    pl.suptitle(starname)
    pl.savefig(outputfolder+"/CombinedSAP.png")
    pl.close('all')

    print time.time()-InitialTime
    outputlightcurvefolder = os.path.join(outputfolder,'lcs/') # one may want to put the LCs in a different folder e.g. to keep all together
    if not os.path.exists(outputlightcurvefolder):
      os.makedirs(outputlightcurvefolder)
    np.savetxt(os.path.join(outputlightcurvefolder, 'centroiddetrended_lightcurve_' + str(starname) + '.txt'),np.transpose([t,f_t]),header='Time, Flux')
  else:
    outputlightcurvefolder = os.path.join(os.path.join(outputpath,str(starname)),'lcs')
    [t,f_t] = np.loadtxt(os.path.join(outputlightcurvefolder, 'centroiddetrended_lightcurve_' + str(starname) + '.txt'),unpack=True,usecols=(0,1))

  print time.time()-InitialTime
  if find_transits:
    folded,f_t_folded,period,freqlist,powers = periodfinder.get_period(t,f_t,outputpath=outputpath,starname=starname,get_mandelagolmodel=True)
    np.savetxt(os.path.join(outputfolder, 'powerspectrum_' + str(starname) + '.txt'),np.transpose([freqlist,powers]),header='Frequencies, Powers')
    periodfinder.make_combo_figure(filepath,t,f_t,period,freqlist,powers,starname=starname,outputpath=outputpath)
  #pl.show() # comment out to keep things running for multiple stars

  pl.close('all')
  print time.time()-InitialTime
