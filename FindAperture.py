import numpy as np
import matplotlib.pyplot as pl
import re
from astropy.io import fits
import os
from scipy.ndimage import measurements

import matplotlib
matplotlib.use('Qt4Agg')


def ApertureOutline(StdAper,AvgFlux, outputfolder, starname):
    #find the outline and save the aperture in the relevant folder
    ver_seg = np.where(StdAper[:,1:] != StdAper[:,:-1])
    hor_seg = np.where(StdAper[1:,:] != StdAper[:-1,:])
    l = []

    for p in zip(*hor_seg):
        l.append((p[1], p[0]+1))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan,np.nan))

    for p in zip(*ver_seg):
        l.append((p[1]+1, p[0]))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan, np.nan))


    segments = np.array(l)
    pl.figure(figsize=(10,10))
    pl.imshow(AvgFlux,cmap='rainbow',interpolation='none')
    pl.colorbar()
    pl.plot(segments[:,0]-0.5, segments[:,1]-0.5, color=(1,0,0,.5), linewidth=3)
    pl.title(starname)
    pl.savefig(outputfolder+"/"+starname+".png")
    pl.close()

def Case1(AvgFlux):
    ExpectedFluxUnder = np.median(AvgFlux)

    #find a standard Aperture
    StdAper = (AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0

    #Finding the background aperture
    BkgAper = 1.0 - StdAper


    BkgFrame = (BkgAper*AvgFlux)
    BkgFrame = BkgFrame[np.nonzero(BkgFrame)]
    BkgStd = np.std(BkgFrame)
    BkgMedian = np.median(BkgFrame) #negative values for background are sometimes seen, which means that will be added to the flux values rather than subtracted
    Sigma = 5.0 #Usual value is 5
    CutoffLower = BkgMedian - Sigma*BkgStd #5 sigma cutoff for excluding really unusual pixel


    #New method
    BkgFrame = BkgFrame[np.nonzero((BkgFrame>CutoffLower)*1.0)]
    #BkgNewMean = np.median(BkgFrame)
    BkgNewMean = np.abs(np.median(BkgFrame))
    BkgNewStd = np.std(BkgFrame)

    Sigma = 2.0 ###Important for determining the aperture
    ExpectedFluxUnder = BkgNewMean+Sigma*BkgNewStd+15.0 #15.0 to consider the case where the background is really small


    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #
    return StdAper

def Case2(AvgFlux):
    ExpectedFluxUnder = 2*np.median(AvgFlux)

    StdAper = (AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0

    return StdAper

def Case3(AvgFlux):
    ExpectedFluxUnder = 175
    StdAper = (AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0
    return StdAper


def Case4(AvgFlux):

    AvgFlux[0,:] = 0
    AvgFlux[-1,:] = 0
    AvgFlux[:,0] = 0
    AvgFlux[:,-1] = 0
    
    return Case3(AvgFlux)

def Case5(AvgFlux):
    '''
    Manual
    '''
    pass


def FindAperture(filepath='',outputpath='',SubFolder=''):
    '''
    Centroid are calculated by center of mass function from scipy
    Background are fitting by spline.
    '''
    outputfolder = outputpath+'/'+SubFolder

    #extracting the starname
    starname = str(re.search('[0-9]{9}',filepath).group(0))

    #if short cadence data
    if "spd" in filepath:
      starname = starname+"_spd"

    #read the FITS file
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')



    #extract the vital information from the fits file
    TotalFlux = FitsFile[1].data['Flux']

    AvgFlux = np.nanmedian(TotalFlux, axis=0)
    AvgFlux[np.isnan(AvgFlux)] = 0

    #StdAper = Case1(AvgFlux) #Run twice
    #StdAper = Case2(AvgFlux) #  Twice the median
    #StdAper = Case3(AvgFlux) #Fixed value of 175
    StdAper = Case4(AvgFlux) #remove the boundary element
    #StdAper = Case5(AvgFlux) #Manual

    ApertureOutline(StdAper,AvgFlux, outputfolder, starname)
    np.savetxt(outputfolder+"/"+starname+".txt",StdAper)

    SummaryFile = open(outputfolder+".csv",'a')
    SummaryFile.write(starname+",1,0 \n")
    SummaryFile.close()
