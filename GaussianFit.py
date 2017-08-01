import numpy as np
from scipy.ndimage import convolve, measurements

from astropy.modeling import models, fitting
from DerivativeCentroid import cntrd
from K2SFF import K2SFF_VJ14
import periodfinder
import os

from photutils import DAOStarFinder
import matplotlib.pyplot as pl
from matplotlib import colors
from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.psf import IntegratedGaussianPRF, DAOGroup, BasicPSFPhotometry
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.io import fits
import re

from scipy.interpolate import UnivariateSpline
from scipy.signal import gaussian
from scipy.ndimage import filters
#from photutils.aperture_core.PixelAperture import do_photometry



def ApertureOutline(StdAper,AvgFlux, X,Y):
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
    YCen, XCen = measurements.center_of_mass(StdAper*AvgFlux)
    pl.figure()
    pl.imshow(AvgFlux,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
    pl.colorbar()
    pl.plot(segments[:,0]-0.5, segments[:,1]-0.5, color=(1,0,0,.5), linewidth=3)
    pl.plot(X,Y, "r+", markersize=10)
    pl.plot(XCen,YCen, "g+", markersize=10)
    pl.title("Aperture Selected")
    pl.gca().invert_yaxis()
    pl.axis('equal')
    pl.tight_layout()
    pl.show()


def Gaussian2DFit(filepath='', SubFolder='', outputpath=''):
    #Read the fits file
    Campaign = re.search('c[0-9]{2}',filepath).group(0)
    Campaign = int(Campaign[1:])
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')

    KeplerID = FitsFile[0].header['KEPLERID']
    print "KEPLERID:", KeplerID
    TotalDate = FitsFile[1].data['Time']
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    RA = FitsFile[0].header['RA_OBJ']
    Dec = FitsFile[0].header['DEC_OBJ']
    KepMag = FitsFile[0].header['Kepmag']
    print "Kepler Magnitude:", KepMag
    X = FitsFile[2].header['CRPIX1']  - 1.0 #-1 to account for the fact indexing begins at 0 in python
    Y = FitsFile[2].header['CRPIX2'] - 1.0


    AvgFlux = np.nanmean(TotalFlux, axis=0)
    AvgFlux[np.isnan(AvgFlux)] = np.nanmedian(AvgFlux)

    #load up the mask
    if SubFolder:
        AperLocation = "Apertures/"+SubFolder+"/"
    else:
        AperLocation = "Apertures/"+"Campaign"+str(Campaign)+"/"

    StarAperName = str(KeplerID)
    if Campaign>8:
        if "1_" in filepath:
            StarAperName = StarAperName+"_1"
        else:
            StarAperName = StarAperName+"_2"
        StdAper = np.loadtxt(AperLocation+StarAperName+".txt")
    else:
        StdAper = np.loadtxt(AperLocation+StarAperName+".txt")

    Fitter = fitting.LevMarLSQFitter()
    GaussInit = models.Gaussian2D(amplitude=1000, x_mean=5, y_mean=5, x_stddev=5, y_stddev=5, theta=0)

    XLen = len(AvgFlux[0])
    YLen = len(AvgFlux)
    Xs = np.linspace(-XLen/2.0,XLen/2.0,XLen)
    Ys = np.linspace(-YLen/2.0,YLen/2.0,YLen)

    XX,YY = np.meshgrid(Xs, Ys)

    COM_Xs = [] #Center of Mass X values
    COM_Ys = [] #Center of Mass Y values

    Gau_Xs = [] #GaussianFit X values
    Gau_Ys = [] #GaussianFit Y values
    time = []
    flux = []
    BkgArray = []

    for i in range(len(TotalDate)):
        CurrentFrame = TotalFlux[i]
        BkgMedian = np.nanmedian(CurrentFrame)
        CurrentFrame[np.isnan(CurrentFrame)] = BkgMedian
        FluxValue = np.sum(StdAper*CurrentFrame)
        QualityPass = bool(Quality[i]==0 or Quality[i]==16384 or Quality[i]==128 or np.isnan(BkgMedian))
        if FluxValue>0 and QualityPass:
            Background = np.sum(StdAper)*BkgMedian
            StarFrame = CurrentFrame*StdAper
            YCen, XCen = measurements.center_of_mass(StarFrame)
            FittedModel = Fitter(GaussInit,XX,YY,CurrentFrame)
            time.append(TotalDate[i])
            flux.append(FluxValue)
            BkgArray.append(Background)
            COM_Xs.append(XCen)
            COM_Ys.append(YCen)
            Gau_Xs.append(FittedModel.x_mean.value)
            Gau_Ys.append(FittedModel.y_mean.value)

    def moving_average(series, sigma=3):
        b = gaussian(39, sigma)
        average = filters.convolve1d(series, b/b.sum())
        var = filters.convolve1d(np.power(series-average,2), b/b.sum())
        return average, var

    _, var = moving_average(BkgArray)

    spl =  UnivariateSpline(time, BkgArray, w=1/np.sqrt(var))

    flux = np.array(flux)
    flux = flux - spl(time)

    #save the background vs time
    pl.figure()
    pl.plot(time, BkgArray, "ko")
    pl.plot(time, spl(time), "r-", lw=2)
    pl.savefig(outputpath+"/"+str(KeplerID)+"/BackgroundGauss.png")

    COM_Xs = np.array(COM_Xs)
    COM_Xs = COM_Xs - np.median(COM_Xs)

    COM_Ys = np.array(COM_Ys)
    COM_Ys = COM_Ys - np.median(COM_Ys)

    Gau_Xs = np.array(Gau_Xs) #GaussianFit X values
    Gau_Xs = Gau_Xs - np.median(Gau_Xs)

    Gau_Ys = np.array(Gau_Ys)
    Gau_Ys = np.array(Gau_Ys)



    return time, flux, Gau_Xs, Gau_Ys
