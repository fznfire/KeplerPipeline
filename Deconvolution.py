import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as colors
import re
from astropy.io import fits
import os
from scipy.ndimage import convolve, measurements
import scipy.stats as stats
import astropy
import glob
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
import operator
from lmfit import Parameter, minimize


def PSF_Model(params):
    A = params['A'] #Amplitude of the PSF function
    B = params['B'] #rotation term
    D = params['D'] #the sky background term
    return A*np.exp(-z - B*(x-Xc)*(y-Yc))+D


def Deconvolution(filepath):
    starname = str(re.search('[0-9]{9}',filepath).group(0))
    #read the FITS file
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')
    RA = FitsFile[0].header['RA_OBJ']
    DEC = FitsFile[0].header['DEC_OBJ']
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    X = FitsFile[2].header['CRPIX1']  - 1.0 #-1 to account for the fact indexing begins at 0 in python
    Y = FitsFile[2].header['CRPIX2'] - 1.0

    Index = np.where(Quality==0)
    GoodFlux = operator.itemgetter(*Index[0])(TotalFlux)

    GoodFluxAvg = np.nanmedian(GoodFlux, axis=0)
    print "Good data percentage::", float(len(Index[0]))/len(TotalFlux)*100.0
    AvgFlux = np.nanmedian(TotalFlux, axis=0)

    LaplacianStencil = np.array([[0, -1, 0],[-1, 4, -1], [0, -1, 0]])
    Laplacian = convolve(AvgFlux, LaplacianStencil)

    LapAperture = Laplacian>25
    NewAvgFlux = AvgFlux>np.nanmedian(AvgFlux)
    CombinedImage = AvgFlux*LapAperture

    LW, Num = measurements.label(CombinedImage)


    #Find you stars and their location
    NumStar = 0
    ApertureFitting = np.zeros(len(GoodFluxAvg)*len(GoodFluxAvg[0])).reshape(len(GoodFluxAvg),len(GoodFluxAvg[0]))

    Locations = []
    for i in range(1,Num+1):
        Aper  = (LW==i)
        YCen, XCen = measurements.center_of_mass(Aper*GoodFluxAvg)
        print X,Y
        print XCen,YCen
        Distance = np.sqrt((X-XCen)**2+(Y-YCen)**2)
        if Distance<6.0:
            ApertureFitting += Aper
            Locations.append([XCen,YCen])



    pl.figure()
    pl.imshow(ApertureFitting,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
    pl.plot(X,Y, "ro")
    pl.title(str(i))
    pl.show()

    #try to minimize the convolving function





    #

    #Median = np.median(np.nonzero(CombinedImage))
    #print "Median...", Median



    #Find Apertures for different stars


    #print CombinedImage

    #Fit the PSF function

    #NumStar =


    #TODO AvgFlux and good flux are different
    '''
    pl.figure()
    pl.subplot(2,2,1)
    pl.imshow(GoodFluxAvg)
    pl.colorbar()

    pl.subplot(2,2,2)
    pl.imshow(AvgFlux)
    pl.colorbar()

    pl.subplot(2,2,3)
    pl.imshow(AvgFlux-GoodFluxAvg)
    pl.colorbar()

    pl.show()
    '''





def Deconvolution1(filepath):
    starname = str(re.search('[0-9]{9}',filepath).group(0))
    #read the FITS file
    print "Stage 1"
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')
    RA = FitsFile[0].header['RA_OBJ']
    DEC = FitsFile[0].header['DEC_OBJ']


#Filepath = "/Volumes/westep/prajwal/Campaign1/ktwo201426232-c01_lpd-targ.fits"
#Filepath = "//home/pniraula/Downloads/PhaseCurveFitsFiles/ktwo201205469-c01_lpd-targ.fits"
Filepath = "/Volumes/westep/prajwal/Campaign1/ktwo201252863-c01_lpd-targ.fits"
Deconvolution(Filepath)
