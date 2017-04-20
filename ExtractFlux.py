from __future__ import division
import numpy as np
import matplotlib.pyplot as pl
import scipy.ndimage
from astropy.io import fits
from scipy.ndimage import convolve, measurements
import os
import re

def StandardAperture(filepath='',outputpath='',plot=False):
    starname = str(re.search('[0-9]{9}',filepath).group(0))
    if "spd" in filepath:
      starname = starname+"_spd"

    outputfolder = os.path.join(outputpath,starname)
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')
    TestPaths = [outputpath,outputpath+"/"+starname+"/"]
    for path in TestPaths:
        if not os.path.exists(path):
            os.system("mkdir %s" %(path))

    #extract the vital information from the fits file
    KeplerID = FitsFile[0].header['KEPLERID']
    TotalDate = FitsFile[1].data['Time']
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    RA = FitsFile[0].header['RA_OBJ']
    Dec = FitsFile[0].header['DEC_OBJ']
    KepMag = FitsFile[0].header['Kepmag']
    Xabs = FitsFile[2].header['CRVAL2P'] # X position of pixel on kepler spacecraft
    Yabs = FitsFile[2].header['CRVAL1P'] # Y position of pixel on kepler spacecraft

    #Doing median stack again average
    AvgFlux = np.nanmedian(TotalFlux, axis=0)
    AvgFlux[np.isnan(AvgFlux)] = 0

    #initiating array to collect values
    FluxArray = []
    DateArray = []
    XArray = []
    YArray = []
    BkgMeanArray = []

    #Find Background Value
    if "spd" in filepath:
      ExpectedFluxUnder = 100.0
    else:
      ExpectedFluxUnder = 400.0 #This value is arbitrarily choosen, and works for most stars

    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0

    #Check if the Aperture is too large
    if np.sum(StdAper)>(0.85*len(AvgFlux[0])*len(AvgFlux)):
        print "The traditional flux did not work."
        ExpectedFluxUnder = 2*np.median(AvgFlux)
        StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
        lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
        area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
        StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
        StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0
        if np.sum(StdAper)>(0.85*len(AvgFlux[0])*len(AvgFlux)):
          raise Exception('Error in finding Aperture. Maximum pixel value is ' + str(np.max(AvgFlux)))

    #Estimating the background
    BkgAper = 1.0 - StdAper

    #simple background subtraction
    BkgFrame = (BkgAper*AvgFlux)
    BkgFrame = BkgFrame[np.nonzero(BkgFrame)]
    BkgStd = np.std(BkgFrame)
    BkgMean = np.mean(BkgFrame)
    Sigma = 5
    CutoffUpper = BkgMean + Sigma*BkgStd #5 sigma cutoff for excluding really unusual pixel
    CutoffLower = BkgMean - Sigma*BkgStd #5 sigma cutoff for excluding really unusual pixel

    #redefining Background frame based on finding of the
    BkgFrame = BkgFrame[np.nonzero((BkgFrame<CutoffUpper)*(BkgFrame>CutoffLower)*1.0)]
    BkgNewMean = np.mean(BkgFrame)
    BkgNewStd = np.std(BkgFrame)

    Sigma = 2.5
    ExpectedFluxUnder = BkgNewMean+Sigma*BkgNewStd

    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #

    #Check if the Aperture is too large
    if np.sum(StdAper)>(0.85*len(AvgFlux[0])*len(AvgFlux)):
        raise Exception('Error in finding Aperture. Second Point. Maximum pixel value is ' + str(np.max(AvgFlux))+" and KepMag is " + str(KepMag))

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
    pl.figure()
    pl.imshow(AvgFlux,cmap='rainbow',interpolation='none')
    pl.colorbar()
    pl.plot(segments[:,0]-0.5, segments[:,1]-0.5, color=(1,0,0,.5), linewidth=3)
    pl.title("Aperture Selected")
    pl.savefig(outputfolder+"/Aperture.png")
    pl.close()


    for i in range(len(TotalFlux)):
        CurrentFrame = TotalFlux[i]
        CurrentFrame[np.isnan(CurrentFrame)] = 0.0 #converting all nan to zero
        Flux = CurrentFrame*StdAper

        FluxValue = np.sum(Flux)
        if FluxValue>0:
            YPos, XPos = measurements.center_of_mass(Flux)
            FluxArray.append(FluxValue)
            DateArray.append(TotalDate[i])
            #Uncomment later
            XArray.append(XPos)
            YArray.append(YPos)


    return np.array(DateArray), np.array(FluxArray), np.array(XArray), np.array(YArray)


def AdaptiveAperture(filepath='',outputpath='',campaign='',plot=False):
    starname = str(re.search('[0-9]{9}',filepath).group(0))
    if "spd" in filepath:
      starname = starname+"_spd"

    outputfolder = os.path.join(outputpath,starname)
    FitsFile = fits.open(filepath) #opening the fits file

    TestPaths = [outputpath,outputpath+"/"+starname+"/"]
    for path in TestPaths:
        if not os.path.exists(path):
            os.system("mkdir %s" %(path))

    #extract the vital information from the fits file
    KeplerID = FitsFile[0].header['KEPLERID']
    print "KEPLERID", KeplerID
    TotalDate = FitsFile[1].data['Time']
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    RA = FitsFile[0].header['RA_OBJ']
    Dec = FitsFile[0].header['DEC_OBJ']
    KepMag = FitsFile[0].header['Kepmag']


    print "Magnitude::",KepMag
    Xabs = FitsFile[2].header['CRVAL2P'] # X position of pixel on kepler spacecraft
    Yabs = FitsFile[2].header['CRVAL1P']

    #Doing median stack again average
    AvgFlux = np.nanmedian(TotalFlux, axis=0)
    AvgFlux[np.isnan(AvgFlux)] = 0

    #initiating array to collect values
    FluxArray = []
    DateArray = []
    XArray = []
    YArray = []
    BkgMeanArray = []



    #Find Background Value
    if "spd" in filepath:
      ExpectedFluxUnder = 100.0
    else:
      ExpectedFluxUnder = 400.0

    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0

    #Check if the Aperture is too large
    if np.sum(StdAper)>(0.85*len(AvgFlux[0])*len(AvgFlux)):
        raise Exception('Error in finding Aperture. Maximum pixel value is ' + str(np.max(AvgFlux)))

    #Estimating the background
    BkgAper = 1.0 - StdAper

    #simple background subtraction
    BkgFrame = (BkgAper*AvgFlux)
    BkgFrame = BkgFrame[np.nonzero(BkgFrame)]
    BkgStd = np.std(BkgFrame)
    BkgMean = np.mean(BkgFrame)
    CutoffUpper = BkgMean + 5*BkgStd #5 sigma cutoff for excluding really unusual pixel
    CutoffLower = BkgMean - 5*BkgStd #5 sigma cutoff for excluding really unusual pixel

    #redefining Background frame based on finding of the
    BkgFrame = BkgFrame[np.nonzero((BkgFrame<CutoffUpper)*(BkgFrame>CutoffLower)*1.0)]
    BkgNewMean = np.mean(BkgFrame)
    BkgNewStd = np.std(BkgFrame)

    ExpectedFluxUnder = BkgNewMean+1.75*BkgNewStd

    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1

    #Check if the Aperture is too large
    if np.sum(StdAper)>(0.85*len(AvgFlux[0])*len(AvgFlux)):
        raise Exception('Error in finding Aperture. Maximum pixel value is ' + str(np.max(AvgFlux)))


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

    pl.figure()
    pl.imshow(AvgFlux,cmap='gray',interpolation="none")
    pl.plot(segments[:,0]-0.5, segments[:,1]-0.5, color=(1,0,0,.5), linewidth=3)
    pl.title("Aperture Selected")
    pl.savefig(outputfolder+"/Aperture.png")

    for i in range(len(TotalFlux)):
        CurrentFrame = TotalFlux[i]
        CurrentFrame[np.isnan(CurrentFrame)] = 0.0 #converting all nan to zero
        Flux = CurrentFrame*StdAper
        XPos, YPos = measurements.center_of_mass(Flux)

        #Estimating the background
        BkgAper = 1.0 - StdAper

        #simple background subtraction
        BkgFrame = (BkgAper*CurrentFrame)
        BkgFrame = BkgFrame[np.nonzero(BkgFrame)]
        BkgStd = np.std(BkgFrame)
        BkgMean = np.mean(BkgFrame)
        CutoffUpper = BkgMean + 5*BkgStd #5 sigma cutoff for excluding really unusual pixel
        CutoffLower = BkgMean - 5*BkgStd #5 sigma cutoff for excluding really unusual pixel

        BkgFrame = BkgFrame[np.nonzero((BkgFrame<CutoffUpper)*(BkgFrame>CutoffLower)*1.0)]
        BkgNewMean = np.mean(BkgFrame) #Calculating the new mean based on the cutoff
        BkgMeanArray.append(BkgNewMean) #Delete this or save the diagram
        Background = np.sum(StdAper)*BkgNewMean
        FluxValue = np.sum(Flux) - Background
        if FluxValue>0:
            FluxArray.append(np.sum(FluxValue))
            DateArray.append(TotalDate[i])
            XArray.append(XPos)
            YArray.append(YPos)
    return np.array(DateArray), np.array(FluxArray), np.array(XArray), np.array(YArray)


def BackupApertureMethod(filepath='',outputpath='',campaign='',plot=False):
    starname = re.search('[0-9]{9}',filepath).group(0)
    outputfolder = os.path.join(outputpath,str(starname))
    FitsFile = fits.open(filepath) #opening the fits file

    #create outputpath if it does not already exist
    TestPaths = [outputpath,outputpath+"/"+starname+"/"]
    for path in TestPaths:
        if not os.path.exists(path):
            os.system("mkdir %s" %(path))


    LaplacianStencil = np.array([[0, 1, 0],[1, -4, 1], [0, 1, 0]])
    KeplerID = FitsFile[0].header['KEPLERID'] #error with the fix it later
    print "KEPLERID", KeplerID
    TotalDate = FitsFile[1].data['Time']

    TotalFlux = FitsFile[1].data['Flux']
    ArrayLength = len(TotalFlux)
    Quality = FitsFile[1].data['Quality']
    RA = FitsFile[0].header['RA_OBJ']
    Dec = FitsFile[0].header['DEC_OBJ']
    KepMag = FitsFile[0].header['Kepmag']
    Xabs = FitsFile[2].header['CRVAL2P'] # X position of pixel on kepler spacecraft
    Yabs = FitsFile[2].header['CRVAL1P']

    SumFlux = np.nansum(TotalFlux, axis=0)
    SumFlux[np.isnan(SumFlux)] = 0
    AvgFlux = SumFlux/len(TotalFlux)
    FluxArray = []
    BkgDateArray = []
    FluxDateArray = []
    MeanBkgArray= []
    ApertureArray = []
    XArray = []
    YArray = []
    BackgroundArray = []
    ExpectedFlux = 400.0
    StdAper = AvgFlux<ExpectedFlux
    StdAper = 1.0*(AvgFlux>ExpectedFlux)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1



    print "Starting Loop"
    for i in range(len(TotalFlux)):
        Value = np.nansum((1.0-StdAper)*TotalFlux[i])
        if Value != 0:
            BkgDateArray.append(TotalDate[i])
            BackgroundArray.append(Value/np.sum(1.0-StdAper))

    coeff = np.polyfit(BkgDateArray,BackgroundArray,3)

    x = np.linspace(min(BkgDateArray),max(BkgDateArray),1000)
    y = np.polyval(coeff,x)

    '''
    pl.figure()
    pl.plot(BkgDateArray,BackgroundArray,"k+")
    pl.plot(x,y,"r--")
    pl.xlabel("Julian Date")
    pl.ylabel("Background Pixel Value")
    pl.title("Background vs DateArray")
    pl.show()
    '''

    for i in range(len(TotalFlux)):
        CurrentDate = TotalDate[i]
        ExpectedBkg = np.polyval(coeff,CurrentDate)
        CurrentFrame = TotalFlux[i]
        StdAper = 1.0*(CurrentFrame>ExpectedBkg+50.0)
        lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
        area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
        StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
        StdAper = (StdAper >= np.max(StdAper))*1
        Flux = CurrentFrame*StdAper
        XPos, YPos = measurements.center_of_mass(Flux)
        Value = np.nansum(StdAper*TotalFlux[i]) - np.sum(StdAper)*ExpectedBkg
        if Value > 0.0:
            StdAper*TotalFlux[i]
            FluxArray.append(Value)
            FluxDateArray.append(CurrentDate)
            XArray.append(XPos)
            YArray.append(YPos)

    return np.array(FluxDateArray), np.array(FluxArray), np.array(XArray), np.array(YArray) #need to return the array
