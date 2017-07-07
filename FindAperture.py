import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as colors
import re
from astropy.io import fits
import os
from scipy.ndimage import convolve, measurements
from scipy.stats import laplace
import itertools
import operator


def ApertureOutline(StdAper,KepMag, AvgFlux, outputfolder, starname, XPos, YPos):
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
    pl.figure(figsize=(10,10))
    pl.imshow(AvgFlux,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
    pl.colorbar()
    pl.plot(segments[:,0]-0.5, segments[:,1]-0.5, color=(1,0,0,.5), linewidth=3)
    pl.plot(XPos, YPos, "ro")
    pl.plot(XCen,YCen,"go")
    pl.title(starname+":"+str(KepMag))
    pl.gca().invert_yaxis()
    pl.axis('equal')
    pl.tight_layout()
    pl.savefig(outputfolder+"/"+starname+".png", bbox_inches='tight')
    pl.close()

def Case1(AvgFlux,X,Y, StdCutOff=1):
    #returns the maximum aperture
    ExpectedFluxUnder = 1.1*np.median(AvgFlux)

    #find a standard Aperture
    AllAper = (AvgFlux>ExpectedFluxUnder)

    BkgAper = 1- AllAper
    BkgArray = AvgFlux[np.nonzero(BkgAper*AvgFlux)]
    BkgMedian = np.median(BkgArray)

    BkgStd = np.std(BkgArray)
    ExpectedFluxUnder = ExpectedFluxUnder+BkgStd*StdCutOff

    #return the biggest aperture
    TotalAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(TotalAper) # this numbers the different apertures distinctly
    area = measurements.sum(TotalAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    TotalAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (TotalAper >= np.max(TotalAper))*1  #backend process
    Distance = 3.5
    for i in range(1,np.max(TotalAper)+1):
        TempAper = (TotalAper==i)*1
        if np.sum(TempAper)>3:
            YCen, XCen = measurements.center_of_mass(TempAper*AvgFlux)
            TempDist = np.sqrt((X-XCen)**2+(Y-YCen)**2)
            if TempDist<Distance:
                Distance = TempDist
                StdAper = TempAper
    YCen, XCen = measurements.center_of_mass(StdAper*AvgFlux)
    return StdAper



def Case2(AvgFlux,X,Y,MedianTimes=1.5):
    #Use laplacian stencil to find all the stars in the scenes
    ExpectedFluxUnder = MedianTimes*np.median(AvgFlux)

    #find a standard Aperture
    AllAper = (AvgFlux>ExpectedFluxUnder)
    AllAper, num = measurements.label(AllAper) # this numbers the different apertures distinctly

    Distance  = 6.0 #Unacceptable distance
    for i in range(1,num+1):
        TempAper = (AllAper == i)
        YCen, XCen = measurements.center_of_mass(TempAper*AvgFlux)
        TempDist = np.sqrt((X-XCen)**2+(Y-YCen)**2)
        if TempDist<Distance:
            Distance = TempDist
            StdAper = TempAper
    if Distance>5.0:
      raise Exception('Failed to find the aperture')
    return StdAper

def Case3(AvgFlux,X,Y,StdCutOff=3.0):
    #Use laplacian stencil to find all the stars in the scenes
    Median = np.median(AvgFlux)

    #find a background
    BkgAper = (AvgFlux<Median)

    FluxArray = AvgFlux[np.nonzero(BkgAper*AvgFlux)]
    Std = np.std(FluxArray)

    AllAper = AvgFlux>(Std*StdCutOff+Median)
    AllAper, num = measurements.label(AllAper) # this numbers the different apertures distinctly

    Distance  = 5.0 #Unacceptable distance
    for i in range(1,num+1):
        TempAper = (AllAper == i)
        YCen, XCen = measurements.center_of_mass(TempAper*AvgFlux)
        TempDist = np.sqrt((X-XCen)**2+(Y-YCen)**2)
        if TempDist<Distance:
            Distance = TempDist
            StdAper = TempAper
    return StdAper




def Case5(AvgFlux, X,Y,Spacing=1):
    Xint, Yint = [int(round(X,0)), int(round(Y,0))]
    BkgVal = np.median(AvgFlux)

    XValues = np.arange(Xint-Spacing,Xint+Spacing+1,1)
    YValues = np.arange(Yint-Spacing,Yint+Spacing+1,1)

    ReferenceValue = 0
    #TODO Change tolerance value
    Dist_tolerance = 1.25
    for i,j in list(itertools.product(XValues,YValues)):
        try:
            TempAper2_2 = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
            TempAper2_2[i:i+2, j:j+2] = 1
            Signal = np.sum(AvgFlux*TempAper2_2)-4.*BkgVal
            Y_Cen, X_Cen = measurements.center_of_mass(AvgFlux*TempAper2_2)
            Distance = np.sqrt((X- X_Cen)**2+(Y- Y_Cen)**2)
            if Distanc2<Dist_tolerance:
                Value2_2 = Signal/np.sqrt(Signal+4.*BkgVal)
            else:
                Value2_2 = 1
        except:
            Value2_2 = 1

        try:
            TempAper2_3 = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
            TempAper2_3[i:i+2, j:j+3] = 1
            Signal = np.sum(AvgFlux*TempAper2_3)-6.*BkgVal
            Y_Cen, X_Cen = measurements.center_of_mass(AvgFlux*TempAper2_3)
            Distance = np.sqrt((X- X_Cen)**2+(Y- Y_Cen)**2)
            if Distance<Dist_tolerance:
                Value2_3 = Signal/np.sqrt(Signal+6.*BkgVal)
            else:
                Value2_3 = 2
        except:
            Value2_3 = 2

        try:
            TempAper3_2 = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
            TempAper3_2[i:i+3, j:j+2] = 1
            Value3_2 = np.sum(AvgFlux*TempAper3_2)-6.*BkgVal
            Y_Cen, X_Cen = measurements.center_of_mass(AvgFlux*TempAper3_2)
            Distance = np.sqrt((X- X_Cen)**2+(Y- Y_Cen)**2)
            if Distance<Dist_tolerance:
                Value3_2 = Signal/np.sqrt(Signal+6.*BkgVal)
            else:
                Value3_2 = 3
        except:
            Value3_2 = 3

        try:
            TempAper3_3 = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
            TempAper3_3[i:i+3, j:j+3] = 1
            Signal = np.sum(AvgFlux*TempAper3_3)-9.*BkgVal
            Y_Cen, X_Cen = measurements.center_of_mass(AvgFlux*TempAper3_3)
            Distance = np.sqrt((X- X_Cen)**2+(Y- Y_Cen)**2)
            if Distance<Dist_tolerance:
                Value3_3 = Signal/np.sqrt(Signal+9.*BkgVal)
            else:
                Value3_3 = 4
        except:
            Value3_3 = 4

        #star like shaped with five selection
        try:
            TempAper_Star = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
            TempAper_Star[i+1:i+2, j:j+3] = 1
            TempAper_Star[i:i+3, j+1:j+2] = 1
            Signal = np.sum(AvgFlux*TempAper_Star)-5.*BkgVal
            Y_Cen, X_Cen = measurements.center_of_mass(AvgFlux*TempAper_Star)
            Distance = np.sqrt((X- X_Cen)**2+(Y- Y_Cen)**2)
            if Distance<Dist_tolerance:
                Value_Star = Signal/np.sqrt(Signal+5.*BkgVal)
            else:
                Value_Star = 5
        except:
            Value_Star = 5

        #See which one is the best fit
        Values = np.array([Value2_2, Value2_3, Value3_2, Value3_3, Value_Star])
        MaxValue = max(Values)
        if MaxValue>ReferenceValue:
            ReferenceValue = MaxValue
            RefX, RefY = [i,j]
            TypeAperture = np.where(MaxValue == Values)[0][0]

    if ReferenceValue<7:
        #Find suitable 4 by 4 aperture based on distance
        print "-"*50
        print "Failed to find a good aperture"
        print "-"*50

        #Aperture just based on the distance
        RefX, RefY = [Xint, Yint]
        DistanceRef = 5
        for i,j in list(itertools.product(XValues,YValues)):
          try:
            TempAper2_2 = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
            TempAper2_2[i:i+2, j:j+2] = 1
            Y_Cen, X_Cen = measurements.center_of_mass(AvgFlux*TempAper2_2)
            Distance = np.sqrt((X- X_Cen)**2+(Y- Y_Cen)**2)
            if Distance<DistanceRef:
                RefX,RefY = [i,j]
                DistanceRef = Distance
          except:
            pass
        TypeAperture = 0


    Aperture = np.zeros(len(AvgFlux[0])*len(AvgFlux)).reshape(len(AvgFlux),len(AvgFlux[0]))
    i,j = [RefX, RefY]
    if TypeAperture == 0:
        Aperture[i:i+2,j:j+2] = 1
    elif TypeAperture == 1:
        Aperture[i:i+2, j:j+3] = 1
    elif TypeAperture == 2:
        Aperture[i:i+3, j:j+2] = 1
    elif TypeAperture == 3:
        Aperture[i:i+3, j:j+3] = 1
    elif TypeAperture == 4:
        Aperture[i+1:i+2, j:j+3] = 1
        Aperture[i:i+3, j+1:j+2] = 1
    else:
        raise Exception('Error finding a good aperture')

    return Aperture





def FindAperture(filepath='',outputpath='',SubFolder='',CampaignNum='1'):
    '''
    Centroid are calculated by center of mass function from scipy
    Background are fitting by spline.
    '''
    outputfolder = outputpath+'/'+SubFolder

    #extracting the starname
    starname = str(re.search('[0-9]{9}',filepath).group(0))

    #read the FITS file
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')


    #extract the vital information from the fits file
    TotalDate = np.array(FitsFile[1].data['Time'])
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    KepMag = FitsFile[0].header['Kepmag']
    X = FitsFile[2].header['CRPIX1']  - 1.0 #-1 to account for the fact indexing begins at 0 in python
    Y = FitsFile[2].header['CRPIX2'] - 1.0
    FitsFile.close()

    Index = np.where(Quality==0)
    GoodFlux = np.array(operator.itemgetter(*Index[0])(TotalFlux))
    TotalDate = np.array(operator.itemgetter(*Index[0])(TotalDate))

    if CampaignNum>8:
        AvgFlux = np.nanmedian(GoodFlux, axis=0)
        MedianValue = np.nanmedian(AvgFlux)-0.5
        AvgFlux[np.isnan(AvgFlux)] = MedianValue

        if "1_lpd" in filepath:
            starname = starname+"_1"
        else:

            starname = starname+"_2"
        if KepMag<16:
            StdAper = Case1(AvgFlux,X,Y,StdCutOff=4.5)*1 #use standard deviation of Background to cut it off
            #StdAper = Case2(AvgFlux,X,Y,MedianTimes=2.0) #Use the median value in the aperture
            #StdAper = Case3(AvgFlux,X,Y,StdCutOff=0.75) #Use the flux value of the star as cut off
            #StdAper = Case5(AvgFlux, X, Y, Spacing=3)*1 #
        else:
            StdAper = Case5(AvgFlux, X, Y, Spacing=2)*1 #
        ApertureOutline(StdAper,KepMag, AvgFlux, outputfolder, starname, X, Y)
        np.savetxt(outputfolder+"/"+starname+".txt",StdAper)

    else:
         #Two aperture method
         DateHalf = (max(TotalDate)+min(TotalDate))/2.0

         FirstHalf = TotalDate<DateHalf
         SecondHalf = TotalDate>DateHalf

         AvgFlux1 = np.nanmedian(GoodFlux[FirstHalf], axis=0)
         MedianValue1 = np.nanmedian(GoodFlux[FirstHalf])-0.5
         AvgFlux1[np.isnan(AvgFlux1)] = MedianValue1

         AvgFlux2 = np.nanmedian(GoodFlux[SecondHalf], axis=0)
         MedianValue2 = np.nanmedian(GoodFlux[SecondHalf])-0.5
         AvgFlux2[np.isnan(AvgFlux2)] = MedianValue2


         if KepMag<17:
             StdAper1 = Case1(AvgFlux1,X,Y,StdCutOff=2.0)*1 #use standard deviation of Background to cut it off
             StdAper2 = Case1(AvgFlux2,X,Y,StdCutOff=2.0)*1
             #StdAper = Case2(AvgFlux,X,Y,MedianTimes=1.25) #Use the median value in the aperture
             #StdAper = Case3(AvgFlux,X,Y,StdCutOff=1.0) #Use the flux value of the star as cut off
             #StdAper = Case4(AvgFlux, X, Y, Spacing=2)*1 #
         else:
             StdAper1 = Case5(AvgFlux1, X, Y, Spacing=2)*1 #
             StdAper2 = Case5(AvgFlux2, X, Y, Spacing=2)*1
         ApertureOutline(StdAper1,KepMag, AvgFlux1, outputfolder, starname+"_1", X, Y)
         ApertureOutline(StdAper2,KepMag, AvgFlux2, outputfolder, starname+"_2", X, Y)

         pl.figure(figsize=(16,16))
         pl.subplot(221)
         pl.imshow(AvgFlux1,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
         pl.title(str(np.sum(FirstHalf)))
         pl.colorbar()
         pl.subplot(222)
         pl.imshow(AvgFlux2,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
         pl.title(str(np.sum(SecondHalf)))
         pl.colorbar()
         pl.subplot(223)
         pl.imshow(StdAper1,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
         pl.title(str(np.sum(FirstHalf)))
         pl.subplot(224)
         pl.imshow(StdAper2,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
         pl.title(str(np.sum(SecondHalf)))
         pl.suptitle(starname+":"+str(KepMag))
         pl.savefig(outputfolder+"/"+starname+"_twohalf_diag.png", bbox_inches='tight')

         np.savetxt(outputfolder+"/"+starname+"_1.txt",StdAper1)
         np.savetxt(outputfolder+"/"+starname+"_2.txt",StdAper2)


    SummaryFile = open(outputfolder+".csv",'a')
    SummaryFile.write(starname+",1,0 \n")
    SummaryFile.close()

    print "Stage 3"
















def OldCase1(AvgFlux):
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

def OldCase2(AvgFlux):
    ExpectedFluxUnder = 2*np.median(AvgFlux)

    StdAper = (AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0

    return StdAper

def OldCase3(AvgFlux):
    ExpectedFluxUnder = 175
    StdAper = (AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0
    return StdAper
