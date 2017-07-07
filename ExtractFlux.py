from __future__ import division
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as colors
from scipy.interpolate import UnivariateSpline
from scipy.signal import gaussian
from scipy.ndimage import filters
from astropy.io import fits
from scipy.ndimage import convolve, measurements
import os
import re
import operator

from pixeltoflux import get_lightcurve

def ApertureOutline(StdAper,AvgFlux, outputfolder, X,Y):
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
    pl.imshow(AvgFlux,cmap='gray',norm=colors.PowerNorm(gamma=1./2.),interpolation='none')
    pl.colorbar()
    pl.plot(segments[:,0]-0.5, segments[:,1]-0.5, color=(1,0,0,.5), linewidth=3)
    pl.plot(X,Y, "ro")
    pl.title("Aperture Selected")
    pl.gca().invert_yaxis()
    pl.axis('equal')
    pl.tight_layout()
    pl.savefig(outputfolder+"/Aperture.png")
    pl.close()


def RawFluxDiagram(Quality,TotalDate,FluxArray,outputfolder):
    IndexArray = [Quality==0]
    FluxArray = np.array(FluxArray)
    for i in range(21):
        IndexArray.append(Quality==(2**i))

    ColorList = ["black","orange","blue","green", "cyan", "magenta", "red"]
    MarkerList = ['o','^','*']
    pl.figure(figsize=(16,8))
    pl.clf()
    for i in range(22):
        pl.plot(TotalDate[IndexArray[i]], FluxArray[IndexArray[i]],label=str(i)+":"+str(np.sum(IndexArray[i])), MarkerSize=3, marker = MarkerList[i%3], color=ColorList[i%7], linestyle='none')
    pl.legend(loc='best')
    pl.savefig(outputfolder+'/QualityIndex.png')
    pl.close()

    QualityIndex = Quality==0
    NoQualityIndex = Quality!=0

    pl.figure(figsize=(15,5))
    pl.plot(TotalDate[QualityIndex], FluxArray[QualityIndex], "r*",MarkerSize=2,label="Qualified")
    pl.plot(TotalDate[NoQualityIndex], FluxArray[NoQualityIndex], "ko",MarkerSize=2, label="Unqualified")
    pl.legend(loc='best')
    pl.tight_layout()
    pl.savefig(outputfolder+"/GoodDatavsBadData.png")
    pl.close()


def PredeterminedAperture(filepath='',outputpath='',plot=False, SubFolder=''):
    '''
    Centroid are calculated by center of mass function from scipy
    Background are fitting by spline.
    '''
    print "Running Predetermined Aperture"
    #extracting the starname
    starname = str(re.search('[0-9]{9}',filepath).group(0))
    Campaign = re.search('c[0-9]{2}',filepath).group(0)
    Campaign = int(Campaign[1:])

    print Campaign

    #if short cadence data
    if "spd" in filepath:
      starname = starname+"_spd"


    #constructing the output folder path
    outputfolder = os.path.join(outputpath,starname)

    #read the FITS file
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')

    #make the directory if the directory does not exist
    TestPaths = [outputpath,outputfolder]
    for path in TestPaths:
        if not os.path.exists(path):
            os.system("mkdir %s" %(path))

    #extract the vital information from the fits file
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



    #Writing the information to the summary
    if (Campaign>8) and ("2_" in filepath):
        pass
    else:
        RecordFile = open(outputpath+"/RunSummary.csv","a")
        RecordFile.write(starname+','+str(RA)+','+str(Dec)+','+str(KepMag)+',')
        RecordFile.close()



    #initiating array to collect values
    FluxArray = []
    X_Pos_Array = []
    Y_Pos_Array = []
    BkgArray = []
    FluxIndex = []


    if starname.endswith('_spd'):
        StarAperName = starname[:-4]
    else:
        StarAperName = starname

    if SubFolder:
        AperLocation = "Apertures/"+SubFolder+"/"
    else:
        AperLocation = "Apertures/"+"Campaign"+str(Campaign)+"/"

    if Campaign>8:
        if "1_" in filepath:
            StarAperName = StarAperName+"_1"
        else:
            StarAperName = StarAperName+"_2"
        StdAper = np.loadtxt(AperLocation+StarAperName+".txt")
    else:
        StdAper_1 = np.loadtxt(AperLocation+StarAperName+"_1.txt")
        StdAper_2 = np.loadtxt(AperLocation+StarAperName+"_2.txt")



    #copy the aperture file from the aperture file to the local file
    os.system("cp %s/*%s*.png %s" %(AperLocation,StarAperName,outputfolder))

    if Campaign>8:
        for i in range(len(TotalFlux)):
            CurrentFrame = TotalFlux[i]
            CurrentFrame[np.isnan(CurrentFrame)] = 0.0 #converting all nan to zero
            Flux = CurrentFrame*StdAper
            BkgMedian = np.median((1-StdAper)*CurrentFrame)

            #Getting the total value of the Flux
            Background = np.sum(StdAper)*BkgMedian
            FluxValue = np.sum(Flux)
            FluxArray.append(FluxValue)
            BkgArray.append(Background)
            if FluxValue>0 and Quality[i]==0:
                YPos, XPos = measurements.center_of_mass(Flux)
                X_Pos_Array.append(XPos)
                Y_Pos_Array.append(YPos)
                FluxIndex.append(True)
            else:
                FluxIndex.append(False)
    else:
        #for campaign 0-8
        DateHalf = (max(TotalDate)+min(TotalDate))/2.0
        for i in range(len(TotalFlux)):
            CurrentFrame = TotalFlux[i]
            CurrentFrame[np.isnan(CurrentFrame)] = 0.0 #converting all nan to zero
            CurrentDate = TotalDate[i]
            if CurrentDate<DateHalf:
                StdAper = StdAper_1
            else:
                StdAper = StdAper_2
            Flux = CurrentFrame*StdAper
            BkgMedian = np.median((1-StdAper)*CurrentFrame)

            #Getting the total value of the Flux
            Background = np.sum(StdAper)*BkgMedian
            FluxValue = np.sum(Flux)
            FluxArray.append(FluxValue)
            BkgArray.append(Background)
            if FluxValue>0 and Quality[i]==0:
                YPos, XPos = measurements.center_of_mass(Flux)
                X_Pos_Array.append(XPos)
                Y_Pos_Array.append(YPos)
                FluxIndex.append(True)
            else:
                FluxIndex.append(False)


    RawFluxDiagram(Quality,TotalDate,FluxArray,outputfolder)
    FluxIndex = np.array(FluxIndex)
    FluxArray = np.array(FluxArray)
    BkgArray = np.array(BkgArray)

    FluxArray = FluxArray[FluxIndex]
    BkgArray = BkgArray[FluxIndex]
    DateArray = TotalDate[FluxIndex]
    X_Pos_Array =  np.array(X_Pos_Array)
    Y_Pos_Array = np.array(Y_Pos_Array)

    def moving_average(series, sigma=3):
        b = gaussian(39, sigma)
        average = filters.convolve1d(series, b/b.sum())
        var = filters.convolve1d(np.power(series-average,2), b/b.sum())
        return average, var

    _, var = moving_average(BkgArray)
    if "spd" in starname:
        factor = 0.75
    else:
        factor = 0.9
    spl =  UnivariateSpline(DateArray, BkgArray, w=factor/np.sqrt(var))
    SplEstimatedBkg = spl(DateArray)



    #saving diagnostic plot
    pl.figure(figsize=(20,10))
    pl.subplot(2,1,1)
    pl.plot(DateArray, BkgArray,"k.", MarkerSize=2)
    pl.plot(DateArray, SplEstimatedBkg,"g-",lw=2)
    pl.xlabel('Time (days)')
    pl.ylabel('Flux Count')
    pl.title('Background')

    pl.subplot(2,1,2)
    pl.plot(DateArray, FluxArray,"ko",MarkerSize=2)
    pl.title('Flux Reading')
    pl.xlabel('Time (days)')
    pl.ylabel('Flux Count')

    pl.suptitle(str(KeplerID)+" Diagnostic")
    pl.savefig(outputfolder+"/Background.png")
    pl.close()

    return DateArray, FluxArray, X_Pos_Array, Y_Pos_Array







def ApertureOnTheRun(filepath='',outputpath='',plot=False,Campaign=None, QualityCheck=False):
    '''
    Centroid are calculated by center of mass function from scipy
    Background are fitting by spline.
    '''

    #extracting the starname
    starname = str(re.search('[0-9]{9}',filepath).group(0))

    #if short cadence data
    if "spd" in filepath:
      starname = starname+"_spd"

    #constructing the output folder path
    outputfolder = os.path.join(outputpath,starname)

    #read the FITS file
    try:
        FitsFile = fits.open(filepath,memmap=True) #opening the fits file
    except:
        raise Exception('Error opening the file')

    #make the directory if the directory does not exist
    TestPaths = [outputpath,outputfolder]
    for path in TestPaths:
        if not os.path.exists(path):
            os.system("mkdir %s" %(path))

    #extract the vital information from the fits file
    KeplerID = FitsFile[0].header['KEPLERID']
    print "KEPLERID:", KeplerID
    TotalDate = FitsFile[1].data['Time']
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    RA = FitsFile[0].header['RA_OBJ']
    Dec = FitsFile[0].header['DEC_OBJ']
    KepMag = FitsFile[0].header['Kepmag']
    print "Kepler Magnitude:", KepMag
    Xabs = FitsFile[2].header['CRVAL2P'] # X position of pixel on kepler spacecraft
    Yabs = FitsFile[2].header['CRVAL1P'] # Y position of pixel on kepler spacecraft

    #Writing the information to the summary
    if (Campaign==9 or Campaign==10) and ("2_" in filepath):
        pass
    else:
        RecordFile = open(outputpath+"/RunSummary.csv","a")
        RecordFile.write(starname+','+str(RA)+','+str(Dec)+','+str(KepMag)+',')
        RecordFile.close()


    #Taking the median of the frames
    AvgFlux = np.nanmedian(TotalFlux, axis=0)
    #take mean first and convert nans to zeros
    AvgFlux[np.isnan(AvgFlux)] = 0



    #initiating array to collect values
    FluxArray = []
    X_Pos_Array = []
    Y_Pos_Array = []
    BkgArray = []
    FluxIndex = []

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

    #Finding the background aperture
    BkgAper = 1.0 - StdAper

    #simple background subtraction
    BkgFrame = (BkgAper*AvgFlux)
    BkgFrame = BkgFrame[np.nonzero(BkgFrame)]
    BkgStd = np.std(BkgFrame)
    BkgMedian = np.median(BkgFrame) #negative values for background are sometimes seen, which means that will be added to the flux values rather than subtracted
    Sigma = 5.0 #Usual value is 5
    CutoffUpper = BkgMedian + Sigma*BkgStd #5 sigma cutoff for excluding really unusual pixel
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

    #Check if the Aperture is too large
    if np.sum(StdAper)>(0.85*len(AvgFlux[0])*len(AvgFlux)):
        raise Exception('Error in finding Aperture. Second Point. Maximum pixel value is ' + str(np.max(AvgFlux))+" and KepMag is " + str(KepMag))


    ####save the aperture
    ApertureOutline(StdAper, AvgFlux, outputfolder)



    for i in range(len(TotalFlux)):
        CurrentFrame = TotalFlux[i]
        CurrentFrame[np.isnan(CurrentFrame)] = 0.0 #converting all nan to zero
        Flux = CurrentFrame*StdAper

        BkgMedian = np.median((1-StdAper)*CurrentFrame)

        #Getting the total value of the Flux
        Background = np.sum(StdAper)*BkgMedian
        FluxValue = np.sum(Flux)

        FluxArray.append(FluxValue)
        BkgArray.append(Background)


        if FluxValue>0 and (Quality[i]==0 or Quality[i]==128):
            YPos, XPos = measurements.center_of_mass(Flux)
            X_Pos_Array.append(XPos)
            Y_Pos_Array.append(YPos)
            FluxIndex.append(True)
        else:
            FluxIndex.append(False)


    RawFluxDiagram(Quality,TotalDate,FluxArray,outputfolder)
    FluxIndex = np.array(FluxIndex)
    FluxArray = np.array(FluxArray)
    BkgArray = np.array(BkgArray)

    FluxArray = FluxArray[FluxIndex]
    BkgArray = BkgArray[FluxIndex]
    DateArray = TotalDate[FluxIndex]
    X_Pos_Array =  np.array(X_Pos_Array)
    Y_Pos_Array = np.array(Y_Pos_Array)

    def moving_average(series, sigma=3):
        b = gaussian(11, sigma)
        average = filters.convolve1d(series, b/b.sum())
        var = filters.convolve1d(np.power(series-average,2), b/b.sum())
        return average, var

    _, var = moving_average(BkgArray)
    if "spd" in starname:
        factor = 0.75
    else:
        factor = 0.9
    spl =  UnivariateSpline(DateArray, BkgArray, w=factor/np.sqrt(var))
    SplEstimatedBkg = spl(DateArray)



    #saving diagnositic plot
    pl.figure(figsize=(20,10))
    pl.subplot(2,1,1)
    pl.plot(DateArray, BkgArray,"k.", MarkerSize=2)
    pl.plot(DateArray, SplEstimatedBkg,"g-",lw=2)
    pl.xlabel('Time (days)')
    pl.ylabel('Flux Count')
    pl.title('Background')

    pl.subplot(2,1,2)
    pl.plot(DateArray, FluxArray,"ko",MarkerSize=2)
    pl.title('Flux Reading')
    pl.xlabel('Time (days)')
    pl.ylabel('Flux Count')

    pl.suptitle(str(KeplerID)+" Diagnostic")
    pl.savefig(outputfolder+"/Background.png")
    pl.close()

    return DateArray, FluxArray, X_Pos_Array, Y_Pos_Array




































def ModifiedStandardAperture(filepath='',outputpath='',plot=False, Campaign=1):
    '''Spline Fitted Background subtraction'''
    print "Running Modified Standard Aperture"
    starname = str(re.search('[0-9]{9}',filepath).group(0))

    if "spd" in filepath:
      starname = starname+"_spd"

    outputfolder = os.path.join(outputpath,starname)
    print "Output folder is :",outputfolder
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
    print "KEPLERID:", KeplerID
    TotalDate = FitsFile[1].data['Time']
    TotalFlux = FitsFile[1].data['Flux']
    Quality = FitsFile[1].data['Quality']
    RA = FitsFile[0].header['RA_OBJ']
    Dec = FitsFile[0].header['DEC_OBJ']
    KepMag = FitsFile[0].header['Kepmag']
    print "Kepler Magnitude:", KepMag
    Xabs = FitsFile[2].header['CRVAL2P'] # X position of pixel on kepler spacecraft
    Yabs = FitsFile[2].header['CRVAL1P'] # Y position of pixel on kepler spacecraft

    #Writing the information to the summary
    if (Campaign==9 or Campaign==10) and ("2_" in filepath):
        pass
    else:
        RecordFile = open(outputpath+"/RunSummary.csv","a")
        RecordFile.write(starname+','+str(RA)+','+str(Dec)+','+str(KepMag)+',')
        RecordFile.close()

    #Doing median stack again average
    AvgFlux = np.nanmean(TotalFlux, axis=0)
    AvgFlux[np.isnan(AvgFlux)] = 0.0

    #initiating array to collect values
    FluxArray = []
    DateArray = []
    XArray = []
    YArray = []
    BkgArray = []

    #Find Background Value
    if "spd" in filepath:
      ExpectedFluxUnder = 10.0 #Changed this value
    else:
      ExpectedFluxUnder = 10.0 #This value is arbitrarily choosen, and works for most stars

    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #make the standard aperture as 1.0


    #Estimating the background
    BkgAper = 1.0 - StdAper
    NewStdAper = BkgAper+StdAper

    #simple background subtraction
    BkgFrame = (BkgAper*AvgFlux)
    BkgFrame = BkgFrame[np.nonzero(BkgFrame)]
    BkgStd = np.std(BkgFrame)
    BkgMean = np.abs(np.mean(BkgFrame)) #negative values for background are sometimes seen
    Sigma = 5
    CutoffUpper = BkgMean + Sigma*BkgStd #5 sigma cutoff for excluding really unusual pixel
    CutoffLower = BkgMean - Sigma*BkgStd #5 sigma cutoff for excluding really unusual pixel


    #New method
    BkgFrame = BkgFrame[np.nonzero((BkgFrame>CutoffLower)*1.0)]
    #BkgNewMean = np.median(BkgFrame)
    BkgNewMean = np.abs(np.mean(BkgFrame))
    BkgNewStd = np.std(BkgFrame)

    Sigma = 2.5 ###Important for determining the aperture
    ExpectedFluxUnder = BkgNewMean+Sigma*BkgNewStd #15.0 to consider the case where the background is really small

    print "Std:",BkgNewStd
    print "Mean:", BkgNewMean
    #find a standard Aperture
    StdAper = 1.0*(AvgFlux>ExpectedFluxUnder)
    lw, num = measurements.label(StdAper) # this numbers the different apertures distinctly
    area = measurements.sum(StdAper, lw, index=np.arange(lw.max() + 1)) # this measures the size of the apertures
    StdAper = area[lw].astype(int) # this replaces the 1s by the size of the aperture
    StdAper = (StdAper >= np.max(StdAper))*1 #
    BkgAper = 1 - StdAper

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

    print "New StandardAperture"
    NewStdAper = BkgAper+StdAper

    for i in range(len(TotalFlux)):
        CurrentFrame = TotalFlux[i]
        CurrentFrame[np.isnan(CurrentFrame)] = 0.0 #converting all nan to zero
        Flux = CurrentFrame*NewStdAper

        BkgMedian = np.median((1-StdAper)*CurrentFrame)
        #Getting the total value of the Flux
        Background = np.sum(NewStdAper)*BkgMedian
        FluxValue = np.sum(Flux)#-Background


        if FluxValue>0:
            YPos, XPos = measurements.center_of_mass(Flux)
            FluxArray.append(FluxValue)
            DateArray.append(TotalDate[i])
            BkgArray.append(Background)
            XArray.append(XPos)
            YArray.append(YPos)

    def moving_average(series, sigma=3):
        b = gaussian(39, sigma)
        average = filters.convolve1d(series, b/b.sum())
        var = filters.convolve1d(np.power(series-average,2), b/b.sum())
        return average, var

    _, var = moving_average(BkgArray)
    if "spd" in starname:
        factor = 0.75
    else:
        factor = 0.9
    spl =  UnivariateSpline(DateArray, BkgArray, w=factor/np.sqrt(var))
    SplEstimatedBkg = spl(DateArray)

    FluxArray = np.array(FluxArray)-SplEstimatedBkg

    pl.figure(figsize=(20,10))
    pl.subplot(2,2,1)
    pl.plot(DateArray, BkgArray,"k.", MarkerSize=2)
    pl.plot(DateArray, SplEstimatedBkg,"g-",lw=2)
    pl.xlabel('Time (days)')
    pl.ylabel('Flux Count')
    pl.title('Background')

    pl.subplot(2,2,2)
    pl.plot(DateArray, FluxArray,"ko",MarkerSize=2)
    pl.title('Flux Reading')
    pl.xlabel('Time (days)')
    pl.ylabel('Flux Count')

    pl.subplot(2,2,3)
    pl.imshow(AvgFlux)
    pl.axis('equal')
    pl.title("Average Flux ")

    pl.subplot(2,2,4)
    pl.imshow(StdAper)
    pl.axis('equal')
    pl.title("Aperture")
    pl.suptitle(str(KeplerID)+" Diagnostic")
    pl.savefig(outputfolder+"/Background.png")
    pl.close()

    pl.figure()
    pl.plot(XArray, YArray,"g.", MarkerSize=2)
    pl.title('Drift in the field')
    pl.savefig(outputfolder+'/Drift.png')
    pl.close()

    return np.array(DateArray), FluxArray, np.array(XArray), np.array(YArray)
