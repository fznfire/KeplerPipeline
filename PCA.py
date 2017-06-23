import matplotlib
matplotlib.use('Qt4Agg')

import numpy as np
import glob
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, KernelPCA
from ExtractFlux import StandardAperture

import centroidfit
import copy
import re
from scipy.interpolate import UnivariateSpline
from scipy.signal import gaussian
from scipy.ndimage import filters
from astropy.io import fits
from scipy.ndimage import convolve, measurements

from TransitFit import IndividualTransitFit, MaskTransit, FindPeriod
from pixeltoflux import get_lightcurve

#inputpath = "/home/pniraula/Downloads/Test"+"/*220403421*.fits" #Target with two apparent transits


inputpath = "/home/pniraula/Downloads/Test"+"/*211841249*lpd*.fits" #sinusoidal variation


filepaths = glob.glob(inputpath)
filepath = filepaths[0]
print filepath

starname = str(re.search('[0-9]{9}',filepath).group(0))
if "spd" in inputpath:
    starname+="_spd"

outputpath = "DeleteMe/"
outputfolder = outputpath+"/"+starname

'''Extract the file here'''

try:
    FitsFile = fits.open(filepath,memmap=True) #opening the fits file
except:
    raise Exception('Error opening the file')

dates = FitsFile[1].data['Time']
fluxes = FitsFile[1].data['Flux']

MedianFlux = np.nanmedian(fluxes,axis=0)
aperture = MedianFlux>400

plt.figure()
plt.subplot(211)
plt.imshow(MedianFlux)
plt.subplot(212)
plt.imshow(aperture)
plt.show()


t,f_t, X,Y = get_lightcurve(dates,fluxes,aperture,starname=starname,plot=True,outputfolder='',Xabs=0,Yabs=0)
t_new, f_new, X_new, Y_new = StandardAperture(filepath,outputpath=outputpath,plot=False,Campaign=5)

t_intersection = np.intersect1d(t,t_new)
Indices1 = np.sum([t==i for i in t_intersection], axis=0).astype(bool)
Indices2 = np.sum([t_new==i for i in t_intersection], axis=0).astype(bool)


X_1_ = X[np.array(Indices1)]
X_2_ = X_new[np.array(Indices2)]
Y_1_ = Y[np.array(Indices1)]
Y_2_ = Y_new[np.array(Indices2)]


plt.figure(figsize=(16,7))
plt.subplot(121)
plt.xlabel("Centroid Value")
plt.ylabel("Center of Mass")
plt.plot(X_1_, X_2_, "k.")
plt.title("X Values")
plt.subplot(122)
plt.plot(Y_1_, Y_2_, "k.")
plt.xlabel("Centroid Value")
plt.ylabel("Center of Mass")
plt.title("Y Values")
plt.suptitle(starname)
plt.show()


plt.figure()
plt.subplot(211)
plt.plot(t, f_t,"k.")
plt.xlabel("Date")
plt.ylabel("Flux Count")
plt.subplot(212)
plt.plot(t_new, f_new,"b." )
plt.xlabel("Date")
plt.ylabel("Flux Count")
plt.show()

plt.figure()
plt.plot()


'''

#t_thruster, f_thruster, X_thruster, Y_thruster = centroidfit.find_thruster_events(t,f_t,Xc,Yc,starname=starname,outputpath=outputfolder)
t_thruster, f_thruster, X_thruster, Y_thruster = [t,f_t,Xc,Yc]



chunksize = 300


[t_spitzer,f_spitzer] = centroidfit.spitzer_fit(t_thruster, f_thruster, X_thruster, Y_thruster,starname=starname,outputpath=outputpath,chunksize=chunksize)
[t_sff,f_t_sff] = centroidfit.sff_fit(t_thruster, f_thruster, X_thruster, Y_thruster,starname=starname,outputpath=outputpath,chunksize=chunksize)

t_clean1, f_clean1 = centroidfit.clean_data(t_spitzer,f_spitzer)
t_clean2, f_clean2 = centroidfit.clean_data(t_sff,f_t_sff)

STD_Original = np.std(f_t)/np.median(f_t)*1e6
STD_Spitzer = np.std(f_spitzer)*1e6
STD_SFF = np.std(f_t_sff)*1e6


plt.close('all')
plt.figure(figsize=(20,10))
plt.subplot(3,1,1)
#plt.plot(t,f_t,"ro",MarkerSize=3, label="Thruster Flagged data")
plt.plot(t_thruster,f_thruster,"ko",MarkerSize=3, label = "Accepted data")
plt.xlabel("Time(days)")
plt.ylabel("Flux Count")
plt.legend()
plt.title(str(round(STD_Original,1)))
plt.subplot(3,1,2)
#plt.plot(t_spitzer,f_spitzer,"ko",MarkerSize=3,label="Cleaning Flagged data")
plt.plot(t_clean1, f_clean1,"bo-",MarkerSize=3,  label="Spitzer detrended data")
plt.xlabel("Time(days)")
plt.ylabel("Flux Count")
plt.title(str(round(STD_Spitzer,1))+":"+str(chunksize))
plt.legend()
plt.subplot(3,1,3)
#plt.plot(t_sff,f_t_sff,"ko",MarkerSize=3,label="Cleaning Flagged data")
plt.plot(t_clean2, f_clean2,"go-",MarkerSize=3, label= "SFF detrended data")
plt.xlabel("Time(days)")
plt.ylabel("Flux Count")
plt.title(str(round(STD_SFF,1))+":"+str(chunksize))
plt.legend()
plt.tight_layout()
plt.show()
'''


'''
#remove the long term variation in the data with fourth order polynomial
NIter = 10

t_detrend,f_detrend = np.copy([t,f_t])


for i in range(NIter):
    PolyDeg = 4
    params = np.polyfit(t_detrend,f_detrend,PolyDeg)
    y_pred = np.polyval(params,t_detrend)
    Std = np.std(f_detrend)
    SigmaCutOff = 3.5
    IDS = np.abs(f_detrend-y_pred)<SigmaCutOff*Std
    t_detrend,f_detrend = [t_detrend[IDS],f_detrend[IDS]]

X_detrend, Y_detrend = np.copy([Xc[IDS],Yc[IDS]])


def moving_average(series, sigma=3):
    b = gaussian(39, sigma)
    average = filters.convolve1d(series, b/b.sum())
    var = filters.convolve1d(np.power(series-average,2), b/b.sum())
    return average, var

_, var = moving_average(f_detrend)
factor = 1.25
spl =  UnivariateSpline(t_detrend, f_detrend, w=factor/np.sqrt(var))

f_pred = spl(t_detrend)

plt.figure(figsize=(25,4))
plt.plot(t, f_t,"k.",label="Raw Data")
plt.plot(t_detrend,np.polyval(params,t_detrend),"r-",label = "Iterative Sigma Fitting")
plt.plot(t_detrend,f_pred, label="Spline Fitting")
plt.legend()
plt.show()
'''


#t_detrend,f_detrend,X_detrend,Y_detrend = centroidfit.find_thruster_events(t,f_t,Xc,Yc)
#sigma clipping to avoid outliers effect






#period = FindPeriod(t_sff, f_t_sff)




#fitting the transit
#t, f = MaskTransit(t_sff,f_t_sff, period*2)
#IndividualTransitFit(outputfolder,t_sff,f_t_sff)




'''
NumPoints = 100
for i in range(20):
    plt.close()
    plt.figure(figsize=(16,7))
    plt.subplot(1,2,1)
    plt.plot(Xc[i*NumPoints:(i+1)*NumPoints], Yc[i*NumPoints:(i+1)*NumPoints],"ko--")

    plt.subplot(1,2,2)
    plt.plot(t[i*NumPoints:(i+1)*NumPoints], f_t[i*NumPoints:(i+1)*NumPoints],"ko--")
    plt.suptitle(str(i*NumPoints)+":"+str((i+1)*NumPoints))
    plt.show()
'''

'''
plt.close()
plt.figure(figsize=(25,4))
plt.plot(t,f_t,"ko",MarkerSize=7,label="Raw Data")
#plt.plot(t2,f_t2,"rd",MarkerSize=3,label="Spitzer Detrended")
plt.plot(t_sff,f_t_sff,"g*",MarkerSize=2,label="SFF Detrended")
#plt.plot(t_detrend,f_detrend/f_pred-1,"rd-", MarkerSize=2,label="Spline Fitting")
plt.legend()
#plt.grid()
plt.show()
'''


'''
#
#X = np.array([[i,j,k,l] for i,j,k,l in zip(f_t,t,Xc,Yc)])

#f_t = si


#mX = X - np.mean(X,axis=0)
#cvr = np.cov(mX.T)

#pca = PCA(n_components=3)
#pca.fit(X)
#YNew = pca.transform(X)
#print YNew
#plt.plot(x,YNew, "ko")
#show()

# = y/(pca.components_*pca.explained_variance_ratio_)

#s = '{:10.2f} '*N
#for i in range(N):
#    print s.format( *cvr[i] )

# calculate the eigenvalues and eigenvectors
#eigenvalue, eigenvector = np.linalg.eig( cvr )

# sort the eigenvectors according to the eigenvalues
#srt = np.argsort( eigenvalue )[::-1]
#eigenvector = np.matrix( eigenvector[:,srt] )
#eigenvalue = eigenvalue[srt]

#print s.format( *eigenvalue )


#semilogy( eigenvalue.real, '-o' )
#title("Log-plot of Eigenvalues")
#xlabel( "Eigenvalue Index" )
#ylabel( "Eigenvalue Magnitude" )
#grid()
#show()
'''

'''
folded,f_t_folded,period,freqlist,powers = periodfinder.get_period(t_sff,f_t_sff,outputpath=outputpath,starname=starname,get_mandelagolmodel=False)
periodfinder.make_combo_figure(filepath,t_sff,f_t_sff,period,freqlist,powers,starname=starname,outputpath=outputpath)
'''
