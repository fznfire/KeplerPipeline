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


#inputpath = "/home/pniraula/Downloads/Test"+"/*220403421*.fits" #Target with two apparent transits


inputpath = "/home/pniraula/Downloads/Test"+"/*211841249*lpd*.fits" #sinusoidal variation


filepaths = glob.glob(inputpath)
filepath = filepaths[0]
print filepath

starname = str(re.search('[0-9]{9}',filepath).group(0))
outputpath = "DeleteMe/"
outputfolder = outputpath+"/"+starname

t,f_t,Xc,Yc = StandardAperture(filepath,outputpath=outputpath,plot=False)



#t_thruster, f_thruster, X_thruster, Y_thruster = centroidfit.find_thruster_events(t,f_t,Xc,Yc,starname=starname,outputpath=outputfolder)
t_thruster, f_thruster, X_thruster, Y_thruster = [t,f_t,Xc,Yc]



Chunksizes = np.arange(20,1020,20)
print Chunksizes
STD_1 = []
STD_2 = []
for chunksize in Chunksizes:
    print chunksize
    [t_spitzer,f_spitzer] = centroidfit.spitzer_fit(t_thruster, f_thruster, X_thruster, Y_thruster,starname=starname,outputpath=outputpath,chunksize=chunksize)
    [t_sff,f_t_sff] = centroidfit.sff_fit(t_thruster, f_thruster, X_thruster, Y_thruster,starname=starname,outputpath=outputpath,chunksize=chunksize)

    STD_1.append(np.std(f_spitzer)*1e6)
    STD_2.append(np.std(f_t_sff)*1e6)

plt.figure()
plt.plot(Chunksizes,STD_1,"bo-",label="Spitzer")
plt.plot(Chunksizes,STD_2, "go-", label="SFF" )
plt.ylabel("PPM")
plt.xlabel("Chunksize")
plt.title("Study on 211841249")
plt.show()
