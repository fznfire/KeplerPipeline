import matplotlib
matplotlib.use('Qt4Agg')

import numpy as np
import glob
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, KernelPCA
from ExtractFlux import StandardAperture
from auxiliaries import *
import centroidfit
import matplotlib
import copy
import re
from scipy.interpolate import UnivariateSpline
from scipy.signal import gaussian
from scipy.ndimage import filters
from astropy.io import fits
from scipy.ndimage import convolve, measurements


x = np.linspace(0,10,50)
y = 0.001*x**1.02+np.sin(x)+ 0.2* np.random.random(len(x))+1000

y_real = 0.001*x**1.02+np.sin(x)+1000
y[10] = 990
def moving_average1(series, sigma=3):
    b = gaussian(101, sigma)
    average = filters.convolve1d(series, b/b.sum())
    var = filters.convolve1d(np.power(series-average,2), b/b.sum())
    return average, var

def moving_average2(series, sigma=3):
    b = gaussian(6, sigma)
    average = filters.convolve1d(series, b/b.sum())
    var = filters.convolve1d(np.power(series-average,2), b/b.sum())
    return average, var

_, var1 = moving_average1(y)
_, var2 = moving_average2(y)
factor = 1.00
spl1 =  UnivariateSpline(x, y, w=factor/np.sqrt(var1))
spl2 =  UnivariateSpline(x, y, w=factor/np.sqrt(var2))

plt.figure()
plt.plot(x,y, "ko")
plt.plot(x,y_real,"k-")
plt.plot(x,spl1(x),"rd-",MarkerSize = 2)
plt.plot(x,spl2(x),"g*-",MarkerSize = 2)
plt.show()
