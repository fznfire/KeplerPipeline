import matplotlib
matplotlib.use('Qt4Agg')

import numpy as np
import glob
import matplotlib.pyplot as plt
from ExtractFlux import StandardAperture

import centroidfit
import copy
import re
from scipy.interpolate import UnivariateSpline
from scipy.signal import gaussian
from scipy.ndimage import filters, convolve, measurements


from pixeltoflux import get_lightcurve


inputpath = "/home/pniraula/Downloads/PhaseCurveFitsFiles/"+"/*.fits" #For testing

filepaths = glob.glob(inputpath)
outputpath = "Spline/"
for filepath in filepaths:
    EPIC_ID = re.search('[0-9]{9}',filepath).group(0)
    t,f_t,Xc,Yc = StandardAperture(filepath,outputpath=outputpath,plot=False)

    def moving_average(series, sigma=3):
        b = gaussian(12, sigma)
        average = filters.convolve1d(series, b/b.sum())
        var = filters.convolve1d(np.power(series-average,2), b/b.sum())
        return average, var

    _, var = moving_average(f_t)
    factor = 1.25
    spl =  UnivariateSpline(t, f_t, w=factor/np.sqrt(var))

    f_pred = spl(t)

    plt.figure()
    plt.plot(t,f_t/f_pred,"k.")
    plt.show()
