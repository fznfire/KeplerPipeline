from __future__ import division
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import UnivariateSpline
from scipy.signal import gaussian
from scipy.ndimage import filters

#Everest library input
from everest import Everest
import everest



def EverestReduction(EPIC_ID):
    star = Everest(EPIC_ID, clobber=False) #clobber -- download the file from the internet
    Time = star.time
    Flux = star.flux

    star.cbv_num = 2
    star.compute()

    Flux1 = star.flux
    #star = everest.rPLD(EPIC_ID, clobber = False, mission = 'k2',
                      #giter = 1, gmaxf = 3, lambda_arr = [1e0, 1e5, 1e10], oiter = 3,
                      #pld_order = 2, get_hires = False, get_nearby = False)



    #star = everest.Detrender(EPIC_ID)
    #nPLD is the best for the low brightness star
    #print "Inside Detrender function"
    #star = everest.nPLD(EPIC_ID, clobber = False, mission = 'k2',
    #                  giter = 1, gmaxf = 3, lambda_arr = [1e0, 1e5, 1e10], oiter = 3,
    #                  pld_order = 2, get_hires = False, get_nearby = False)



    pl.figure()
    pl.subplot(211)
    pl.plot(Time,Flux,"ko")
    pl.xlabel("Time")
    pl.ylabel("Flux")
    pl.title("Raw Light Curve")

    pl.subplot(212)
    pl.plot(Time,Flux1,"go")
    pl.xlabel("Time")
    pl.ylabel("Flux")
    pl.title("GP reduced light curve")

    pl.savefig(str(EPIC_ID)+".png")
    pl.close('all')
    raise Exception('Incomplete Buildup yet')
    return 0


def Detrend(t,f_t,Xc,Yc):

    print "Local Detrender"

    #Normalize the flux
    f_t = f_t/np.median(f_t)
    std = np.std(f_t)

    #find the regions of discontinuity and divide the chunk
    Distance = (Xc**2+Yc**2)
    Distance = Distance/np.median(Distance)
    StopPoints = [0]

    #Try polyfitting the data
    ChunkSize = 20
    StepSize = 1/ChunkSize
    ErrorTol = 2e-2
    DegPoly = 2
    x = np.arange(0,1,1/(ChunkSize))
    y = Distance[StopPoints[-1]:StopPoints[-1]+ChunkSize]

    params = np.polyfit(x, y,DegPoly)

    i = 0
    while((Distance[StopPoints[-1]+ChunkSize+1]-np.polyval(params,1+i*StepSize))<ErrorTol):
        i=i+1

    print i

    pl.figure()
    pl.plot(x,np.polyval(params,x),"g-")
    pl.plot(x,y, "ko")
    pl.show()



    #Remove the first clunk of data
    #find the clunks of data


    #first fit fourth order polynomial

    #remove any sinusoidal term
    #t,f_t,Xc,Yc = Detrend()

    pl.figure()
    pl.plot(Distance,"ro", MarkerSize=2)
    #pl.plot(Difference,"g+",MarkerSize=2)
    #pl.plot(f_t, "ko", MarkerSize=2)
    #pl.plot(Xc**2+Yc**2)
    pl.show()

def SplineDetrend(t,f_t,Xc,Yc):
    #Performing detrending using spline
    t = t/100.


    def moving_average(self, series, sigma=3):
        b = gaussian(39, sigma)
        average = filters.convolve1d(series, b/b.sum())
        var = filters.convolve1d(np.power(series-average,2), b/b.sum())
        return average, var*10

    _, var = moving_average(t,f_t)
    sp = UnivariateSpline(t, f_t, w=1/np.sqrt(var))

    #sp =

    #residual = spl(t).get_residual()
    pl.figure()
    #pl.plot(t[300:500],f_t[300:500],"ko",MarkerSize=2)
    #pl.plot(t[300:500],sp(t)[300:500],'r-')
    #pl.plot(residual)
    pl.plot(t,f_t/sp(t))
    pl.show()

    return [t,f_t]
