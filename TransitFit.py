import model_transits
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
import batman
from auxiliaries import savitzky_golay
from scipy import stats
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import bls

def fold_data(t,y,period):
  folded = t % period
  inds = np.array(folded).argsort()
  t_folded = folded[inds]
  y_folded = y[inds]

  return t_folded,y_folded

def MaskTransit(t_detrended, f_detrended, period, TWidth=0.25):
    '''
    This function is to mask the transit around T0 with width
    '''
    t_new = []
    f_new = []
    T0 = ((t_detrended[np.where(f_detrended==min(f_detrended))] - t_detrended[0])%period + t_detrended[0])[0]

    MaxValue = (int((t_detrended[-1]-t_detrended[0])/period)+1)
    ForbiddenIndex = []

    t_start = t_detrended[0]
    t_end = t_detrended[0]+period

    for i in range(MaxValue):
      Index = np.array(np.array(t_detrended>=t_start)*1+np.array(t_detrended<=t_end)*1)==2
      T_temp = t_detrended[Index]
      F_temp = f_detrended[Index]
      ForbiddenIndex = np.array(abs(T_temp-T0)<TWidth)

      T0+=period
      t_start += period
      t_end += period

      t_new.extend(T_temp[~ForbiddenIndex])
      f_new.extend(F_temp[~ForbiddenIndex])

    return t_new, f_new

    
def SegmentedTest():
    pass


def PlanetModel(LMFitparams,t,Data):
    '''This method fits the transit based on the '''
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = LMFitparams['t0']         #time of inferior conjunction
    params.per = LMFitparams['per']       #orbital period
    params.rp = LMFitparams['rp']         #planet radius (in units of stellar radii)
    params.a = LMFitparams['a']           #semi-major axis (in units of stellar radii)
    params.inc = LMFitparams['inc']       #orbital inclination (in degrees)
    params.ecc = LMFitparams['ecc']       #eccentricity
    params.w = LMFitparams['w']           #longitude of periastron (in degrees)
    #ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
    #ld_coefficients = [[], [0.3], [0.1, 0.3], [0.5, 0.1, 0.1, -0.1]]

    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [LMFitparams['q1'],LMFitparams['q2']]#

    #params.limb_dark = "nonlinear"        #limb darkening model
    #params.u = [LMFitparams['u1'],LMFitparams['u2'],LMFitparams['u3'],LMFitparams['u4']]      #limb darkening coefficients


    m = batman.TransitModel(params, t)    #initializes model
    ModelFlux = m.light_curve(params)     #calculates light curve

    return ModelFlux - Data

def PlanetFit(t,Flux,Per, t0):
    #Initializing the parameters
    LMFitparams = Parameters()
    LMFitparams.add('t0', value= t0, min=t[0], max=t[-1])
    LMFitparams.add('per', value= Per, vary=False)
    LMFitparams.add('rp', value= 0.01, min=0.0001, max=0.35)
    LMFitparams.add('a', value= 12., min= 0)
    LMFitparams.add('inc', value= 90., min=86,max=90.)
    LMFitparams.add('ecc', value=0., min=0, max=0.33)
    LMFitparams.add('w', value= 90)

    #for quadratic fitting
    LMFitparams.add('q1', value= 0.99997, min=0, max=2)#.vary=False)
    LMFitparams.add('q2', value= 0.10319,min=-2,max=2)#,vary=False)

    #For non linear fit
    #LMFitparams.add('u1', value= 0.5, min=-1, max=1)
    #LMFitparams.add('u2', value= 0.1, min=-1, max=1)
    #LMFitparams.add('u3', value= 0.1, min=-1, max=1)
    #LMFitparams.add('u4', value= -0.1, min=-1, max=1)
    fitter = Minimizer(PlanetModel, LMFitparams, fcn_args=(t, Flux))
    result = fitter.minimize(method='nelder') #'leastsq', 'nelder'
    final = Flux + result.residual
    report_fit(result)
    return final

def FindPeriod(t,f_t):
    fmin = 2.0/((t[len(t)-1]-t[0])) # minimum frequency. we can't find anything longer than 90 days obviously
    nf = 2e5 # amount of frequencies to try
    df = 1e-4#0.00001 # frequency step

    qmi = 0.0005 # min relative length of transit (in phase unit)
    qma = 0.1 # max relative length of transit (in phase unit)
    nb = 200 # number of bins in folded LC

    t = np.array(t)
    f_t = np.array(f_t)
    results = bls.eebls(t,f_t,t,f_t,nf,fmin,df,nb,qmi,qma)
    period = results[1]
    return period



def SingleFit(outputfolder,counter,Time,Flux,period,T0, TWidth):
    print "Min Time:: ",Time[0]
    print "Max Time:: ",Time[-1]
    print "T0:: ",T0
    T_Selected = Time[np.where(abs(Time-T0)<=TWidth)]
    F_Selected = Flux[np.where(abs(Time-T0)<=TWidth)]

    ModelB = PlanetFit(T_Selected,F_Selected,period,T0)

    pl.figure()
    pl.plot(T_Selected*24 - 24*T0, F_Selected,'k.',label="Kepler Data")
    pl.plot(T_Selected*24 - 24*T0, ModelB,'r-',lw=3, label="Leastsq Fitting")
    pl.axhline(y=1,color='k',lw=1)
    pl.ylabel('Normalized Flux')
    pl.xlabel('Hours from mid-transit')
    pl.legend(loc='best')
    pl.title(str(counter)+"_Transit")
    pl.savefig(outputfolder+"/"+str(counter)+"_Transit.png")

def IndividualTransitFit(outputfolder,t_detrended, f_detrended):
    t_detrended = np.array(t_detrended)
    f_detrended = np.array(f_detrended)+ 1.0

    period = FindPeriod(t_detrended, f_detrended)

    period = period*2
    T_Fold, F_Fold = fold_data(t_detrended,f_detrended,period)

    T0 = ((t_detrended[np.where(f_detrended==min(f_detrended))] - t_detrended[0])%period + t_detrended[0])[0]

    t_start = t_detrended[0]
    t_end = t_detrended[0]+period

    print t_start
    print t_end

    TWidth = 0.75
    for i in range(int((t_detrended[-1]-t_detrended[0])/period)+1):
        Index = np.array(np.array(t_detrended>t_start)*1+np.array(t_detrended<t_end)*1)==2
        SingleFit(outputfolder,i,t_detrended[Index],f_detrended[Index],period,T0, TWidth)
        T0 += period
        t_start, t_end = t_start+period, t_end+period

    SingleFit(outputfolder,"folded",T_Fold,F_Fold,period, T0%period, TWidth)
