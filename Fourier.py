from __future__ import division
import numpy as np
#import seaborn
import matplotlib.pyplot as plt
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit

def SingleSineModel(params, x, data):
    """Single sine curve"""
    amp = params['amp']
    shift = params['shift']
    omega = params['omega']
    offset = params['offset']
    model = amp*np.sin(x*omega+shift)+offset
    return model - data

def DoubleSineModel(params, x, data):
    """Double sine curve"""
    amp1 = params['amp1']
    shift1 = params['shift1']
    omega1 = params['omega1']
    amp2 = params['amp2']
    shift2 = params['shift2']
    omega2 = params['omega2']
    offset = params['offset']
    model = amp1*np.sin(x*omega1+shift1) + amp2*np.sin(x*omega2+shift2)+offset
    return (model - data)

def ellipsoidalVariation(u1, t1, r1, a, q, i, phi):
  '''
  parameters:
  q = m2/m1
  t1 is the gravity-darkening coefficients for primary
  a is the semimajor axis of the system
  i is the orbital inclination
  r1 is the radius of the primary
  u1 is the linear limb darkening law (values between 0.1 and 0.5 from Parsons et. al (2011))
  '''
  return -3*(15+u1)*(1+t1)*(r1/a)**3*q*np.sin(i)**2/(20*(3-u1))*np.cos(2*phi)


def DoubleSineFit(Time,Flux,Omega=1,ReportFlag = False):
    '''Method for fitting two Sinusoidal'''
    MinFreq = 1/(Time[-1]-Time[0])
    MaxFreq = 1/(Time[1]-Time[0])
    AmpGuess = np.std(Flux)*2
    params = Parameters()
    params.add('amp1', value= AmpGuess, min=0)
    params.add('omega1', value= 10, min = MinFreq, max=MaxFreq)
    params.add('shift1', value= 0, min=-np.pi/2., max=np.pi/2)
    params.add('amp2',   value= AmpGuess/10, min=0)
    params.add('omega2', value= 5, min = MinFreq, max=MaxFreq)
    params.add('shift2', value= 0.0, min=-np.pi/2., max=np.pi/2)
    params.add('offset', value= 1.0, min=0.)
    minner = Minimizer(DoubleSineModel, params, fcn_args=(Time, Flux))
    result = minner.minimize()
    final = Flux + result.residual
    if ReportFlag:
        report_fit(result)
    return final, params

def SingleSineFit(Time,Flux,Omega=1.0, ReportFlag = False):
    '''Method for fitting Single Sinusoidal'''
    MedianFlux = np.median(Flux)
    AmpGuess = np.std(Flux)*2
    print AmpGuess

    print "MedianFlux: ", MedianFlux
    params = Parameters()
    params.add('amp',   value=AmpGuess,  min=0)
    params.add('omega', value= Omega, min=0)
    params.add('shift', value= 0.0, min=-np.pi/2., max=np.pi/2)
    params.add('offset', value= MedianFlux, vary=False)
    minner = Minimizer(SingleSineModel, params, fcn_args=(Time, Flux))
    result = minner.minimize()
    final = Flux + result.residual
    if ReportFlag:
        report_fit(result)
    return final, params

def FourierFit(Time,Flux, Threshold=0.3):
    dt = Time[1]-Time[0]
    PowerSpectrum = np.abs(np.fft.fft(Flux))**2
    Freq = np.fft.fftfreq(len(Time), dt)
    idx = np.where(Freq>0)
    PowerSpectrum = PowerSpectrum[idx]
    Freq = Freq[idx]*2*np.pi
    Index  = np.argsort(PowerSpectrum)
    FreqList =  Freq[Index][-5:][::-1]
    return Freq, PowerSpectrum
