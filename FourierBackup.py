import numpy as np
#import seaborn
import matplotlib.pyplot as plt
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit

def SingleSine(params, x, data):
    """Single sine curve"""
    amp = params['amp']
    shift = params['shift']
    omega = params['omega']
    offset = params['offset']
    model = amp*np.sin(x*omega+shift)+offset
    return model - data

def DoubleSine(params, x, data):
    """Double sine curve"""
    amp1 = params['amp1']
    shift1 = params['shift1']
    omega1 = params['omega1']
    amp2 = params['amp2']
    shift2 = params['shift2']
    omega2 = params['omega2']
    offset = params['offset']
    model = amp1*np.sin(x*omega1+shift1) + amp2*np.sin(x*omega2+shift2)+offset
    return (model - data)**2

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

t = np.linspace(0,25,500)
Flux = 100*np.sin(t)+0.5*np.random.random(len(t))+3000
dt = t[1] - t[0]

PowerSpectrum = np.abs(np.fft.fft(Flux))**2
Freq = np.fft.fftfreq(len(t), dt)
idx = np.where(Freq>0)
PowerSpectrum = PowerSpectrum[idx]
Freq = Freq[idx]*2*np.pi
Index  = np.argsort(PowerSpectrum)
FreqList =  Freq[Index][-5:][::-1]

if 1==2:#len(Freq[np.where(PowerSpectrum>=2e5)])>1:
    print "Case 2"
    params = Parameters()
    params.add('amp1', value= 100, min=0)
    params.add('omega1', value= FreqList[0])
    params.add('shift1', value= 0, min=-np.pi/2., max=np.pi/2)
    params.add('amp2',   value= 2, min=0)
    params.add('omega2', value= FreqList[1])
    params.add('shift2', value= 0.0, min=-np.pi/2., max=np.pi/2)
    params.add('offset', value= 1.0, min=0.)
    minner = Minimizer(DoubleSine, params, fcn_args=(t, Flux))
    result = minner.minimize()

else:
    print "Case 1"
    params = Parameters()
    params.add('amp',   value= 10,  min=0)
    params.add('omega', value= FreqList[0])
    params.add('shift', value= 0.0, min=-np.pi/2., max=np.pi/2)
    params.add('offset', value= 1.0, min=0.)
    minner = Minimizer(SingleSine, params, fcn_args=(t, Flux))
    result = minner.minimize()

final = Flux + result.residual
report_fit(result)

Normalized = Flux/final

plt.figure()
plt.subplot(2,1,1)
plt.plot(t,Flux,'k.',label="Data")
plt.plot(t,final,'r-',label="Model")
plt.legend()
plt.subplot(2,1,2)
plt.plot(t,Flux/final,'k.',label="Data")
plt.show()

print "The selected frequencies:", Freq[np.where(PowerSpectrum>=max(PowerSpectrum)/10.0)]
