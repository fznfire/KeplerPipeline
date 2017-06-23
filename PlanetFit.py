from __future__ import division
import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import  model_transits
from lmfit import minimize, Parameters
import bls


params = Parameters()
params.add('Tc', value = 0,vary=False) #time of conjunction for each individual transit
params.add('b', value = 0.1,vary=False) #the impact parameter (b = a cos i/Rstar)
params.add('Rs', value = 0.01,vary=False) #the stellar radius in units of orbital distance (Rstar/a),
params.add('Rp', value = 0.2,vary=False) #planet-to-star radius ratio (Rp/Rstar),
params.add('F', value = 1.0,vary=False) #stellar flux (F0),

#print params['Tc']



P = 10 #Period in days
ModelParams = [0.,0.1,0.008,0.2,1.0, 0.2,0.2]
#ModelParams = [params['Tc'],params['b'],params['Rs'],params['Rp'],params['F']]
t = np.linspace(9.8,10.2,4000)

print "Calculating the model"
#Model1 = (model_transits.modeltransit(ModelParams, model_transits.occultuniform, P, t))#-1.0
Model2 = (model_transits.modeltransit(ModelParams, model_transits.occultquad, P, t))-1.0

'''
print "Model Calculation finished"
fmin = 1 # minimum frequency. we can't find anything longer than 90 days obviously
nf = 4000 # amount of frequencies to try
df = 1e-3# frequency step

qmi = 0.0005 # min relative length of transit (in phase unit)
qma = 0.1 # max relative length of transit (in phase unit)
nb = 20 # number of bins in folded LC

u = np.zeros(len(t))
v = np.zeros(len(t))

results = bls.eebls(t,Model1,u,v,nf,fmin,df,nb,qmi,qma)
power, best_period, best_power, depth, q, in1, in2  = results

print best_power
print depth
'''
#print len(results[0])
#print results[1]
#PowerSpectrum = results[0]
#period = results[1]


pl.figure()
pl.plot(t-10,Model2,"ko-", MarkerSize='2')
#pl.plot(t,PowerSpectrum,"ko",MarkerSize='2')
pl.ylim([min(Model2)-0.05, max(Model2)+0.025])
pl.xla
pl.show()
