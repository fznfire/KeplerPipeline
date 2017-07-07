import model_transits
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
import batman
from auxiliaries import savitzky_golay
from scipy import stats
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import bls

'''
import matplotlib as mpl
mpl.rc('font',**{'family':'sans-serif', 'serif':['Computer Modern Serif'],'sans-serif':['Helvetica'], 'size':15,'weight':500, 'variant':'normal'})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':12, 'color':'k'})
mpl.rc('xtick',**{'major.pad':8,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})
'''

def fold_data(t,y,period):
  folded = t % period
  inds = np.array(folded).argsort()
  t_folded = folded[inds]
  y_folded = y[inds]
  return t_folded,y_folded
  
def MaskTransit():
    pass

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
    LMFitparams.add('t0', value= t0, min=0)
    LMFitparams.add('per', value= Per, vary=False)
    LMFitparams.add('rp', value= 0.15, min=0)
    LMFitparams.add('a', value= 12., min= 0)
    LMFitparams.add('inc', value= 90., min=86,max=90.)
    LMFitparams.add('ecc', value=0., min=0, max=0.33)
    LMFitparams.add('w', value= 90, vary=False)

    #for quadratic fitting
    LMFitparams.add('q1', value= 0.99997, min=0, max=2)#.vary=False)
    LMFitparams.add('q2', value= 0.10319,min=-2,max=2)#,vary=False)

    #For non linear fit
    #LMFitparams.add('u1', value= 0.5, min=-1, max=1)
    #LMFitparams.add('u2', value= 0.1, min=-1, max=1)
    #LMFitparams.add('u3', value= 0.1, min=-1, max=1)
    #LMFitparams.add('u4', value= -0.1, min=-1, max=1)
    fitter = Minimizer(PlanetModel, LMFitparams, fcn_args=(t, Flux))
    result = fitter.minimize(method='leastsq') #'leastsq', 'nelder'
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

def EvenOddTest(EPIC_ID,Time,Flux,period,T0, TWidth):
    '''This function fits the parameters for even and odd period of the orbit'''
    print "Starting Even Odd Test"
    T_Fold, F_Fold = fold_data(Time,Flux,period)
    T_2Period,F_2Period = fold_data(Time,Flux,2.0*period)
    T_Even = T_2Period[:int(len(T_2Period)/2.0)]
    F_Even = F_2Period[:int(len(T_2Period)/2.0)]
    T_Odd = T_2Period[int(len(T_2Period)/2.0):]
    F_Odd = F_2Period[int(len(T_2Period)/2.0):]
    del T_2Period, F_2Period


    T_FoldBoth = T_Fold[np.where(abs(T_Fold-T0)<=TWidth)]
    F_FoldBoth = F_Fold[np.where(abs(T_Fold-T0)<=TWidth)]

    #Odd Even Planet fitting
    T_EvenFit = T_Even[np.where(abs(T_Even-T0)<=TWidth)]# & (T_Even<T_High))]
    F_EvenFit = F_Even[np.where(abs(T_Even-T0)<=TWidth)]# & (T_Even<T_High))]

    T_OddFit = T_Odd[np.where(abs(T_Odd-(T0+period))<=TWidth)]# & (T_Odd<period+T_High))]
    F_OddFit = F_Odd[np.where(abs(T_Odd-(T0+period))<=TWidth)]# & (T_Odd<period+T_High))]

    print "Fitting the Models"
    ModelB = PlanetFit(T_FoldBoth,F_FoldBoth,period, T0)
    print "Second Model"
    ModelE = PlanetFit(T_EvenFit,F_EvenFit,period, T0)
    print "Third Model"
    ModelO = PlanetFit(T_OddFit,F_OddFit,period, T0)

    print "Plotting"
    pl.figure(figsize=(18,10))
    pl.subplot(3,1,1)
    pl.plot(T_FoldBoth,F_FoldBoth,'k.',label="Kepler Data")
    pl.plot(T_FoldBoth,ModelB,'r-',lw=3, label="Leastsq Fitting")
    #pl.axhline(y=1,color='b',lw=2)
    pl.ylabel('Normalized Flux')
    #pl.xlabel('Time(hours)')
    pl.legend(loc='best')

    pl.subplot(3,1,2)
    pl.plot(T_EvenFit,F_EvenFit,'k.',label="Even Period Data")
    pl.plot(T_EvenFit,ModelE,'g-',lw=3, label="Leastsq Fitting")
    #pl.axhline(y=1,color='b',lw=2)
    pl.ylabel('Normalized Flux')
    #pl.xlabel('Time(hours)')
    pl.legend(loc='best')

    pl.subplot(3,1,3)
    pl.plot((T_OddFit),F_OddFit,'k.',label="Odd Period Data")
    pl.plot((T_OddFit),ModelO,'b-',lw=3, label="Leastsq Fitting")
    #pl.axhline(y=1,color='b',lw=2)
    pl.ylabel('Normalized Flux')
    pl.xlabel('Time(hours)')
    pl.legend(loc='best')
    #pl.suptitle('EPIC:'+EPIC_ID)
    pl.savefig(EPIC_ID+'__Fit.png')

def EvenOddTest_histogram(EPIC_ID,T_Fold, F_Fold, period,T_Even, F_Even):
    BinSize = 30
    bin_means, bin_edges, binnumber = stats.binned_statistic(T_Fold,F_Fold,statistic='mean',bins=BinSize)
    bin_meansE, bin_edgesE, binnumberE = stats.binned_statistic(T_Even,F_Even,statistic='mean',bins=BinSize)
    bin_meansO, bin_edgesO, binnumberO = stats.binned_statistic(T_Odd,F_Odd,statistic='mean',bins=BinSize)

    pl.figure(figsize=(15,12))
    pl.subplot(3,1,1)
    pl.hlines(bin_means-1.0, bin_edges[:-1], bin_edges[1:], colors='g', lw=2)
    pl.vlines(bin_edges[1:-1],bin_means[:-1]-1.0,bin_means[1:]-1.0, colors='g', lw=2, label='Phase Folded Light Curve')
    #pl.plot(T_Fold,Model, label="Leastsq fitting")
    pl.axhline(y=0,color='k',lw=1)
    pl.ylabel('Normalized Flux')
    pl.xlabel('Time(hours)')
    pl.legend(loc='center left')
    if DayFlag:
        pl.xlabel('Time [d]')
    else:
        pl.xlabel('Time [hr]')
    pl.ylabel('Normalized Flux')


    pl.subplot(3,1,2)
    pl.hlines(bin_meansE-1.0, bin_edgesE[:-1], bin_edgesE[1:], colors='r', lw=2)
    pl.vlines(bin_edgesE[1:-1],bin_meansE[:-1]-1.0,bin_meansE[1:]-1.0, colors='r', lw=2, label='Even Phase')
    #pl.plot(T_Fold,Model, label="Leastsq fitting")
    pl.axhline(y=0,color='k',lw=1)
    pl.ylabel('Normalized Flux')
    pl.xlabel('Time(hours)')
    pl.legend(loc='center left')
    if DayFlag:
        pl.xlabel('Time [d]')
    else:
        pl.xlabel('Time [hr]')
    pl.ylabel('Normalized Flux')

    pl.subplot(3,1,3)
    pl.hlines(bin_meansO-1.0, bin_edgesO[:-1], bin_edgesO[1:], colors='r', lw=2)
    pl.vlines(bin_edgesO[1:-1],bin_meansO[:-1]-1.0,bin_meansO[1:]-1.0, colors='r', lw=2, label='Odd Phase')
    #pl.plot(T_Fold,Model, label="Leastsq fitting")
    pl.axhline(y=0,color='k',lw=1)
    pl.ylabel('Normalized Flux')
    pl.xlabel('Time(hours)')
    pl.legend(loc='center left')
    if DayFlag:
        pl.xlabel('Time [d]')
    else:
        pl.xlabel('Time [hr]')
    pl.ylabel('Normalized Flux')
    pl.suptitle("EPIC ID:"+str(EPIC_ID)+" Period: "+str('%0.3f' %(period)))
    pl.savefig(str(EPIC_ID)+'_Even_Odd.png')


def SingleFit(EPIC_ID,Time,Flux,period,T0, TWidth):
    T_Fold, F_Fold = fold_data(Time,Flux,period)
    T_FoldBoth = T_Fold[np.where(abs(T_Fold-T0)<=TWidth)]
    F_FoldBoth = F_Fold[np.where(abs(T_Fold-T0)<=TWidth)]
    ModelB = PlanetFit(T_FoldBoth,F_FoldBoth,period)

    pl.figure()
    pl.plot(T_FoldBoth,F_FoldBoth,'k.',label="Kepler Data")
    pl.plot(T_FoldBoth,ModelB,'r-',lw=3, label="Leastsq Fitting")
    pl.axhline(y=1,color='k',lw=1)
    pl.ylabel('Normalized Flux')
    pl.xlabel('Time(hours)')
    pl.legend(loc='best')
    pl.show()



def RunStar(EPIC_ID,RUNID):
    #open the file from the Run
    #filepath = "Output/%s/%s/RawLightCurve.txt" %(str(RUNID),str(EPIC_ID))
    filepath = "Output/%s/%s/Cleaned.txt" %(str(RUNID),str(EPIC_ID))
    Data = np.loadtxt(filepath,delimiter=" ", skiprows=1)
    Time = Data[:,0]
    Flux = Data[:,1]

    period = FindPeriod(Time, Flux)

    print period
    DayFlag = True

    if period<1.0:
        Time = Time*24.0
        period = period*24.0
        DayFlag = False
        print "Strongest Period :", period, " hours"
    else:
        print "Strongest Period :", period, " days or ", period*24.0, " hours"

    #period = 2.3968



    Offset = 1.0-np.median(Flux)
    Flux = Flux+Offset


    T_Fold, F_Fold = fold_data(Time,Flux,period)
    #T_2Period,F_2Period = fold_data(Time,Flux,2.0*period)

    pl.figure()
    pl.plot(T_Fold, F_Fold,"ko",MarkerSize=2,alpha=0.3)
    pl.show()

    #Determing T0 and TWidth
    T0 = 1.0 #Find the time since conjunction
    TWidth = 0.6

    #input("Do you have T0 and TWidth right?")
    EvenOddTest(EPIC_ID,Time,Flux,period, T0, TWidth)
    #SingleFit(EPIC_ID,Time,Flux,period, T0, TWidth)



#RunStar('201563164','1000001')

#RunStar('229228879','1000001') #Interesting target
#RunStar('220155299_spd','1000001')
RunStar('210659779_spd','1000002')
#RunStar('210659779','1000002')
#RunStar('210659779','1000002')
#RunStar('201273498_spd','1000004')
#RunStar('210659779','1000002')

#RunStar('220403421','1000001')
#RunStar('210412981_spd','1000002')

#RunStar('212410755_spd', '1000006')
#RunStar('212628097_spd', '1000006')
