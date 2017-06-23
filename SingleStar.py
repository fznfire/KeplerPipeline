import model_transits
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
import batman
from auxiliaries import savitzky_golay
from scipy import stats
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import bls
from Fourier import *


import matplotlib as mpl
mpl.rc('font',**{'family':'sans-serif', 'serif':['Computer Modern Serif'],'sans-serif':['Helvetica'], 'size':15,'weight':500, 'variant':'normal'})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':12, 'color':'k'})
mpl.rc('xtick',**{'major.pad':8,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})


def fold_data(t,y,period):
  folded = t % period
  inds = np.array(folded).argsort()
  t_folded = folded[inds]
  y_folded = y[inds]
  return t_folded,y_folded

def FindPeriod(t,f_t):
    fmin = 2.0/((t[len(t)-1]-t[0])) # minimum frequency. we can't find anything longer than 90 days obviously
    nf = 1e5 # amount of frequencies to try
    df = 1e-4#0.00001 # frequency step

    qmi = 0.0005 # min relative length of transit (in phase unit)
    qma = 0.1 # max relative length of transit (in phase unit)
    nb = 200 # number of bins in folded LC

    t = np.array(t)
    f_t = np.array(f_t)
    results = bls.eebls(t,f_t,t,f_t,nf,fmin,df,nb,qmi,qma)

    powers = results[0]
    period = results[1]
    #freq = np.arange(fmin, df*nf,nf)
    freq = np.linspace(fmin, df*nf,nf)
    print len(freq),len(powers)

    #pl.figure()
    #pl.plot(freq,powers,"k.",MarkerSize=2)
    #pl.title("Period")
    #pl.show()
    return period

def RunStar(EPIC_ID,RUNID):
    #open the file from the Run

    #filepath = "Output/%s/CentroidDetrended.txt" %(str(EPIC_ID))
    #filepath = "Output/%s/%s/RawLightCurve.txt" %(str(RUNID),str(EPIC_ID))
    filepath = "Output/%s/%s/CentroidDetrended.txt" %(str(RUNID),str(EPIC_ID))
    #filepath = "Output/212410755_spd/RawLightCurve.txt"
    #filepath = "Output/%s/%s/Cleaned.txt" %(str(RUNID),str(EPIC_ID))
    Data = np.loadtxt(filepath,delimiter=" ", skiprows=1)
    Time = Data[:,0]
    Flux = Data[:,1]

    #Time = Time[90:]
    #Time = Flux[90:]
    #Time_Folded,Period,freqlist,powers = periodfinder.get_period(Time,Flux,outputpath='',starname='',get_mandelagolmodel=False)
    RealPeriod = FindPeriod(Time, Flux)
    period = 2*RealPeriod
    print "Strongest Period :", period, " days or ", period*24.0, " hours"



    print "Std deviation 1111111    sss:", np.std(Flux)*1e6
    Offset = 1.0-np.median(Flux)
    Flux = Flux+Offset
    T_Fold, F_Fold = fold_data(Time,Flux,period)
    F_Smooth = savitzky_golay(F_Fold,29,1)

    #Freq, PowerSpectrum = FourierFit(Time,Flux)

    #pl.figure()
    #pl.plot(Freq, PowerSpectrum)
    #pl.show()

    Model, param = SingleSineFit(T_Fold,F_Fold,Omega=0.827)
    Leftover = F_Fold - Model



    #Model_New,_ = SingleSineFit(T_Fold,Leftover,Omega=0.424)



    fig = pl.figure(figsize=(15,8))
    sub1 = pl.subplot(2,1,1)
    sub1.patch.set_facecolor('#DBFFFF')
    pl.plot(T_Fold/RealPeriod, F_Fold,"k.")
    #pl.plot(T_Fold/max(T_Fold), Model_New,"r-", lw=3)
    plt.ylim([0.9993,1.0007])
    pl.ylabel('Normalized Flux')
    pl.xlabel('Phase')
    pl.title('Phase Folded Data')

    sub2 = pl.subplot(2,1,2)
    sub2.patch.set_facecolor('#DBFFFF')
    pl.plot(Time,Flux,"k.",MarkerSize=3)
    plt.ylim([0.9993,1.0007])
    pl.ylabel('Normalized Flux')
    pl.xlabel('Time(day)')
    pl.title("Centroid Detrended K2 Data")
    fig.tight_layout()
    fig.patch.set_facecolor('#085000')
    pl.savefig(EPIC_ID+".png")

    dt = Time[1]- Time[0]
    DataPoints = int(period/dt)
    print DataPoints
    StartPoint = 2600
    EndPoint = DataPoints



    #single plot
    '''
    pl.figure(figsize=(20,10))

    for i in range(12):
        pl.subplot(2,6,i+1)
        pl.plot(Time[StartPoint:EndPoint], Flux[StartPoint:EndPoint],"ko",MarkerSize=6)
        StartPoint += DataPoints
        EndPoint += DataPoints
        if i==0 or i==6:
            pl.ylabel("Time(hours)")
        if i>5:
            pl.xlabel("Flux Count")
    pl.show()
    '''
    #subplots

    '''
    ColorName = ["red","blue","green","grey","orange","gold"]
    LabelType = [".","^", "+","1","2","3","4"]
    pl.figure(figsize=(20,6))
    #StartPoint
    TransitNumber = 1
    StartPoint = DataPoints*(TransitNumber-1)+2600
    EndPoint = DataPoints*TransitNumber
    for i in range(5):
        pl.plot(Time[StartPoint+80:EndPoint-30]%period*24, Flux[StartPoint+80:EndPoint-30], Marker=LabelType[i], MarkerSize=10,color=ColorName[i], label=str(TransitNumber+i))
        StartPoint += DataPoints
        EndPoint += DataPoints
    pl.legend()
    pl.xlabel("Time(hours)")
    pl.ylabel("Flux Count")
    pl.show()
    '''

RunStar('212410755', '')
#RunStar('210659779_spd','1000002')
