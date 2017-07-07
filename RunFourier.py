import numpy as np
import matplotlib.pyplot as pl
import glob
import re
import pandas as pd
import os
from bls import eebls
from auxiliaries import *
#get the list of high SNR stars

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
    power, best_period, best_power, depth, q, in1, in2  = eebls(t,f_t,t,f_t,nf,fmin,df,nb,qmi,qma)
    period = best_period
    SNR = best_power/np.median(power)
    return period, SNR

#InputPaths

Folders = ['Output/1000000','Output/1000001','Output/1000002','Temporary']

for folder in Folders:
    #path = folder+"/figs/high_sn/*.png"
    path = folder+"/figs/low_sn/*.png"
    filenames = glob.glob(path)
    for filename in filenames:
        starname= re.search('[0-9]{9}',filename).group(0) #extracting epic ID number
        if "spd" in filename:
            starname=starname+"_spd"
        print starname

        if not("spd"in starname):
            #LightCurvePath = folder+"/"+starname+"/Cleaned.txt"
            #LightCurvePath = folder+"/"+starname+"/CentroidDetrended.txt"
            LightCurvePath = folder+"/"+starname+"/RawLightCurve.txt"
            Data = pd.read_csv(LightCurvePath,delimiter=' ',skiprows=1, header=None, names=['Time', 'Flux'])
            Time = np.array(Data['Time'])
            Flux = np.array(Data['Flux'])
            MidFlag = (max(Time)+min(Time))/2.0
            T1 = Time[Time<MidFlag]
            F1 = Flux[Time<MidFlag]

            T2 = Time[Time>MidFlag]
            F2 = Flux[Time>MidFlag]

            period1, SNR1 = FindPeriod(T1, F1)
            period2, SNR2 = FindPeriod(T2, F2)

            T1_Fold, F1_Fold = fold_data(T1, F1, period1)
            T2_Fold, F2_Fold = fold_data(T2, F2, period2)

            f_t_smooth1 = savitzky_golay(F1,29,1)
            f_t_smooth2 = savitzky_golay(F2,29,1)

            abs_diff_per = np.abs(period1-period2)/(period1+period2)*200
            if period1<1:
                period1 = period1*24.0
                periodText1 = " Hours"
            else:
                periodText1 = " Days"

            if period2<1:
                period2 = period2*24.0
                periodText2 = " Hours"
            else:
                periodText2 = " Days"

            print"p1: ",  period1, " p2: ", period2, " abs_diff: ", abs_diff_per
            if abs_diff_per<2.0 and not(abs(period1-2.45)<0.1):
                print "Saving Diagram"
                pl.figure(figsize=(16,7))
                pl.subplot(121)
                pl.plot(T1_Fold, F1_Fold, "b.")
                pl.plot(T1_Fold, f_t_smooth1, "r-", lw=2)
                #pl.plot(T1_Fold, F1_Fold,'k.')
                pl.xlabel('Time')
                pl.ylabel('Flux')
                pl.title("Period::"+str(round(period1,2))+periodText1+" SNR::"+str(round(SNR1,2)))
                pl.subplot(122)
                pl.plot(T2_Fold, F2_Fold, "g.")
                pl.plot(T2_Fold, f_t_smooth2, "r-", lw=2)
                pl.xlabel('Time')
                pl.ylabel('Flux')
                pl.title("Period::"+str(round(period2,2))+periodText2+" SNR::"+str(round(SNR2,2)))
                pl.suptitle(starname)
                pl.savefig('RawLightCurve/'+starname+'.png')





#Extract the name of the stars
