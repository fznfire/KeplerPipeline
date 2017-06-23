import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',**{'family':'sans-serif', 'serif':['Computer Modern Serif'],'sans-serif':['Helvetica'], 'size':15,'weight':500, 'variant':'normal'})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':12, 'color':'k'})
mpl.rc('xtick',**{'major.pad':8,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})

Runs = [1000001,1000002,1000003,1000004]

SuccessMag = []
FailMag = []
for RUN in Runs:

    #Construct filepath
    filepath = "Output/%s/RunSummary.csv" %(str(RUN))


    with open(filepath) as f:
        for line in f:
            #parsing the text
            try:
                Data = line.split(",")
                try:
                    Mag = float(Data[4])
                    if Data[8] == "Success":
                        SuccessMag.append(Mag)
                    elif Data[8] == "Failed":
                        FailMag.append(Mag)
                    else:
                        raise Exception('Something Funky is going on.')
                except:
                    pass
            except:
                pass


plt.figure(figsize=(15,10))
plt.hist(SuccessMag,label="Successfully Reduced: %s" %(str(len(SuccessMag))))
plt.hist(FailMag, label="Reduction Failed: %s" %(str(len(FailMag))))
plt.xlabel("Kepler Magnitude")
plt.ylabel("Number of targets")
plt.legend()
plt.savefig('Histogram.png')
