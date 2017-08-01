import glob
import numpy as np
import os

Folder = "PhaseCurves"

CSVFile = "Apertures/"+Folder+".csv"
ReadFile = np.loadtxt(CSVFile, delimiter=',', skiprows=1)

NewConfirmationStatus = []

for target in ReadFile:
    KEPID = int(target[0])
    Worked = int(target[1])
    Confirmation = int(target[2])
    if Confirmation == 0 :#S
        print "Confirming ", str(KEPID)
        Location = "Apertures/"+Folder+"/"+str(KEPID)+".png"
        os.system("eog %s" %Location)
        InputVariable = input("Enter 1 for good/0 for bad:")
        if int(InputVariable)==1:
            NewConfirmationStatus.append(1)
        else:
            NewConfirmationStatus.append(0)
    else:
        NewConfirmationStatus.append(1)


KEPID_Array = ReadFile[:,0].astype(np.int32)
Worked = ReadFile[:,1]
NewConfirmationStatus = np.array(NewConfirmationStatus)

#Saving the file
np.savetxt(CSVFile, np.transpose([KEPID_Array, Worked, NewConfirmationStatus]), delimiter=',', header="EPIC_ID,Run,Verified")
