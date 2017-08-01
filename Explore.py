import glob
from os import system
import re

Folders = ["1000004"]
Folders = ["Output/"+i for i in Folders]

for folder in Folders:
    path_L = folder+"/figs/low_sn/*.png"
    path_H = folder+"/figs/high_sn/*.png"
    file_L = glob.glob(path_L)
    file_H = glob.glob(path_H)
    filenames = file_L+file_H
    print filenames
    for filename in filenames:
        starname= re.search('[0-9]{9}',filename).group(0)
        print starname
        if "spd" in filename:
            starname=starname+"_spd"
        FilePath = folder+"/"+starname
        system("open %s" %(FilePath+"/*.png"))

        input()
