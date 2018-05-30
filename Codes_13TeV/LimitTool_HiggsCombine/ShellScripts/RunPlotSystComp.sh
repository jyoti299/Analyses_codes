#!/bin/bash

model=Bstar ##Qstar
coupling=f1p0

boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

mMin=1000  #1000  #5000 #1000
mMax=5000 #6500 #6000 #5500 #6000

xMin=1e-5 #1e-5
xMax=1e-1 #1e-1 ##5e-4 ##1e-1

python python/PlotSystComp.py -m ${model} -b ${boxName} -l 35.866 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --coup ${coupling}
