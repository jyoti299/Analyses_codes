#!/bin/bash

model=Qstar ##Qstar
coupling=f-1p0
coup=f1p0
boxName=ExcitedQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

mMin=1000  #5000 #1000
mMax=6500 #6500 #5500 #6000

xMin=1e-4
xMax=1e-1 ##5e-4 ##1e-1

python python/PlotBiasLimit.py -m ${model} -b ${boxName} -l 35.866 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} 
