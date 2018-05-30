#!/bin/bash

model=Bstar ##Qstar
coupling=f0p1
boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData/

mMin=1000 #5400 #1000  #5000 #1000
mMax=5500 #5800 #6500 #6500 #5500 #6000

xMin=1e-4 #8e-5 #1e-5
xMax=1e-0 #3e-4 #1e-1 ##5e-4 ##1e-1

python python/Plot1DLimit_SigxBR_only.py -d ${outdir} -m ${model} -b ${boxName} -l 35.866 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --coup ${coupling}
