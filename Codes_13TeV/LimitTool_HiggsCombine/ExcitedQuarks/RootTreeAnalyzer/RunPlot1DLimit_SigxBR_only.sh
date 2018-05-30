#!/bin/bash

model=Bstar ##Qstar
coupling=f1p0
boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData_LatestBinning/

mMin=1.0 #5400 #1000  #5000 #1000
mMax=7.0 #5800 #6500 #6500 #5500 #6000

xMin=5e-7 #5e-6 #3e-6
xMax=5e-1 #3e-1 #2e-1 ##5e-4 ##1e-1

python python/Plot1DLimit_SigxBR_only.py -d ${outdir} -m ${model} -b ${boxName} -l 35.866 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --coup ${coupling}
#python python/Plot1DLimit_SigxBR_only_alltheroyCurves.py -d ${outdir} -m ${model} -b ${boxName} -l 35.866 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --coup ${coupling}
