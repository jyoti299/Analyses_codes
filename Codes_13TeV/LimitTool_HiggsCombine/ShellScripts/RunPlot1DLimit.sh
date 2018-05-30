#!/bin/bash

model=Bstar ##Qstar
coupling=f1p0
boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks201

outdir=Cards_${model}_${coupling}_${boxName}_Significance/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyStats/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyBkgSyst/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsAllExceptBkg/

mMin=1000 #1000  #5000 #1000
mMax=5000 #5800 #6500 #6500 #5500 #6000

xMin=2e-3 #8e-5 #1e-5
xMax=1e-1 #3e-4 #1e-1 ##5e-4 ##1e-1

python python/Plot1DLimit.py -d ${outdir} -m ${model} -b ${boxName} -l 35.866 --massMin ${mMin} --massMax ${mMax} --xsecMin ${xMin} --xsecMax ${xMax} --coup ${coupling} --signif
