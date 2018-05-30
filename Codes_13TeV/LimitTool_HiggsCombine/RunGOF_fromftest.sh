#!/bin/bash

model=Qstar
coupling=f-1p0
boxName=ExcitedQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

CardName=ExcitedQuarks_combine_${model}_1000_lumi-36.813_${boxName}.txt

Func=Dijet

outdir=Cards_GOF/${model}/${Func}BinnedFit/

python python/ftest.py -t 500 --datacard ${outdir}/${CardName} -M GoodnessOfFit -o ${outdir}/ --data --lumi 36.813 #--algo KS