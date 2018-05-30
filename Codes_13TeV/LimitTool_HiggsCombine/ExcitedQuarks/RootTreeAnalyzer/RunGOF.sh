#!/bin/bash

model=Qstar

fitfun=dijet   ##dijet  ## f1  f2  f3  f4  f5

#config=config/Qstar_${fitfun}.config  ### Qstar_F1.config  Qstar_F2.config  Qstar_F3.config  Qstar_F4.config  Qstar_F5.config
config=config/Qstar.config  ### Qstar_F1.config  Qstar_F2.config  Qstar_F3.config  Qstar_F4.config  Qstar_F5.config

boxname=ExcitedQuarks2016

#outdir=Cards_GOF/${model}/${fitfun}BinnedFit/
outdir=Cards_GOFCheck/TotalMC/
#xoutdir=Cards_GOFCheck/Data/
#outdir=Cards_test/${fitfun}BinnedFit/PseudoData/

python python/RunToys.py -b ${boxname} --freq -c ${config} --lumi 35866 --fit-region Full  -d ${outdir}  -i ${outdir}/FitResults_${boxname}.root -t 1000 -s 0

python python/PlotGOF.py -b ${boxname} -c ${config} -d ${outdir} -t ${outdir}/toys_Freq_s0_${boxname}.root -l 35866 --data