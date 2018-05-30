#!/bin/bash

model=Bstar
coupling=f1p0
boxName=Excited0btagQuarks2016  ##ExcitedQuarks2016 ##Excited1btagQuarks2016

config=config/Qstar.config

#outdir=Cards_${model}_${coupling}_${boxName}_Significance/
outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyStats/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData_LatestBinning/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyPERSyst/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyLumiSyst/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsAllExceptBSF/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsAllExceptBSF/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyBkgSyst/
#outdir=Cards_PToverMCheck/NoDEta_NoPToverM/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

inFitFile=${outdir}/FitResults_${boxName}.root
#inFitFile=Cards_${model}_${coupling}_${boxName}_LimitsFromData/FitResults_${boxName}.root

##for Qstar
#python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(1000,7500,50\) -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l 35.866 --coup ${coupling} #--signif  #--partial-signal-sys --freezeNuc ,jes,jer,pes,per,bsf,trig #--deco
#python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(1000,6500,100\) -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l 35.866 --coup ${coupling}  --partial-signal-sys --freezeNuc ,jes,jer,pes,per,bsf,trig #--deco
#python python/RunCombine.py -m ${model} -d ${outdir} --mass 2800 -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l 35.866 --coup ${coupling} #--partial-signal-sys --deco #--freezeNuc jes,jer,pes,per,bsf,trig --deco

##for Bstar
python python/RunCombine.py -m ${model} -d ${outdir} --mass 4500 -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l 35.866 --coup ${coupling} --no-sys #--no-signal-sys #--partial-signal-sys --freezeNuc jes,jer,pes,per,bsf,trig #--deco  
#python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(1000,5000,50\) -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l 35.866 --coup ${coupling} #--signif  #--no-signal-sys #--partial-signal-sys --freezeNuc jes,jer,pes,per,bsf,trig #--deco  


###Options required for including different syst:
### For no syst: --no-sys
###For no signal syst: --no-signal-sys
###For some of signal syst: --partial-signal-sys --freezeNuc (comma seperated list of nuisances which are not to include) --deco
## The --deco option is taken, as for adding only some of the signal syst, the best way is to first do a S+B fit and then freeze the nuisances to
## the best fit values, and --deco options does (1) S+B fit (2) decorrelate bkg parameters, we do not need (2) but needs (1) so we use this option.
##But with experience I find out that on using --deco option, the expected limit with syst comes out to be smaller that expected limit wihtout any
##syst for some mass points. so This option should not be used. So not using this --deco option in final results.
