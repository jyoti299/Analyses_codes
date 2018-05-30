#!/bin/bash

model=Bstar
coupling=f0p1
boxName=Excited1btagQuarks2016_Excited0btagQuarks2016  ##ExcitedQuarks2016 ##Excited1btagQuarks2016
lumi=35.866_35.866

config=config/Qstar.config

outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData_LatestBinning/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsOnlyStats/
#outdir=Cards_${model}_${coupling}_${boxName}_Significance/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

inFitFile1Tag=Cards_${model}_${coupling}_Excited1btagQuarks2016_LimitsFromData_LatestBinning/FitResults_Excited1btagQuarks2016.root
inFitFile0Tag=Cards_${model}_${coupling}_Excited0btagQuarks2016_LimitsFromData_LatestBinning/FitResults_Excited0btagQuarks2016.root

hadd ${outdir}/FitResults_combined.root ${inFitFile1Tag} ${inFitFile0Tag}

inFitFile=${outdir}/FitResults_combined.root

##for Bstar
python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(1000,5050,50\) -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l ${lumi} --coup ${coupling}  #--no-sys  #--signif
#python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(1000,5000,50\) -c ${config} -i ${inFitFile} -b ${boxName} --rMax 2 --xsec 0.01 -l 35.866 --coup ${coupling}  

###Options required for including different syst:
### For no syst: --no-sys
###For no signal syst: --no-signal-sys
###For some of signal syst: --partial-signal-sys --freezeNuc (comma seperated list of nuisances which are not to include) --deco
## The --deco option is taken, as for adding only some of the signal syst, the best way is to first do a S+B fit and then freeze the nuisances to
## the best fit values, and --deco options does (1) S+B fit (2) decorrelate bkg parameters, we do not need (2) but needs (1) so we use this option.
##But with experience I find out that on using --deco option, the expected limit with syst comes out to be smaller that expected limit wihtout any
##syst for some mass points. so This option should not be used. So not using this --deco option in final results.
