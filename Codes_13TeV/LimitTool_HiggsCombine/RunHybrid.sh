#!/bin/bash

model=Qstar
coupling=f1p0
boxName=ExcitedQuarks2016  ##ExcitedQuarks2016 ##Excited1btagQuarks2016

config=config/Qstar.config

#outdir=Cards_${model}_${coupling}_${boxName}_Significance/
#outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromHybridNew/
outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromHybridNew_Check/
DataCarddir=Cards_${model}_${coupling}_${boxName}_LimitsFromData_LatestBinning/
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

inFitFile=${DataCarddir}/FitResults_${boxName}.root

for (( MASS=5400; MASS<=5600; MASS+=200 )); do
#for (( MASS=4500; MASS<=4500; MASS+=100 )); do
    echo "WORKING FOR $MASS"

    lumi=35.866
    mass=${MASS}

    combine -M HybridNew --testStat LHC --frequentist -d ${DataCarddir}/ExcitedQuarks_combine_${model}_${mass}_lumi-${lumi}_${boxName}.txt -T 30000 -n ${model}_${mass}_lumi-${lumi}_${boxName} -H Asymptotic --fullBToys --fork 8
#    combine -M HybridNew --testStat LHC --frequentist -d ${DataCarddir}/ExcitedQuarks_combine_${model}_${mass}_lumi-${lumi}_${boxName}.txt -T 1000 --saveHybridResult --setPhysicsModelParameterRanges r=0,1 -n ${model}_${mass}_lumi-${lumi}_${boxName} -H Asymptotic

    mv higgsCombine${model}_${mass}_lumi-${lumi}_${boxName}.HybridNew.mH120.root ${outdir}/higgsCombine${model}_${mass}_lumi-${lumi}_${boxName}.HybridNew.mH120.root

done