#!/bin/bash

model=Bstar
coupling=f1p0
boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016 ##Excited1btagQuarks2016

#outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData/
outdir=Cards_${model}_${coupling}_${boxName}_Significance/


#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyBkgSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyJERSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyJESSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyLumiSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyPERSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyPESSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyStats/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyTrigSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsAllExceptBkg/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyJESnJERSyst/
#outdir=Cards_Qstar_f1p0_ExcitedQuarks2016_LimitsOnlyPESnPERSyst/

#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsAllExceptBSF/               
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsAllExceptBkg/
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsFromData/                   
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyBSFSyst/
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyBkgSyst/                
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyJERSyst/
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyJESSyst/                
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyJESnJERSyst/
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyLumiSyst/               
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyPERSyst/
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyPESSyst/                
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyPESnPERSyst/
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyStats/                  
#outdir=Cards_Bstar_f1p0_Excited1btagQuarks2016_LimitsOnlyTrigSyst/




##for Qstar
#python python/GetCombine.py -d ${outdir} -m ${model} --mass range\(1000,6500,50\) -b ${boxName} --xsec 0.01 -l 35.866 --signif
#python python/GetCombine.py -d ${outdir} -m ${model} --mass 6000 -b ${boxName} --xsec 0.01 -l 36.813

##for Bstar
python python/GetCombine.py -d ${outdir} -m ${model} --mass range\(1000,5000,50\) -b ${boxName} --xsec 0.01 -l 35.866 --signif
#python python/GetCombine.py -d ${outdir} -m ${model} --mass 3200 -b ${boxName} --xsec 0.01 -l 35.866

