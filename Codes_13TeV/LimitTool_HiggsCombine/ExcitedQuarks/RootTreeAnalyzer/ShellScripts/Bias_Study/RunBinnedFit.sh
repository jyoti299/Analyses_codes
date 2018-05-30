#!/bin/bash

model=Bstar
coupling=f-1p0
#tag=0bTag
#tagging=0tag
boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

Func=F6_5Param

config=config/Qstar_${Func}.config
#config=config/Qstar_3Param.config
#config=config/Qstar.config

#outdir=Cards_BiasStudy_Re-miniAOD/${Func}BinnedFit/
#outdir=Cards_BiasStudy/${Func}BinnedFit/
#outdir=Cards_ftest_withPD/${Func}BinnedFit/
outdir=Cards_ftest/${model}/${Func}BinnedFit/
#outdir=Cards_GOF/${model}/${Func}BinnedFit/
#outdir=Cards_${model}_${coupling}_${boxName}_ExpLimitFromData/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

#Sigfile=inputs/Bias_study/ResonanceShapes_Qstar_f-1p0_13TeV_PhLID-JetTID-Pt190-nodeta-nodphi-CSVM_mass695_80X_Spring16.root
##Change this for b*
#Sigfile=inputs/Bias_study/ResonanceShapes_Bstar_1bTag_f-1p0_13TeV_PhLID-JetTID-Pt190-nodeta-nodphi-CSVM_mass695_80X_Spring16.root
Sigfile=inputs/Bias_study/ResonanceShapes_Bstar_1bTag_f1p0_13TeV_PhMID-JetTID-Pt200_170-DEta1p5-noDPhi-CSVL_mass700_80X_Summer16.root

#Datafile=inputs/Bias_study/TotalMC_Run2016_ReReco-BCDEFG_PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_Nodeta_Nodphi_CSVM_Mass695_QstarInvtMass.root
##Change this for b*
#Datafile=inputs/Bias_study/PseudoData_FromMC_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_1tagBstarInvtMass.root
Datafile=inputs/Bias_study/PseudoData_FromMC_Qstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_1tagBstarInvtMass.root


#Datafile=inputs/Bias_study/Data_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_1tagBstarInvtMass.root
#Check for GOF
#Datafile=inputs/ForGOF/TotalMC_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root
#Datafile=inputs/ForGOF/PseudoData_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root
#Datafile=inputs/ForGOF/PseudoData_FromMC_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root
#Datafile=inputs/ForGOF/Data_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root

##Change this for b*
#Datafile=inputs/TotalMC_${model}_${tag}_${coupling}_Spring16_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_${tagging}${model}InvtMass.root

#python python/BinnedFit.py -c ${config} -l 36813 --mass 1000 -m ${model} --xsec 10 -s ${Sigfile} ${Datafile} -b ${boxName} -d ${outdir} --fit-spectrum
python python/BinnedFit.py -c ${config} -l 35866 --mass 1000 -m ${model} --xsec 10 -s ${Sigfile} ${Datafile} -b ${boxName} -d ${outdir} --fit-spectrum

