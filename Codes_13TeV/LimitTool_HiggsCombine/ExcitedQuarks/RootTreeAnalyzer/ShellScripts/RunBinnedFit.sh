#!/bin/bash

model=Bstar

coupling=f0p1
tag=0bTag

boxName=Excited0btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

config=config/Qstar.config

outdir=Cards_${model}_${coupling}_${boxName}_LimitsFromData/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

##for q*
#Sigfile=inputs/ResonanceShapes_Qstar_${coupling}_13TeV_PhMID-JetTID-Pt200_170-DEta1p5-noDPhi-CSVL_mass700_80X_Summer16.root 
##for b*
Sigfile=inputs/ResonanceShapes_Bstar_${tag}_${coupling}_13TeV_PhMID-JetTID-Pt200_170-DEta1p5-noDPhi-CSVL_mass700_80X_Summer16.root

##for q*
#Datafile=inputs/Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Mass700_QstarInvtMass.root
#for b*
#Datafile=inputs/Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Mass700_1tagBstarInvtMass.root
Datafile=inputs/Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Mass700_0tagBstarInvtMass.root

python python/BinnedFit.py -c ${config} -l 35866 --mass 1000 -m ${model} --xsec 10 -s ${Sigfile} ${Datafile} -b ${boxName} -d ${outdir} --fit-spectrum

