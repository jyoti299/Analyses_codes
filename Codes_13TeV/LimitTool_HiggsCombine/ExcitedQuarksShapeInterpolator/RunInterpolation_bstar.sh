#!/bin/bash

model=Bstar
coupling=f0p1
btag=
Syst=

for sys in NOMINAL JESUP JESDOWN PESUP PESDOWN JER PER BSFUP BSFDOWN
do

for tag in 0bTag #0bTag
do

echo "### Running for ${sys} and ${tag} ##"

if [ "${sys}" = "NOMINAL" ]; then
    Syst=
else
    Syst=_${sys}
fi

NameShape=${model}_${tag}_${coupling}_13TeV_PhMID-JetTID-Pt200_170-DEta1p5-noDPhi-CSVL_mass700_80X_Summer16${Syst}

Xdist=MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/${model}_${tag}_XDists_${coupling}${Syst}.root

shapeDir=inputs/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/

if [ ! -d ${shapeDir} ]; then
echo "--------- Making Directory ${shapeDir} -----------"
mkdir -pv ${shapeDir}
chmod 775 ${shapeDir}
fi

InterpolateDir=Interpolated_ResShapes/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/

if [ ! -d ${InterpolateDir} ]; then
echo "--------- Making Directory ${InterpolateDir} -----------"
mkdir -pv ${InterpolateDir}
chmod 775 ${InterpolateDir}
fi

./extractShapes.py -i ${Xdist} > ${shapeDir}/input_shapes_${NameShape}.py

./getResonanceShapes.py -i ${shapeDir}/input_shapes_${NameShape}.py -f ${model} --massrange 500 5100 50 -o ${InterpolateDir}/ResonanceShapes_${NameShape}.root --fineBinning

done
done