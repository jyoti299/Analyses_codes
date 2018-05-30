#!/bin/bash

model=Qstar
coupling=f-1p0
tag=
Syst=

for opti in DEta_1.0  DEta_1.2  DEta_1.5  DEta_1.8  DEta_2.0  DEta_2.2  DEta_2.5  NoDEta
do

for sys in NOMINAL JESUP JESDOWN PESUP PESDOWN JER PER
do

echo "### Running for ${opti} ##"
echo "### Running for ${sys} ##"

if [ "${sys}" = "NOMINAL" ]; then
    Syst=
else
    Syst=_${sys}
fi

NameShape=${model}_${coupling}_13TeV_PhLID-JetTID-Pt190-nodeta-nodphi-CSVM_mass695_80X_Spring16${Syst}

##For normal
#Xdist=MassXDist/Spring16_36813pb_80X/${model}_XDists_${coupling}${Syst}.root
## For optimization
Xdist=MassXDist/Spring16_36813pb_80X/Optimization/DEta_Opti/${opti}/${model}_XDists_${coupling}${Syst}.root

shapeDir=inputs/Optimization/DEta_Opti/${opti}/

if [ ! -d ${shapeDir} ]; then
echo "--------- Making Directory ${shapeDir} -----------"
mkdir -pv ${shapeDir}
chmod 775 ${shapeDir}
fi

InterpolateDir=Interpolated_ResShapes/Optimization/DEta_Opti/${opti}/

if [ ! -d ${InterpolateDir} ]; then
echo "--------- Making Directory ${InterpolateDir} -----------"
mkdir -pv ${InterpolateDir}
chmod 775 ${InterpolateDir}
fi

##For normal
#./extractShapes.py -i ${Xdist} > inputs/input_shapes_${NameShape}.py
## For optmization
./extractShapes.py -i ${Xdist} > ${shapeDir}/input_shapes_${NameShape}.py

## For normal
##./getResonanceShapes.py -i inputs/input_shapes_${NameShape}.py -f Qstar --massrange 500 6900 100 -o Interpolated_ResShapes/ResonanceShapes_${NameShape}.root --fineBinning
## For optimization
./getResonanceShapes.py -i ${shapeDir}/input_shapes_${NameShape}.py -f Qstar --massrange 500 8000 100 -o ${InterpolateDir}/ResonanceShapes_${NameShape}.root --fineBinning

done
done