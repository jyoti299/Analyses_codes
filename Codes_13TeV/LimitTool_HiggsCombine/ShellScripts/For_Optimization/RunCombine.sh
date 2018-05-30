#!/bin/bash

model=Bstar
coupling=f1p0
Opti=CSVDisc_Opti  ##DEta_Opti  ##DPhi_Opti  ##PhId_Opti  ##JetId_Opti  ##CSVDisc_Opti 

#for cut in  DEta_1.0  DEta_1.2  DEta_1.5  DEta_1.8  DEta_2.0  DEta_2.2  DEta_2.5  NoDEta
#for cut in  DPhi_1.5  DPhi_2.0  DPhi_2.5  NoDPhi
#for cut in PhId_Loose  PhId_Medium  PhId_Tight
#for cut in JetId_Tight  JetId_TightLepVeto
for cut in CSVL  CSVM  CSVT
do

boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016 ##Excited1btagQuarks2016
FitDir=inputs/Optimization_80X/${Opti}/${cut}/
outdir=Cards_Optimization/${Opti}/${cut}/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

Datafile=$(find inputs/Optimization_80X/${Opti}/${cut}/ -name "Data*.root")

##for Qstar
##python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(500,6600,100\) -c config/Qstar.config -i ${FitDir}/FitResults_${boxName}.root -b ${boxName} --rMax 20 --xsec 0.01 -l 24.487 --Opttype ${Opti} --OptCut ${cut} --datafile ${Datafile}

##for Bstar
python python/RunCombine.py -m ${model} -d ${outdir} --mass range\(500,5000,100\) -c config/Qstar.config -i ${FitDir}/FitResults_${boxName}.root -b ${boxName} --rMax 20 --xsec 0.01 -l 24.487 --Opttype ${Opti} --OptCut ${cut} --datafile ${Datafile}


echo "DONE FOR ${cut}"
done
