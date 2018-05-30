#!/bin/bash

model=Qstar
coupling=f1p0
#Opti=CSVDisc_Opti  ##DEta_Opti  ##DPhi_Opti  ##PhId_Opti  ##JetId_Opti  ##CSVDisc_Opti 

#for cut in  DEta_1.0  DEta_1.2  DEta_1.5  DEta_1.8  DEta_2.0  DEta_2.2  DEta_2.5  NoDEta
#for cut in  DPhi_1.5  DPhi_2.0  DPhi_2.5  NoDPhi
#for cut in PhId_Loose  PhId_Medium  PhId_Tight
#for cut in JetId_Tight  JetId_TightLepVeto
#for cut in CSVL  CSVM  CSVT
#do

boxName=ExcitedQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

##NO change below
outdir=inputs/Optimization_80X/${Opti}/${cut}/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

Sigfile=$(find inputs/Optimization_80X/${Opti}/${cut}/ -name "*${cut}.root")
Datafile=$(find inputs/Optimization_80X/${Opti}/${cut}/ -name "Data*.root")


python python/BinnedFit.py -c config/Qstar.config -l 24487 --mass 1000 -m ${model} -s ${Sigfile} ${Datafile} -b ${boxName} -d ${outdir} --fit-spectrum


#echo "DONE FOR ${cut}"
#done