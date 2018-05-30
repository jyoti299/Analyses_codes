#!/bin/bash

model=Bstar
coupling=f1p0
Opti=CSVDisc_Opti   ##DEta_Opti  ##DPhi_Opti  ##PhId_Opti  ##JetId_Opti  ##CSVDisc_Opti 

#for cut in  DEta_1.0  DEta_1.2  DEta_1.5  DEta_1.8  DEta_2.0  DEta_2.2  DEta_2.5  NoDEta
#for cut in  DPhi_1.5  DPhi_2.0  DPhi_2.5  NoDPhi
#for cut in PhId_Loose  PhId_Medium  PhId_Tight
#for cut in JetId_Tight  JetId_TightLepVeto
for cut in CSVL  CSVM  CSVT
do

boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016 ##Excited1btagQuarks2016
outdir=Cards_Optimization/${Opti}/${cut}/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

##for Qstar
#python python/GetCombine.py -d ${outdir} -m ${model} --mass range\(500,6600,100\) -b ${boxName} --xsec 0.01 -l 24.487

##for Bstar
python python/GetCombine.py -d ${outdir} -m ${model} --mass range\(500,5000,100\) -b ${boxName} --xsec 0.01 -l 24.487


echo "DONE FOR ${cut}"
done
