#!/bin/bash

model=Bstar
coupling=f-1p0
boxName=Excited1btagQuarks2016  ##ExcitedQuarks2016  ##Excited1btagQuarks2016

CardDir=Cards_ftest/${model}/
#CardDir=Cards_test/
CardName=ExcitedQuarks_combine_${model}_1000_lumi-35.866_${boxName}.txt

Func=F2
###Func2_par should always be greater than Func1_par
Func1=F2
Func1_par=4
Func1_card=${CardDir}/${Func1}BinnedFit/${CardName}

Func2=F2_5Param
Func2_par=5
Func2_card=${CardDir}/${Func2}BinnedFit/${CardName}

outdir=${CardDir}/${Func}_${Func1_par}${Func2_par}/

if [ ! -d ${outdir} ]; then
echo "--------- Making Directory ${outdir} -----------"
mkdir -pv ${outdir}
chmod 775 ${outdir}
fi

#for q*
#python python/ftest.py -n 61 --p1 ${Func1_par} --p2 ${Func2_par} -t 500 --datacard ${Func1_card} --datacard-alt ${Func2_card} -M FTest -o ${outdir}   --data --lumi 36.813

#for b*
python python/ftest.py -n 56 --p1 ${Func1_par} --p2 ${Func2_par} -t 500 --datacard ${Func1_card} --datacard-alt ${Func2_card} -M FTest -o ${outdir}   --lumi 35.866

#python python/ftest.py -t 500 --datacard ${Func1_card} -M GoodnessOfFit -o ${CardDir}/${Func1}BinnedFit/ --data --lumi 36.813