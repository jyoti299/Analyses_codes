#!/bin/bash

pwd=$PWD
echo "${pwd}"

File_Dir=/uscms_data/d3/rocky86/13TeV/PostAnalyzer_Qstar13TeV/LimitCode/CMSSW_7_0_9/src/Significance/Sig_f0p1_global

COUNT=5000

for ((i=1;i<=COUNT;i++)); do

Inseed=7629
JobID=$[Inseed*i]

for M in  1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000
do 

coupling=f0p1

###Don't Forget to Complile first

./stats ${M} qq ${JobID} > stats_${M}_${coupling}_${i}.log

echo "****************"
echo "****************"
echo "Done for ${i}"
echo "****************"
echo "****************"

rm stats*.root

done



done




