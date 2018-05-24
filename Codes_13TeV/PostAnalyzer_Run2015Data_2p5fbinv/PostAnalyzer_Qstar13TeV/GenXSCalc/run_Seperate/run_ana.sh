#!/bin/bash

pwd=${PWD}
ana_Dir=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/GenXSCalc/CMSSW_7_5_0/src
DASDir=${pwd}/DASOut_Dir

OutDir=XSOutput_Dir
if [ ! -d ${OutDir} ]; then
echo "Making Directory ${OutDir}"
mkdir ${OutDir}
chmod 775 ${OutDir}
fi

#Change only this line for different datasets (take it from corresponding run_das.sh)
for i in GJets_HT-40To100  GJets_HT-100To200  GJets_HT-200To400  GJets_HT-400To600  GJets_HT-600ToInf
#for i in GJets_HT-40To100
do

#InputFile=${DASDir}/DAS_${i}.txt
InputFile=DAS_${i}.txt

FileLen=$(wc -l ${InputFile} | cut -c1-2)

q=0
#Initializing an empty variable
inFile=
while [ $q != ${FileLen} ]; do

((q++))

file=$(tail -n +$q ${InputFile} | head -1)

if [ $q != ${FileLen} ]; then
  Addfile=root://xrootd.unl.edu/${file},
else
  Addfile=root://xrootd.unl.edu/${file}
fi

inFile=${inFile}${Addfile}

done

echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "Dataset= ${Dataset}" >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "inputFiles= ${inFile}" >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt
echo "++++++++++++++LOG STARTS HERE FOR ${i}" >> XSOutput_${i}.txt
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
cmsRun ${ana_Dir}/ana.py inputFiles="file:${inFile}" maxEvents=-1 &>> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "                                                     " >> XSOutput_${i}.txt
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt
echo "++++++++++++++LOG ENDS HERE FOR ${i}" >> XSOutput_${i}.txt
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${i}.txt


done
