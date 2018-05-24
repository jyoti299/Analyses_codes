#!/bin/bash

#generating proxy required for das_client
voms-proxy-init --voms cms

pwd=${PWD}
ana_Dir=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/GenXSCalc/CMSSW_7_5_0/src

OutDir=DASOut_Dir
if [ ! -d ${OutDir} ]; then
echo "Making Directory ${OutDir}"
mkdir ${OutDir}
chmod 775 ${OutDir}
fi

#Change only these two lines for different datasets.
MCFile=GJets_Datasets.txt
for i in GJets_HT-40To100  GJets_HT-100To200  GJets_HT-200To400  GJets_HT-400To600  GJets_HT-600ToInf

do

Dataset=$(grep "${i}" ${MCFile})

nFiles=30

${ana_Dir}/das_client.py --query="file dataset=${Dataset}" --limit=${nFiles} > ${OutDir}/temp.txt

tail -n +4 ${OutDir}/temp.txt > ${OutDir}/DAS_${i}.txt

done

if [ -f ${OutDir}/temp.txt ]; then
echo "++++++++++++++ Deleting temp.txt ++++++++++++++"
rm ${OutDir}/temp.txt
fi
