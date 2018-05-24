#!/bin/bash

##Don't forget to generate proxy first

DatasetFile=${1}

echo "Dataset = ${DatasetFile}"
ana_Dir=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Analyzer_13TeV/GenXSCalc/CMSSW_7_5_0/src

if [ -f XSOutput_${DatasetFile} ]; then
echo "++++++++++++++ Deleting XSOutput_${DatasetFile} ++++++++++++++"
rm XSOutput_${DatasetFile}
fi

FileLen=$(wc -l ${DatasetFile} | cut -c1-2)

q=0
while [ $q != ${FileLen} ]; do

((q++))
Dataset=$(tail -n +$q ${DatasetFile} | head -1)

nFiles=30

if [ -f DAS_out.txt ]; then
echo "++++++++++++++ Deleting DAS_out.txt ++++++++++++++"
rm DAS_out.txt
fi

if [ -f Input.txt ]; then
echo "++++++++++++++ Deleting Input.txt ++++++++++++++"
rm Input.txt
fi

${ana_Dir}/das_client.py --query="file dataset=${Dataset}" --limit=${nFiles} > DAS_out.txt

tail -n +4 DAS_out.txt > Input.txt

InputFile=Input.txt

fileLength=$(wc -l ${InputFile} | cut -c1-2)

if [ ${fileLength} != ${nFiles} ]; then
echo "ERROR: FILELENGTH NOT MATCHING TO THE NFILES"
fi

p=0
#Initializing an empty variable
inFile=
while [ $p != $fileLength ]; do

((p++))

file=$(tail -n +$p ${InputFile} | head -1)

if [ $p != ${fileLength} ]; then
  Addfile=root://xrootd.unl.edu/${file},
else
  Addfile=root://xrootd.unl.edu/${file}
fi

inFile=${inFile}${Addfile}

done


echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "Dataset= ${Dataset}" >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "inputFiles= ${inFile}" >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "++++++++++++++LOG STARTS HERE FOR ${Dataset}" >> XSOutput_${DatasetFile}
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
cmsRun ${ana_Dir}/ana.py inputFiles="file:${inFile}" maxEvents=-1 &>> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "++++++++++++++LOG ENDS HERE FOR ${Dataset}" >> XSOutput_${DatasetFile}
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++" >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}
echo "                                                     " >> XSOutput_${DatasetFile}


done