#!/bin/bash

cmsenv

MyJson_Path=/uscms_data/d3/rocky86/slc5_amd64_gcc462/Analyzer/CMSSW_5_3_8_patch3/src/work/Analyzer/test/Data_Run2012ABCD_LumiSummary

if [ ! -d Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v1 ]; then
echo "Making Directory Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v1"
mkdir Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v1
chmod 775 Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v1
fi

if [ ! -d Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v2 ]; then
echo "Making Directory Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v2"
mkdir Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v2
chmod 775 Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v2
fi

if [ ! -d Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v3 ]; then
echo "Making Directory Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v3"
mkdir Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v3
chmod 775 Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v3
fi

if [ ! -d Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v4 ]; then
echo "Making Directory Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v4"
mkdir Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v4
chmod 775 Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v4
fi

for i in Data_Run2012A  Data_Run2012B  Data_Run2012B_missingLumi  Data_Run2012C  Data_Run2012D  Data_Run2012D_missingLumi
do

pixelLumiCalc.py recorded -i ${MyJson_Path}/${i}_lumiSummary.json --hltpath "HLT_Photon150_v1" > ${i}_pixel_RecordedLumi.txt
echo "Moving ${i}_pixel_RecordedLumi.txt to Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v1 folder"
mv ${i}_pixel_RecordedLumi.txt Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v1/${i}_pixel_RecordedLumi.txt
echo "Lumi Calculation done for ${i} for v1"

pixelLumiCalc.py recorded -i ${MyJson_Path}/${i}_lumiSummary.json --hltpath "HLT_Photon150_v2" > ${i}_pixel_RecordedLumi.txt
echo "Moving ${i}_pixel_RecordedLumi.txt to Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v2 folder"
mv ${i}_pixel_RecordedLumi.txt Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v2/${i}_pixel_RecordedLumi.txt
echo "Lumi Calculation done for ${i} for v2"

pixelLumiCalc.py recorded -i ${MyJson_Path}/${i}_lumiSummary.json --hltpath "HLT_Photon150_v3" > ${i}_pixel_RecordedLumi.txt
echo "Moving ${i}_pixel_RecordedLumi.txt to Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v3 folder"
mv ${i}_pixel_RecordedLumi.txt Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v3/${i}_pixel_RecordedLumi.txt
echo "Lumi Calculation done for ${i} for v3"

pixelLumiCalc.py recorded -i ${MyJson_Path}/${i}_lumiSummary.json --hltpath "HLT_Photon150_v4" > ${i}_pixel_RecordedLumi.txt
echo "Moving ${i}_pixel_RecordedLumi.txt to Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v4 folder"
mv ${i}_pixel_RecordedLumi.txt Data_Run2012ABCD_RecordedLumi_ForHLTPhoton150_v4/${i}_pixel_RecordedLumi.txt
echo "Lumi Calculation done for ${i} for v4"

done



