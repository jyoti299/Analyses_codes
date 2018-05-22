#!/bin/bash

cmsenv

MyJson_Path=/uscms_data/d3/rocky86/slc5_amd64_gcc462/Analyzer/CMSSW_5_3_8_patch3/src/work/Analyzer/test/Data_Run2012ABCD_LumiSummary

if [ ! -d Data_Run2012ABCD_RecordedLumi ]; then
echo "Making Directory Data_Run2012ABCD_RecordedLumi"
mkdir Data_Run2012ABCD_RecordedLumi
chmod 775 Data_Run2012ABCD_RecordedLumi
fi

for i in Data_Run2012A  Data_Run2012B  Data_Run2012B_missingLumi  Data_Run2012C  Data_Run2012D  Data_Run2012D_missingLumi
do

pixelLumiCalc.py -i ${MyJson_Path}/${i}_lumiSummary.json overview > ${i}_pixel_RecordedLumi.txt
echo "Moving ${i}_pixel_RecordedLumi.txt to Data_Run2012ABCD_RecordedLumi folder"
mv ${i}_pixel_RecordedLumi.txt Data_Run2012ABCD_RecordedLumi/${i}_pixel_RecordedLumi.txt
echo "Lumi Calculation done for ${i}"

done   
