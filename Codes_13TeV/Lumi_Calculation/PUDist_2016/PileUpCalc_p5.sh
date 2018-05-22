#!/bin/bash

cmsenv

#MyJson_Path=/uscms_data/d3/varun/13TeV/Bstar/CMSSW_8_0_11/src/ggAnalysis/ggNtuplizer/test/forGJanalysis/Data/processedJSONs/21Sept2016_Cert_271036-279931_13TeV_PromptReco_Collisions16_JSON_NoL1T
MyJson_Path=Data_ReReco-BCDEFG_PromptReco-H-Fulll2016/

if [ ! -d Data_ReReco-BCDEFG_PromptReco-H_PileUpDist_p5 ]; then
echo "Making Directory Data_ReReco-BCDEFG_PromptReco-H_PileUpDist_p5"
mkdir Data_ReReco-BCDEFG_PromptReco-H_PileUpDist_p5
chmod 775 Data_ReReco-BCDEFG_PromptReco-H_PileUpDist_p5
fi

for i in  DataB_processedLumis  DataC_processedLumis  DataD_processedLumis  DataE_processedLumis  DataF_processedLumis  DataG_processedLumis  DataH-v2_processedLumis  DataH-v3_processedLumis

do

pileupCalc.py -i ${MyJson_Path}/${i}.json --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 74200 --maxPileupBin 50 --numPileupBins 50 ${i}_DataPileupHist.root
echo "Moving ${i}_DataPileupHist.root to Data_ReReco-BCDEFG_PromptReco-H_PileUpDist_p5 folder"
mv ${i}_DataPileupHist.root Data_ReReco-BCDEFG_PromptReco-H_PileUpDist_p5/${i}_DataPileupHist.root
echo "PileUp Distribution Calculation done for ${i}"

done
############To Here#####################################

