#!/bin/bash

cmsenv

#MyJson_Path=/uscms_data/d3/varun/13TeV/Bstar/CMSSW_8_0_11/src/ggAnalysis/ggNtuplizer/test/forGJanalysis/Data/processedJSONs/21Sept2016_Cert_271036-279931_13TeV_PromptReco_Collisions16_JSON_NoL1T
#MyJson_Path=Data_ReReco-BCDEFG_PromptReco-H-Fulll2016/
MyJson_Path=Data_Re-miniAOD-3Feb_Full2016/

if [ ! -d Data_Re-miniAOD-3Feb_Full2016_PileUpDist ]; then
echo "Making Directory Data_Re-miniAOD-3Feb_Full2016_PileUpDist"
mkdir Data_Re-miniAOD-3Feb_Full2016_PileUpDist
chmod 775 Data_Re-miniAOD-3Feb_Full2016_PileUpDist
fi

#for i in  DataB_processedLumis  DataC_processedLumis  DataD_processedLumis  DataE_processedLumis  DataF_processedLumis  DataG_processedLumis  DataH-v2_processedLumis  DataH-v3_processedLumis
for i in  singlePhoton_2016B_processedLumis  singlePhoton_2016C_processedLumis  singlePhoton_2016D_processedLumis  singlePhoton_2016E_processedLumis  singlePhoton_2016F_processedLumis  singlePhoton_2016G_processedLumis  singlePhoton_2016Hv2_processedLumis  singlePhoton_2016Hv3_processedLumis

do

pileupCalc.py -i ${MyJson_Path}/${i}.json --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 50 --numPileupBins 50 ${i}_DataPileupHist.root
echo "Moving ${i}_DataPileupHist.root to Data_Re-miniAOD-3Feb_Full2016_PileUpDist folder"
mv ${i}_DataPileupHist.root Data_Re-miniAOD-3Feb_Full2016_PileUpDist/${i}_DataPileupHist.root
echo "PileUp Distribution Calculation done for ${i}"

done

