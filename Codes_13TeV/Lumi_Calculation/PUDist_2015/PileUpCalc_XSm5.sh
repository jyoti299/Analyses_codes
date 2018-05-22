#!/bin/bash

cmsenv

#For Silver Json########################################
MyJson_Path=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Lumi_Calculation/CMSSW_7_4_14/src/RecoLuminosity/LumiDB/LumiSummaryJSON/forQstar_SilverJson/Data/crabJSON

if [ ! -d Data_Run2015D_PileUpDist_SilverJson_XSm5 ]; then
echo "Making Directory Data_Run2015D_PileUpDist_SilverJson_XSm5"
mkdir Data_Run2015D_PileUpDist_SilverJson_XSm5
chmod 775 Data_Run2015D_PileUpDist_SilverJson_XSm5
fi

for i in  lumiSummary_2015D-PromptReco_v3_SilverJSON  lumiSummary_2015D-PromptReco_v4_SilverJSON  lumiSummary_2015D-PromptReco_v4_SilverJSON_Extra

do

pileupCalc.py -i ${MyJson_Path}/${i}.json --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 65550 --maxPileupBin 50 --numPileupBins 50 ${i}_DataPileupHist_XSm5.root  ## for XS - 5% of XS. 

echo "Moving ${i}_DataPileupHist_XSm5.root to Data_Run2015D_PileUpDist_SilverJson_XSm5 folder"

mv ${i}_DataPileupHist_XSm5.root Data_Run2015D_PileUpDist_SilverJson_XSm5/${i}_DataPileupHist_XSm5.root
echo "PileUp Distribution Calculation done for ${i}"

done



