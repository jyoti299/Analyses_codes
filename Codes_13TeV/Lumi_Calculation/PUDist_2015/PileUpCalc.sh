#!/bin/bash

cmsenv

#For Golden Json#######################################
#################From Here#############################
#MyJson_Path=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Lumi_Calculation/CMSSW_7_4_14/src/RecoLuminosity/LumiDB/LumiSummaryJSON/forQstar_GoldenJson/Data/crabJSON

#if [ ! -d Data_Run2015D_PileUpDist_GoldenJson ]; then
#echo "Making Directory Data_Run2015D_PileUpDist_GoldenJson"
#mkdir Data_Run2015D_PileUpDist_GoldenJson
#chmod 775 Data_Run2015D_PileUpDist_GoldenJson
#fi

#for i in lumiSummary_2015D-PromptReco_v3  lumiSummary_2015D-PromptReco_v4

#do

#pileupCalc.py -i ${MyJson_Path}/${i}.json --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 ${i}_DataPileupHist.root
#echo "Moving ${i}_DataPileupHist.root to Data_Run2015D_PileUpDist_GoldenJson folder"
#mv ${i}_DataPileupHist.root Data_Run2015D_PileUpDist_GoldenJson/${i}_DataPileupHist.root
#echo "PileUp Distribution Calculation done for ${i}"

#done
#################To Here################################


#For Silver Json########################################
##########From Here#####################################
MyJson_Path=/uscms_data/d3/rocky86/slc6_amd64_gcc491/Lumi_Calculation/CMSSW_7_4_14/src/RecoLuminosity/LumiDB/LumiSummaryJSON/forQstar_SilverJson/Data/crabJSON

if [ ! -d Data_Run2015D_PileUpDist_SilverJson ]; then
echo "Making Directory Data_Run2015D_PileUpDist_SilverJson"
mkdir Data_Run2015D_PileUpDist_SilverJson
chmod 775 Data_Run2015D_PileUpDist_SilverJson
fi

for i in  lumiSummary_2015D-PromptReco_v3_SilverJSON  lumiSummary_2015D-PromptReco_v4_SilverJSON  lumiSummary_2015D-PromptReco_v4_SilverJSON_Extra

do

pileupCalc.py -i ${MyJson_Path}/${i}.json --inputLumiJSON pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 ${i}_DataPileupHist.root
echo "Moving ${i}_DataPileupHist.root to Data_Run2015D_PileUpDist_SilverJson folder"
mv ${i}_DataPileupHist.root Data_Run2015D_PileUpDist_SilverJson/${i}_DataPileupHist.root
echo "PileUp Distribution Calculation done for ${i}"

done
############To Here#####################################

