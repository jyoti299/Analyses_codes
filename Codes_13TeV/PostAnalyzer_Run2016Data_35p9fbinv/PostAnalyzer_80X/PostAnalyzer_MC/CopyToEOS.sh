#!/bin/bash

#rm condor* Job_* analysis_* *exe PostAnalyzer_MC.*

eosPath=/eos/uscms/store/user/rocky86/13TeV/PA_Results_2015+2016/80X
Dir=Checks/GJQCD_Overlap_Removal/

mkdir -p ${eosPath}/${Dir}
mkdir -p ${eosPath}/${Dir}/MC/GJets_MG
mkdir -p ${eosPath}/${Dir}/MC/DiJet

mv GJets_HT*.root ${eosPath}/${Dir}/MC/GJets_MG/
mv QCD_Pt*.root ${eosPath}/${Dir}/MC/DiJet/
