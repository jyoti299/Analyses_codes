#!/bin/bash

#rm condor* Job_* analysis_* *exe PostAnalyzer_MC.*

eosPath=/eos/uscms/store/user/rocky86/13TeV/PA_Results_2015+2016/80X
Dir=PhLID_JetTID_Pt190_NoDEta_NoDPhi_NoMassCut

mkdir -p ${eosPath}/${Dir}
mkdir -p ${eosPath}/${Dir}/MC/GJets_MG
mkdir -p ${eosPath}/${Dir}/MC/DiJet
mkdir -p ${eosPath}/${Dir}/MC/EWK
mkdir -p ${eosPath}/${Dir}/MC/Signal_Bstar
mkdir -p ${eosPath}/${Dir}/MC/Signal_Bstar/f1p0
mkdir -p ${eosPath}/${Dir}/MC/Signal_Bstar/f0p5
mkdir -p ${eosPath}/${Dir}/MC/Signal_Bstar/f0p1
mkdir -p ${eosPath}/${Dir}/MC/Signal_Qstar
mkdir -p ${eosPath}/${Dir}/MC/Signal_Qstar/f1p0
mkdir -p ${eosPath}/${Dir}/MC/Signal_Qstar/f0p5
mkdir -p ${eosPath}/${Dir}/MC/Signal_Qstar/f0p1


mv GJets_HT*.root ${eosPath}/${Dir}/MC/GJets_MG/
mv QCD_Pt*.root ${eosPath}/${Dir}/MC/DiJet/
mv BstarToGJ_M*_f1p0.root ${eosPath}/${Dir}/MC/Signal_Bstar/f1p0/
mv BstarToGJ_M*_f0p5.root ${eosPath}/${Dir}/MC/Signal_Bstar/f0p5/
mv BstarToGJ_M*_f0p1.root ${eosPath}/${Dir}/MC/Signal_Bstar/f0p1/
mv QstarToGJ_M*_f1p0.root ${eosPath}/${Dir}/MC/Signal_Qstar/f1p0/
mv QstarToGJ_M*_f0p5.root ${eosPath}/${Dir}/MC/Signal_Qstar/f0p5/
mv QstarToGJ_M*_f0p1.root ${eosPath}/${Dir}/MC/Signal_Qstar/f0p1/
mv W*.root ${eosPath}/${Dir}/MC/EWK/
mv ZGTo2LG.root ${eosPath}/${Dir}/MC/EWK/