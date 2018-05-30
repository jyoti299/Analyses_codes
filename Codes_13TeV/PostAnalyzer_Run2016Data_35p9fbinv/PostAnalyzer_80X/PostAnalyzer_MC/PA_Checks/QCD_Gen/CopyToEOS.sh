#!/bin/bash

#rm condor* Job_* analysis_* *exe PostAnalyzer_MC.*

eosPath=/eos/uscms/store/user/rocky86/13TeV/PA_Results_2015+2016/80X
Dir=Checks/QCD_GenDists

mkdir -p ${eosPath}/${Dir}/MC/DiJet/

mv QCD_Pt*.root ${eosPath}/${Dir}/MC/DiJet/
