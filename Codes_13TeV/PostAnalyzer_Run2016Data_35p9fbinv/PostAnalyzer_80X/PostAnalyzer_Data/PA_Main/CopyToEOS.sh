#!/bin/bash

rm condor* Job_* analysis_* *exe PostAnalyzer_Data.*

eosPath=/eos/uscms/store/user/rocky86/13TeV/PA_Results_2015+2016/80X
Dir=PhLID_JetTID_Pt190_NoDEta_NoDPhi_NoMassCut


mkdir -p ${eosPath}/${Dir}
mkdir -p ${eosPath}/${Dir}/Data

xrdcp Run2016*.root ${eosPath}/${Dir}/Data/