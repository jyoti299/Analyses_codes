#!/bin/bash

hadd JetHT_Run2016B_23Sep2016v3.root JetHT_Run2016B_23Sep2016v3_*.root
rm JetHT_Run2016B_23Sep2016v3_*.root
hadd JetHT_Run2016C_23Sep2016v1.root JetHT_Run2016C_23Sep2016v1_*.root
rm JetHT_Run2016C_23Sep2016v1_*.root
hadd JetHT_Run2016D_23Sep2016v1.root JetHT_Run2016D_23Sep2016v1_*.root
rm JetHT_Run2016D_23Sep2016v1_*.root
hadd JetHT_Run2016E_23Sep2016v1.root JetHT_Run2016E_23Sep2016v1_*.root
rm JetHT_Run2016E_23Sep2016v1_*.root
hadd JetHT_Run2016F_23Sep2016v1.root JetHT_Run2016F_23Sep2016v1_*.root
rm JetHT_Run2016F_23Sep2016v1_*.root
hadd JetHT_Run2016G_23Sep2016v1.root JetHT_Run2016G_23Sep2016v1_*.root
rm JetHT_Run2016G_23Sep2016v1_*.root
hadd JetHT_Run2016H_PromptRecov2.root JetHT_Run2016H_PromptRecov2_*.root
rm JetHT_Run2016H_PromptRecov2_*.root
mv JetHT_Run2016H_PromptRecov3_1_25.root JetHT_Run2016H_PromptRecov3.root

hadd JetHT_Run2016Full_23Sep2016.root JetHT_Run2016*.root

mkdir Data_JetHT
mv JetHT_Run2016*.root Data_JetHT/