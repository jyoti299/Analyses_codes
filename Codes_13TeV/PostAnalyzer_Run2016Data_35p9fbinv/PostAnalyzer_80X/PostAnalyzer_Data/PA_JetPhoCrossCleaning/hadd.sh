#!/bin/bash

hadd Run2016B-23Sep2016_ReReco.root  Run2016B-23Sep2016_ReReco-v3_*.root
rm Run2016B-23Sep2016_ReReco-v3_*.root
hadd Run2016C-23Sep2016_ReReco.root Run2016C-23Sep2016_ReReco-v1_*.root
rm Run2016C-23Sep2016_ReReco-v1_*.root
hadd Run2016D-23Sep2016_ReReco.root Run2016D-23Sep2016_ReReco-v1_*.root
rm Run2016D-23Sep2016_ReReco-v1_*.root
hadd Run2016E-23Sep2016_ReReco.root  Run2016E-23Sep2016_ReReco-v1_*.root
rm Run2016E-23Sep2016_ReReco-v1_*.root
hadd Run2016F-23Sep2016_ReReco.root  Run2016F-23Sep2016_ReReco-v1_*.root
rm Run2016F-23Sep2016_ReReco-v1_*.root
hadd Run2016G-23Sep2016_ReReco.root Run2016G-23Sep2016_ReReco-v1_*.root
rm Run2016G-23Sep2016_ReReco-v1_*.root
hadd Run2016H-PromptReco-v2.root Run2016H-PromptReco-v2_*.root
rm Run2016H-PromptReco-v2_*.root
hadd Run2016H-PromptReco-v3.root Run2016H-PromptReco-v3_*.root
rm Run2016H-PromptReco-v3_*.root

hadd Run2016Full-23Sep2016_ReReco.root Run2016*.root
mkdir Data
mv Run2016*.root Data/
