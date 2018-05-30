#!/bin/bash

hadd singleMuRun2016Bv2.root singleMuRun2016Bv2_*.root
rm singleMuRun2016Bv2_*.root
hadd singleMuRun2016Cv1.root singleMuRun2016Cv1_*.root 
rm singleMuRun2016Cv1_*.root
hadd singleMuRun2016Dv1.root singleMuRun2016Dv1_*.root
rm singleMuRun2016Dv1_*.root
hadd singleMuRun2016Ev1.root singleMuRun2016Ev1_*.root
rm singleMuRun2016Ev1_*.root
hadd singleMuRun2016Fv1.root singleMuRun2016Fv1_*.root
rm singleMuRun2016Fv1_*.root
hadd singleMuRun2016Gv1.root singleMuRun2016Gv1_*.root
rm singleMuRun2016Gv1_*.root
hadd singleMuRun2016Hv2.root singleMuRun2016Hv2_*.root
rm singleMuRun2016Hv2_*.root
mv singleMuRun2016Hv3_1_25.root singleMuRun2016Hv3.root

hadd singleMuRun2016Full.root singleMuRun2016*.root
mkdir Data
mv singleMuRun2016*.root Data/
