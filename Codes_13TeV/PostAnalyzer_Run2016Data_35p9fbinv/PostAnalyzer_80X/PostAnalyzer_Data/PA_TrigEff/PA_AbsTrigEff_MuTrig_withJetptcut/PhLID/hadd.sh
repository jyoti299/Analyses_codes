#!/bin/bash

hadd SingleMuonRun2016Bv2.root SingleMuonRun2016Bv2_*.root
rm SingleMuonRun2016Bv2_*.root
hadd SingleMuonRun2016Cv1.root SingleMuonRun2016Cv1_*.root 
rm SingleMuonRun2016Cv1_*.root
hadd SingleMuonRun2016Dv1.root SingleMuonRun2016Dv1_*.root
rm SingleMuonRun2016Dv1_*.root
hadd SingleMuonRun2016Ev1.root SingleMuonRun2016Ev1_*.root
rm SingleMuonRun2016Ev1_*.root
hadd SingleMuonRun2016Fv1.root SingleMuonRun2016Fv1_*.root
rm SingleMuonRun2016Fv1_*.root
hadd SingleMuonRun2016Gv1.root SingleMuonRun2016Gv1_*.root
rm SingleMuonRun2016Gv1_*.root
hadd SingleMuonRun2016Hv2.root SingleMuonRun2016Hv2_*.root
rm SingleMuonRun2016Hv2_*.root
mv SingleMuonRun2016Hv3_1_25.root SingleMuonRun2016Hv3.root

hadd singleMuRun2016Full.root SingleMuonRun2016*.root
mkdir Data
mv *.root Data/
