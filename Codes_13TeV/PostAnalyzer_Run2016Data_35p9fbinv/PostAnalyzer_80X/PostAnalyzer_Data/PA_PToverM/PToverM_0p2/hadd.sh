#!/bin/bash

hadd singlePhotonRun2016Bv2.root singlePhotonRun2016Bv2_*.root
rm singlePhotonRun2016Bv2_*.root
hadd singlePhotonRun2016Cv1.root singlePhotonRun2016Cv1_*.root 
rm singlePhotonRun2016Cv1_*.root
hadd singlePhotonRun2016Dv1.root singlePhotonRun2016Dv1_*.root
rm singlePhotonRun2016Dv1_*.root
hadd singlePhotonRun2016Ev1.root singlePhotonRun2016Ev1_*.root
rm singlePhotonRun2016Ev1_*.root
hadd singlePhotonRun2016Fv1.root singlePhotonRun2016Fv1_*.root
rm singlePhotonRun2016Fv1_*.root
hadd singlePhotonRun2016Gv1.root singlePhotonRun2016Gv1_*.root
rm singlePhotonRun2016Gv1_*.root
hadd singlePhotonRun2016Hv2.root singlePhotonRun2016Hv2_*.root
rm singlePhotonRun2016Hv2_*.root
mv singlePhotonRun2016Hv3.root singlePhotonRun2016Hv3_*.root
rm singlePhotonRun2016Hv3_*.root

hadd singlePhotonRun2016Full.root singlePhotonRun2016*.root
mkdir Data
mv singlePhotonRun2016*.root Data/
