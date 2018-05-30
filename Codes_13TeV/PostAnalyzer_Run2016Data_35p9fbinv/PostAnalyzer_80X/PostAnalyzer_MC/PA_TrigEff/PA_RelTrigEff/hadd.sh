#!/bin/bash

hadd GJets_HT40To100.root GJets_HT40To100_*.root
rm GJets_HT40To100_*.root
hadd GJets_HT100To200.root GJets_HT100To200_*.root
rm GJets_HT100To200_*.root
hadd GJets_HT200To400.root GJets_HT200To400_*.root
rm GJets_HT200To400_*.root
hadd GJets_HT400To600.root GJets_HT400To600_*.root
rm GJets_HT400To600_*.root
hadd GJets_HT600ToInf.root GJets_HT600ToInf_*.root
rm GJets_HT600ToInf_*.root

hadd GJets_HT40ToInf.root GJets_HT*.root
mkdir GJ_madgraph
mv GJets_HT*.root GJ_madgraph/
