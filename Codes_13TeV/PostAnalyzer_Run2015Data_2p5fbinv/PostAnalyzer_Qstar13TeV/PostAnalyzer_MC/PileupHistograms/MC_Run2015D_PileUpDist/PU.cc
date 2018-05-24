#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
//#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
//#include "TMinuit.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TArrow.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>


void pu_GJets(){

  //The hPUTrue distribution is present in the Ntuples.
  //Here this input file is hadd of all the files present in GJets_HT100To200 folder, i did the hadd and saved this file in check dir and then delete it after getting the pileup histogram to save extra space taken unnecessarily by the hadded file. considered the pu dist for one sample only as it will be same for all the samples. 
  TFile *fin = new TFile("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/GJets/check/AOD_GJets_HT100To200.root", "READ");
  TFile *fout = new TFile("GJets_PileupHist.root", "RECREATE");

  fin->cd();
  TDirectory *dir = fin->GetDirectory("ggNtuplizer");

  TH1F *h = (TH1F*)dir->Get("hPUTrue");

  fout->cd();
  h->Write();
  fout->Write();
  fout->Close();

}

void pu_DiJet(){

  //Here this input file is hadd of all the files present in QCD_Pt_300to470 folder, i did the hadd and saved this file in check dir and then delete it after getting the pileup histogram to save extra space taken unnecessarily by the hadded file. considered the pu dist for one sample only as it will be same for all the samples. 
  TFile *fin = new TFile("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/check/AOD_QCD_Pt_300to470.root", "READ");
  TFile *fout = new TFile("DiJet_PileupHist.root", "RECREATE");

  fin->cd();
  TDirectory *dir = fin->GetDirectory("ggNtuplizer");

  TH1F *h = (TH1F*)dir->Get("hPUTrue");

  fout->cd();
  h->Write();
  fout->Write();
  fout->Close();

}

void pu_Qstar(){

  //Here this input file is hadd of all the files present in QstarToGJ_M1000_f1p0 folder, i did the hadd and saved this file in check dir and then delete it after getting the pileup histogram to save extra space taken unnecessarily by the hadded file. considered the pu dist for one sample only as it will be same for all the samples.
  TFile *fin = new TFile("/eos/uscms/store/user/lpcqstar/13TeV/NTuples/MC/check/AOD_QstarToGJ_M1000_f1p0.root", "READ");
  TFile *fout = new TFile("Qstar_PileupHist.root", "RECREATE");

  fin->cd();
  TDirectory *dir = fin->GetDirectory("ggNtuplizer");

  TH1F *h = (TH1F*)dir->Get("hPUTrue");

  fout->cd();
  h->Write();
  fout->Write();
  fout->Close();

}
