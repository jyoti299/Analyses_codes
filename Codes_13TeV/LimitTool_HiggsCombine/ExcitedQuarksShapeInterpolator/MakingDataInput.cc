#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h" 
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMinuit.h"                           
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include <boost/filesystem.hpp>

using namespace std;
using namespace ROOT;


void FromMCRandom_forOpti(){

  const int n_opti = 8;
  TString Opti[n_opti] = {"DEta_1.0", "DEta_1.2", "DEta_1.5", "DEta_1.8", "DEta_2.0", "DEta_2.2", "DEta_2.5", "NoDEta"};
 
  for(unsigned int ip = 0; ip < n_opti; ip++){

    cout << "Running for " << Opti[ip] << endl;

    ostringstream InPath, OutPath;
    InPath.str("");  OutPath.str("");
    InPath.clear();  OutPath.clear();
    InPath << "PA_Results/Spring16_36813pb_80X/Optimization/" << Opti[ip] << "/MC/";
    OutPath << "DataInvtMassDists/Spring16_36813pb_80X/Optimization/" << Opti[ip] << "/";
    //Making dir to store output files
    boost::filesystem::create_directories(OutPath.str().c_str());

    ostringstream OutFile;
    ostringstream InFile;
    InFile.str("");  OutFile.str("");
    InFile.clear();  OutFile.clear();
    InFile << InPath.str() << "Total_BkgMC/Total_BkgMC.root";
    OutFile << OutPath.str() << "PseudoData_80X_36813pb_Cut-PhLID_JetTID_nodphi_CSVL_Mass695_QstarInvtMass_Spring16.root";

    ostringstream InHist;
    ostringstream OutHist;
    InHist.str("");  OutHist.str("");
    InHist.clear();  OutHist.clear();
    InHist << "h_GJetInvtMass_UnitBin_noMassCut";
    OutHist << "h_mgj_13TeV_1GeVBin";

    TFile *fIn = new TFile(InFile.str().c_str(), "READ");
    TFile *fOut = new TFile(OutFile.str().c_str(), "RECREATE");

    fIn->cd();
    TH1F *hIn = (TH1F*)fIn->Get(InHist.str().c_str());

    TH1F *hOut =  (TH1F*)hIn->Clone("hOut");
    hOut->Reset();

    TRandom *rand = new TRandom(0);
    for(Int_t i = 0; i < hIn->GetNbinsX(); i++){
      Double_t N = hIn->GetBinContent(i+1);
      Double_t N_rdm = (double)rand->Poisson(N);
      hOut->SetBinContent(i+1, N_rdm);
    }

    hOut->GetXaxis()->SetTitle("");
    hOut->GetYaxis()->SetTitle("");

    hOut->SetName(OutHist.str().c_str());

    fOut->cd();
    hOut->Write(OutHist.str().c_str(), TObject::kWriteDelete);

    fOut->Close();
  }

}

void DataInputRootFile(){

  //TString InDir = "PA_Results/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/";
  TString InDir = "PA_Results/PToverM/NoDeta_NoPToverM/";
  TString InFile = "Data/singlePhotonRun2016Full.root";

  //TString Out = "DataInvtMassDists/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/";
  TString Out = "DataInvtMassDists/PToverM/NoDeta_NoPToverM/";
  TString Outpath = Out;

  //boost::filesystem::create_directories("DataInvtMassDists/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/");
  boost::filesystem::create_directories("DataInvtMassDists/PToverM/NoDeta_NoPToverM/");

  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Json(Golden/Silver)Lumi(2p67fb or 24490pb)_boxName.root
  //TString OutFile = "Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Mass700_QstarInvtMass.root"; //for q*
  //TString OutFile = "Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Mass700_1tagBstarInvtMass.root"; //for 1tag b*
  //TString OutFile = "Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Mass700_0tagBstarInvtMass.root"; //for 0tag b*
  TString OutFile = "Data_ReminiAOD_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_NoPTM_NoDEta_NoDPhi_CSVL_Mass700_QstarInvtMass.root"; //for q*
  
  ostringstream InHist;
  InHist << "h_GJetInvtMass_UnitBin_noMassCut"; //for q*
  //InHist << "h_GbJetInvtMass_UnitBin_1BTag_noMassCut"; //for 1tag b*
  //InHist << "h_GbJetInvtMass_UnitBin_0BTag_noMassCut"; //for 0tag b*

  //No Change below this
  ostringstream OutHist;
  OutHist << "h_mgj_13TeV_1GeVBin"; // Do not change this name as it is used in python scripts.

  TFile *fIn = new TFile(InDir+InFile, "READ");
  TFile *fOut = new TFile(Outpath+OutFile, "RECREATE");
  //TFile *fIn = new TFile("singlePhotonRun2016Full.root", "READ");
  //TFile *fOut = new TFile("GOF_check.root", "RECREATE");

  fIn->cd();

  TH1F *hIn = (TH1F*)fIn->Get(InHist.str().c_str());

  TH1F *hOut =  (TH1F*)hIn->Clone("hOut");
  hOut->Reset();

  for(Int_t i = 0; i < hIn->GetNbinsX(); i++){
    Double_t N = hIn->GetBinContent(i+1);
    hOut->SetBinContent(i+1, N);
    hOut->SetBinError(i+1, TMath::Sqrt(N));
  }

  hOut->GetXaxis()->SetTitle("");
  hOut->GetYaxis()->SetTitle("");

  hOut->SetName(OutHist.str().c_str());

  fOut->cd();
  hOut->Write(OutHist.str().c_str(), TObject::kWriteDelete);

  fOut->Close();

}



