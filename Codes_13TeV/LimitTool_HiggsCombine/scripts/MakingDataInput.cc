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
//#include "TMinuit.h"                           
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace ROOT;

void DataInputRootFile(){

  TString InPath = "../PostAna_Results/";
  TString InDir = "PA_Results_2015+2016/Re_miniAOD/PhMID_JetTID_Pt200_170_DEta1p5_NoDphi_CSVL_Mass695/";
  TString InFile = "MC/Total_MC/TotalBkg_MC.root";
  //TString InFile = "Data/singlePhotonRun2016Full.root";

  //  TString InDir = "PA_Results_2015+2016/80X/ReReco-BCDEFG_PromptReco-H/PhLID_JetTID_Pt190_NoDEta_NoDphi_CSVM_Mass695/";
  //  TString InFile = "Data/Run2016Full--23Sep2016_ReReco.root";

  TString Out = "../inputs/";
  TString Outpath = Out;
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Json(Golden/Silver)Lumi(2p67fb or 24490pb)_boxName.root
  //TString OutFile = "Data_RunBtoH_rereco_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root"; //for q*
  //TString OutFile = "Data_ReminiAOD_80X_35866pb_Cut-PhLID_JetTID_deta1p5_nodphi_CSVL_Mass695_QstarInvtMass.root"; //for 1tag b*
  TString OutFile = "TotalMC_ReminiAOD_80X_35866pb_Cut-PhLID_JetTID_deta1p5_nodphi_CSVL_Mass695_QstarInvtMass.root"; //for 1tag b*
  //  TString OutFile = "Data_RunBtoH_rereco_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_0tagBstarInvtMass.root"; //for 0tag b*

  ostringstream InHist;
  InHist << "h_GJetInvtMass_UnitBin_noMassCut"; //for q*
  //InHist << "h_GbJetInvtMass_UnitBin_1BTag_noMassCut"; //for 1tag b*
  //InHist << "h_GbJetInvtMass_UnitBin_0BTag_noMassCut"; //for 0tag b*

  //No Change below this
  ostringstream OutHist;
  OutHist << "h_mgj_13TeV_1GeVBin"; // Do not change this name as it is used in python scripts.

  TFile *fIn = new TFile(InPath+InDir+InFile, "READ");
  TFile *fOut = new TFile(Outpath+OutFile, "RECREATE");

  fIn->cd();

  TH1F *hIn = (TH1F*)fIn->Get(InHist.str().c_str());

  TH1F *hOut =  (TH1F*)hIn->Clone("hOut");
  hOut->Reset();

  for(Int_t i = 0; i < hIn->GetNbinsX(); i++){
    Double_t N = hIn->GetBinContent(i+1);
    Double_t N_err = hIn->GetBinError(i+1);
    hOut->SetBinContent(i+1, N);
    hOut->SetBinError(i+1, 0);
  }

  hOut->GetXaxis()->SetTitle("");
  hOut->GetYaxis()->SetTitle("");

  hOut->SetName(OutHist.str().c_str());

  fOut->cd();
  hOut->Write(OutHist.str().c_str(), TObject::kWriteDelete);

  fOut->Close();

}

