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
#include "TRandom.h"
#include "TRandom3.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace ROOT;

void DataInputRootFile(){

  TString InPath = "../PostAna_Results/";
  //TString InDir = "PA_Results_2015+2016/80X/ReReco-BCDEFG_PromptReco-H/PhLID_JetTID_Pt190_NoDEta_NoDphi_CSVM_Mass695/";
  //TString InFile = "MC/Total_BkgMC/Total_BkgMC.root";
  //TString InFile = "Data/Run2016Full--23Sep2016_ReReco.root";
  TString InDir = "PA_Results_2015+2016/Re_miniAOD/PhMID_JetTID_Pt200_170_DEta1p5_NoDphi_CSVL_Mass695/";
  //TString InFile = "MC/Total_MC/TotalBkg_MC.root";
  TString InFile = "Data/singlePhotonRun2016Full.root";

  TString Out = "../inputs/";
  TString path = "Bias_study/";
  //  TString path = "ForGOF/";
  TString Outpath = Out+path;
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Json(Golden/Silver)Lumi(2p67fb or 24490pb)_boxName.root
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Lumi(2p67fb or 24490pb)_Cut_Opti_InvtMass.root

  //TString OutFile = "TotalMC_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root"; //To_cha 
  //TString OutFile = "PseudoData_FromMC_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root" 
  // TString OutFile = "PseudoData_FromMC_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_1tagBstarInvtMass.root";
  // TString OutFile = "Data_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_1tagBstarInvtMass.root"; 
  // TString OutFile = "Data_Run2016_ReReco-BCDEFG-PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_nodeta_nodphi_CSVM_Mass695_QstarInvtMass.root"; //To_change//for csv

  //TString OutFile = "TotalMC_Qstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_QstarInvtMass.root"; //for q*
  //TString OutFile = "PseudoData_FromMC_Qstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_1tagBstarInvtMass.root"; //for b*
  TString OutFile = "Data_Bstar_f-1p0_Spring16_80X_35866pb_Cut-PhMID_JetTID_Pt200_170_DEta1p5_noDPhi_CSVL_Mass695_1tagBstarInvtMass.root"; //for b*

  ostringstream InHist;
  //InHist << "h_GJetInvtMass_UnitBin_noMassCut"; //for q*
  InHist << "h_GbJetInvtMass_UnitBin_1BTag_noMassCut"; //for 1tag b*
  //  InHist << "h_GbJetInvtMass_UnitBin_0BTag_noMassCut"; //for 0tag b*

  //No Change below this
  ostringstream OutHist;
  OutHist << "h_mgj_13TeV_1GeVBin";

  TFile *fIn = new TFile(InPath+InDir+InFile, "READ");
  TFile *fOut = new TFile(Outpath+OutFile, "RECREATE");

  fIn->cd();
  TH1F *hIn = (TH1F*)fIn->Get(InHist.str().c_str());

  /*
  TH1F *hOut =  (TH1F*)hIn->Clone("hOut");
  hOut->Reset();

  for(Int_t i = 0; i < hIn->GetNbinsX(); i++){
    Double_t N = hIn->GetBinContent(i+1);
    hOut->SetBinContent(i+1, N);
    hOut->SetBinError(i+1, TMath::Sqrt(N));
  }
  */

  TH1F *hOut = new TH1F(OutHist.str().c_str(), OutHist.str().c_str(), 14000, 0, 14000);
  hOut->Sumw2();
   
  //  hOut->FillRandom(hIn);
  
  //For pseudo data
  TRandom3 *rand = new TRandom3(0);
  for(Int_t i = 0; i < hIn->GetNbinsX(); i++){
    Double_t N = hIn->GetBinContent(i+1);
    Double_t N_rdm = (double)rand->Poisson(N);
    Double_t Nerr = TMath::Sqrt(N_rdm);
    //cout << "N = " << N_rdm << ", N_err = " << Nerr << endl;
    hOut->SetBinContent(i+1, N_rdm);
    hOut->SetBinError(i+1, Nerr);
  }
  
  hOut->GetXaxis()->SetTitle("");
  hOut->GetYaxis()->SetTitle("");

  //hOut->SetName(OutHist.str().c_str());

  fOut->cd();
  hOut->Write();
  //hOut->Write(OutHist.str().c_str(), TObject::kWriteDelete);

  fOut->Close();

}


void DataInputRootFile_forOpti(){

  TString InPath = "../PostAna_Results/";
  TString InDir = "PA_Results_2015+2016/80X/Optimization/";
  TString OptiDir = "CSVDisc_Opti/"; //To_change
  TString Opti = "CSVT/"; //To_change
  TString InFile = "MC/Total_BkgMC/Total_BkgMC.root";

  TString Out = "../inputs/Optimization_80X/";
  TString Outpath = Out+OptiDir+Opti;
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Json(Golden/Silver)Lumi(2p67fb or 24490pb)_boxName.root
  //  TString OutFile = "Data_Run2015D-16DecReReco_76X_Silver2p67fb_ExcitedQuarks2015.root";
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Lumi(2p67fb or 24490pb)_Cut_Opti_InvtMass.root
  //  TString OutFile = "Data_Run2016BCDEFG-PromptReco_80X_24487pb_Cut-PhLID_JetTID_deta1.8_dphi1.5_CSVM_Mass695_PhIdOpti_QstarInvtMass.root"; //To_change //For deta, dphi, phid ,jetid
  TString OutFile = "Data_Run2016BCDEFG-PromptReco_80X_24487pb_Cut-PhLID_JetTID_deta1.8_dphi1.5_CSVT_Mass695_CSVOpti_1tagBstarInvtMass.root"; //To_change//for csv

  ostringstream InHist;
  //InHist << "h_GJetInvtMass_UnitBin_noMassCut"; //To_change // for deta, dphi, phid, jet id
  InHist << "h_GbJetInvtMass_UnitBin_1BTag_noMassCut"; //To_change //for csv

  //No Change below this
  ostringstream OutHist;
  OutHist << "h_mgj_13TeV_1GeVBin";

  TFile *fIn = new TFile(InPath+InDir+OptiDir+Opti+InFile, "READ");
  TFile *fOut = new TFile(Outpath+OutFile, "RECREATE");

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

void TotalMCInputRootFile(){

  TString InPath = "../PostAna_Results/";
  //  TString InDir = "PA_Results_2015+2016/80X/PhLID_JetTID_Pt190_Nodeta_Nodphi_NoMassCut_CSVM/";
  TString InDir = "PA_Results_2015+2016/80X/ReReco-BCDEFG_PromptReco-H/PhLID_JetTID_Pt190_NoDEta_NoDphi_CSVM_Mass695/";
  //  TString InFile = "Total_MC/Total_MC.root";
  TString InFile = "MC/Total_BkgMC_WithJPSF/Total_BkgMC.root";

  TString Out = "../inputs/Bias_study/";
  TString Outpath = Out;
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Json(Golden/Silver)Lumi(2p67fb or 24490pb)_boxName.root
  //  TString OutFile = "Data_Run2015D-16DecReReco_76X_Silver2p67fb_ExcitedQuarks2015.root";
  ///Output root file name: Data_RunEra(like Run2016BCDEFG)_Prompt/ReReco_CMSSWversion_Lumi(2p67fb or 24490pb)_Cut_Opti_InvtMass.root
  //  TString OutFile = "Data_Run2016BCDEFG-PromptReco_80X_24487pb_Cut-PhLID_JetTID_deta1.8_dphi1.5_CSVM_Mass695_PhIdOpti_QstarInvtMass.root"; //To_change //For deta, dphi, phid ,jetid
  //  TString OutFile = "TotalMC_Run2016_ReReco-BCDEFG_PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_Nodeta_Nodphi_CSVM_Mass695_QstarInvtMass.root"; //To_cha
  TString OutFile = "TotalMC_WithJPSF_Run2016_ReReco-BCDEFG_PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_Nodeta_Nodphi_CSVM_Mass695_QstarInvtMass.root"; //To_change//for csv

  ostringstream InHist;
  InHist << "h_GJetInvtMass_UnitBin_noMassCut";

  //No Change below this
  ostringstream OutHist;
  OutHist << "h_mgj_13TeV_1GeVBin";

  TFile *fIn = new TFile(InPath+InDir+InFile, "READ");
  TFile *fOut = new TFile(Outpath+OutFile, "RECREATE");

  fIn->cd();

  TH1F *hIn = (TH1F*)fIn->Get(InHist.str().c_str());
  //  hIn->Scale(36.813/24.487);

  TH1F *hOut =  (TH1F*)hIn->Clone("hOut");
  hOut->Reset();

  for(Int_t i = 0; i < hIn->GetNbinsX(); i++){
    Double_t N = hIn->GetBinContent(i+1);
    hOut->SetBinContent(i+1, N);
  }

  hOut->GetXaxis()->SetTitle("");
  hOut->GetYaxis()->SetTitle("");

  hOut->SetName(OutHist.str().c_str());

  fOut->cd();
  hOut->Write(OutHist.str().c_str(), TObject::kWriteDelete);

  fOut->Close();

}

