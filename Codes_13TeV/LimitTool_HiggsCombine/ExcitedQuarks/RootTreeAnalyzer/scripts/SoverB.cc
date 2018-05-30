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

void Qstar_SoverB(){

  TString InPath80 = "/eos/uscms/store/user/rocky86/13TeV/PA_Results_2015+2016/80X/Optimization/";

  TString Opti1 = "JetId_Opti/"; //"DEta_Opti/",  "DPhi_Opti/",  "JetId_Opti/",  "PhId_Opti/" 
  //TString Opti2[8] = {"DEta_1.0", "DEta_1.2", "DEta_1.5", "DEta_1.8", "DEta_2.0", "DEta_2.2", "DEta_2.5", "NoDEta"};
  //TString Opti2[4] = {"DPhi_1.5",  "DPhi_2.0",  "DPhi_2.5",  "NoDPhi"}; 
  //TString Opti2[3] = {"PhId_Loose",  "PhId_Medium",  "PhId_Tight"};     
  TString Opti2[2] = {"JetId_Tight",  "JetId_TightLepVeto"}; 

  for(Int_t i = 0; i < 2; i++){

      TFile *f_Total = new TFile(InPath80+Opti1+Opti2[i]+"/MC/Total_BkgMC/Total_BkgMC.root", "READ");
      TH1F *h_Total = (TH1F*)f_Total->Get("h_CutFlowWt_qstar");

      TFile *fsig = new TFile(InPath80+Opti1+Opti2[i]+"/MC/Signal_Qstar/QstarToGJ_M-3000_f1p0.root", "READ");
      TH1F *hsig = (TH1F*)fsig->Get("h_CutFlowWt_qstar");

      cout << "---------------------------------------------------------------------------------------" << endl;

      Double_t binValue_Total, binValue_sig;
      std::string binName = h_Total->GetXaxis()->GetBinLabel(11);
      binValue_Total = h_Total->GetBinContent(11);
      binValue_sig = hsig->GetBinContent(11);
      Double_t sb = binValue_sig/sqrt(binValue_Total);

      cout << "| For " << Opti2[i]  << "    |    " << binName << "    |    Total = "  <<  binValue_Total << "    |   sig = " << binValue_sig << "   |  s/sqrt(b) = " << sb << "   |  "<< endl;
      cout << "" << endl;

    }
  cout << "---------------------------------------------------------------------------------------" << endl;

}


void Bstar_SoverB(){

  TString InPath80 = "/eos/uscms/store/user/rocky86/13TeV/PA_Results_2015+2016/80X/Optimization/";

  TString Opti1 =  "CSVDisc_Opti/";
  TString Opti2[3] = {"CSVL", "CSVM", "CSVT"};

  for(Int_t i = 0; i < 3; i++){

    TFile *f_Total = new TFile(InPath80+Opti1+Opti2[i]+"/MC/Total_BkgMC/Total_BkgMC.root", "READ");
    TH1F *h_Total = (TH1F*)f_Total->Get("h_CutFlowTotalWt_bstar");

    TFile *fsig = new TFile(InPath80+Opti1+Opti2[i]+"/MC/Signal_Bstar/BstarToGJ_M-2000_f1p0.root", "READ");
    TH1F *hsig = (TH1F*)fsig->Get("h_CutFlowTotalWt_bstar");

    cout << "---------------------------------------------------------------------------------------" << endl;

    Double_t binValue_Total, binValue_sig;
    std::string binName = h_Total->GetXaxis()->GetBinLabel(12);
    binValue_Total = h_Total->GetBinContent(12);
    binValue_sig = hsig->GetBinContent(12);
    Double_t sb = binValue_sig/sqrt(binValue_Total);

    cout << "| For " << Opti2[i]  << "    |    " << binName << "    |    Total = "  <<  binValue_Total << "    |   sig = " << binValue_sig << "   |  s/sqrt(b) = " << sb << "   |  "<< endl;
    cout << "" << endl;

  }
  cout << "---------------------------------------------------------------------------------------" << endl;

}
