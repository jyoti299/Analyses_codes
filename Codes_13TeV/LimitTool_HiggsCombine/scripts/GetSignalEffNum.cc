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


void SignalQstar_eff(){

  TString model = "Qstar"; //"Qstar"
  TString coupling = "f1p0";

  Double_t Signal_mass[10] = {500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000};
  TString hist = "h_CutFlowWt_qstar";
 
  TString InPath = "../PostAna_Results/PA_Results_2015+2016/80X/Optimization/";
  TString Opti1 = "JetId_Opti/"; //DEta_Opti/  //DPhi_Opti/  //PhId_Opti/  //JetId_Opti/  //CSVDisc_Opti/  
  //TString Opti2[8] = {"DEta_1.0", "DEta_1.2", "DEta_1.5", "DEta_1.8", "DEta_2.0", "DEta_2.2", "DEta_2.5", "NoDEta"};
  //TString Opti2[4] = {"DPhi_1.5/", "DPhi_2.0/", "DPhi_2.5/", "NoDPhi/"};
  //TString Opti2[3] = {"PhId_Loose/", "PhId_Medium/", "PhId_Tight/"};
  TString Opti2[2] = {"JetId_Tight/", "JetId_TightLepVeto/"};
  //TString Opti2[3] = {"CSVL/", "CSVM/", "CSVT/"};
  
  for(Int_t ii = 0; ii < 2; ii++){

    std::vector<double> Mass_vec;
    std::vector<double> Eff_vec;

    Mass_vec.clear();
    Eff_vec.clear();
    for(int m = 0; m < 9; m++){
      
      std::ostringstream mass1;
      mass1 << Signal_mass[m];
      std::string mass1_str = mass1.str();
      std::ostringstream mass2;
      mass2 << Signal_mass[m+1];
      std::string mass2_str = mass2.str();

      TString InFile1 = InPath+Opti1+Opti2[ii]+"/MC/Signal_"+model+"/"+model+"ToGJ_M-"+mass1_str+"_"+coupling+".root";
      TString InFile2 = InPath+Opti1+Opti2[ii]+"/MC/Signal_"+model+"/"+model+"ToGJ_M-"+mass2_str+"_"+coupling+".root";

      TFile *f_sig1 = new TFile(InFile1, "READ");
      TH1F *h_sig1 = (TH1F*)f_sig1->Get(hist);

      TFile *f_sig2 = new TFile(InFile2, "READ");
      TH1F *h_sig2 = (TH1F*)f_sig2->Get(hist);

      //Efficiency for qstar before mass cut = deta/Total (10/1)                                                                        
      Double_t Eff1 = h_sig1->GetBinContent(10)/h_sig1->GetBinContent(1);
      Double_t Eff2 = h_sig2->GetBinContent(10)/h_sig2->GetBinContent(1);

      Double_t eff_frac, mass_frac;
      if(m == 0){
	eff_frac = (Eff2 - Eff1)/5;
        mass_frac = (Signal_mass[m+1] - Signal_mass[m])/5;
      }else{
	eff_frac = (Eff2 - Eff1)/10;
        mass_frac = (Signal_mass[m+1] - Signal_mass[m])/10;
      }

      Double_t cur_mass, cur_eff;
      if(m == 0){
	for(Int_t j = 0; j < 5; j++){
	  cur_mass = Signal_mass[m] + mass_frac*j;
	  cur_eff = Eff1 + eff_frac*j;

	  Mass_vec.push_back(cur_mass);
	  Eff_vec.push_back(cur_eff);
	}
      }else{
	for(Int_t j = 0; j < 10; j++){
	  cur_mass = Signal_mass[m] + mass_frac*j;
	  cur_eff = Eff1 + eff_frac*j;

	  Mass_vec.push_back(cur_mass);
	  Eff_vec.push_back(cur_eff);
	}
      }
    }
    
    //For last mass point
    std::ostringstream mass3;
    mass3 << Signal_mass[9];
    std::string mass3_str = mass3.str();
    TString InFile3 = InPath+Opti1+Opti2[ii]+"/MC/Signal_"+model+"/"+model+"ToGJ_M-"+mass3_str+"_"+coupling+".root";
    
    TFile *f_sig3 = new TFile(InFile3, "READ");
    TH1F *h_sig3 = (TH1F*)f_sig3->Get(hist);
 
    Double_t Eff3 = h_sig3->GetBinContent(10)/h_sig3->GetBinContent(1);
    Mass_vec.push_back(Signal_mass[9]);
    Eff_vec.push_back(Eff3);
    
    cout << "*************************************************************************" << endl;
    cout << "for " << Opti1 << " ===>>>>> " << Opti2[ii] << endl;
    cout << "*************************************************************************" << endl;
    cout << "" << endl;
    cout << "Mass    Eff" << endl;;
    for(unsigned int k = 0; k < Mass_vec.size(); k++){
      cout << Mass_vec[k] << "    " << Eff_vec[k] << endl;
    }

    cout << "" << endl;
    cout << "*************************************************************************" << endl;
    cout << "" << endl;
  }

}

void SignalBstar_eff(){

  TString model = "Bstar";
  TString coupling = "f1p0";

  Double_t Signal_mass[10] = {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
  TString hist = "h_CutFlowTotalWt_bstar";
 
  TString InPath = "../PostAna_Results/PA_Results_2015+2016/80X/Optimization/";
  TString Opti1 = "CSVDisc_Opti/"; //DEta_Opti/  //DPhi_Opti/  //PhId_Opti/  //JetId_Opti/  //CSVDisc_Opti/  
  //TString Opti2[8] = {"DEta_1.0", "DEta_1.2", "DEta_1.5", "DEta_1.8", "DEta_2.0", "DEta_2.2", "DEta_2.5", "NoDEta"};
  //TString Opti2[4] = {"DPhi_1.5/", "DPhi_2.0/", "DPhi_2.5/", "NoDPhi/"};
  //TString Opti2[3] = {"PhId_Loose/", "PhId_Medium/", "PhId_Tight/"};
  //TString Opti2[2] = {"JetId_Tight/", "JetId_TightLepVeto/"};
  TString Opti2[3] = {"CSVL/", "CSVM/", "CSVT/"};
  
  for(Int_t ii = 0; ii < 3; ii++){

    std::vector<double> Mass_vec;
    std::vector<double> Eff_vec;

    Mass_vec.clear();
    Eff_vec.clear();
    for(int m = 0; m < 9; m++){

      std::ostringstream mass1;
      mass1 << Signal_mass[m];
      std::string mass1_str = mass1.str();
      std::ostringstream mass2;
      mass2 << Signal_mass[m+1];
      std::string mass2_str = mass2.str();

      TString InFile1 = InPath+Opti1+Opti2[ii]+"/MC/Signal_"+model+"/"+model+"ToGJ_M-"+mass1_str+"_"+coupling+".root";
      TString InFile2 = InPath+Opti1+Opti2[ii]+"/MC/Signal_"+model+"/"+model+"ToGJ_M-"+mass2_str+"_"+coupling+".root";

      TFile *f_sig1 = new TFile(InFile1, "READ");
      TH1F *h_sig1 = (TH1F*)f_sig1->Get(hist);

      TFile *f_sig2 = new TFile(InFile2, "READ");
      TH1F *h_sig2 = (TH1F*)f_sig2->Get(hist);

      //Efficiency for qstar before mass cut = deta/Total (10/1)                                                                        
      Double_t Eff1 = h_sig1->GetBinContent(11)/h_sig1->GetBinContent(1);
      Double_t Eff2 = h_sig2->GetBinContent(11)/h_sig2->GetBinContent(1);

      Double_t eff_frac = (Eff2 - Eff1)/5;
      Double_t mass_frac = (Signal_mass[m+1] - Signal_mass[m])/5;

      for(Int_t j = 0; j < 5; j++){
        Double_t cur_mass = Signal_mass[m] + mass_frac*j;
        Double_t cur_eff = Eff1 + eff_frac*j;

        Mass_vec.push_back(cur_mass);
	Eff_vec.push_back(cur_eff);
      }
    }
    
    //For last mass point
    std::ostringstream mass3;
    mass3 << Signal_mass[9];
    std::string mass3_str = mass3.str();
    TString InFile3 = InPath+Opti1+Opti2[ii]+"/MC/Signal_"+model+"/"+model+"ToGJ_M-"+mass3_str+"_"+coupling+".root";
    
    TFile *f_sig3 = new TFile(InFile3, "READ");
    TH1F *h_sig3 = (TH1F*)f_sig3->Get(hist);
 
    Double_t Eff3 = h_sig3->GetBinContent(11)/h_sig3->GetBinContent(1);
    Mass_vec.push_back(Signal_mass[9]);
    Eff_vec.push_back(Eff3);
    
    cout << "*************************************************************************" << endl;
    cout << "for " << Opti1 << " ===>>>>> " << Opti2[ii] << endl;
    cout << "*************************************************************************" << endl;
    cout << "" << endl;
    cout << "Mass    Eff" << endl;;
    for(unsigned int k = 0; k < Mass_vec.size(); k++){
      cout << Mass_vec[k] << "    " << Eff_vec[k] << endl;
    }

    cout << "" << endl;
    cout << "*************************************************************************" << endl;
    cout << "" << endl;
  }

}






