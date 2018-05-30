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

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include <boost/filesystem.hpp>

using namespace std;
using namespace ROOT;

void QstarXdists_oneFile(){

  TString model = "Qstar";
  TString coupling = "f-1p0";

  const int n = 8;
  Double_t Mass[n] = {500, 1000, 2000, 3000, 4000, 5000, 6000, 7000};//, 8000, 9000};

  TString InPath = "PA_Results/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866Pb/ResShapes_n_Syst/";

  TString OutPath = "MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/";
  //Making dir to store output files
  boost::filesystem::create_directories("MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_35866Pb/");

  ostringstream OutFile;
  ostringstream InHist;

  InHist.str("");  OutFile.str("");
  InHist.clear();  OutFile.clear();
  InHist << "h_mass_X_bin1";
  OutFile << OutPath << model << "_XDists_" << coupling << ".root";

  TFile *fout = new TFile(OutFile.str().c_str(), "RECREATE");

  for(unsigned int im = 0; im < n; im++){

    ostringstream InFile;
    InFile.str("");
    InFile.clear();
    InFile << InPath << model << "/" << model << "ToGJ_M-" << Mass[im]   << "_" << coupling << ".root";

    TFile *f_sig = new TFile(InFile.str().c_str(), "READ");
    TH1F *h_sig = (TH1F*)f_sig->Get(InHist.str().c_str());

    h_sig->Rebin(10);

    ostringstream OutHist;
    OutHist << "h_" << model << "ToGJ_" << "_M" << Mass[im] << "_XDists";

    h_sig->GetXaxis()->SetTitle("");
    h_sig->GetYaxis()->SetTitle("");
    h_sig->SetName(OutHist.str().c_str());

    fout->cd();
    h_sig->Write(OutHist.str().c_str(), TObject::kWriteDelete);
  }
  fout->Close();
}





void QstarXdists(){

  TString model = "Qstar";
  TString coupling = "f1p0";

  const int n_syst = 7;
  TString SystHist[n_syst] = {"", "_JES_up", "_JES_down", "_PES_up", "_PES_down", "_JER", "_PER"};
  TString SystFile[n_syst] = {"", "_JESUP", "_JESDOWN", "_PESUP", "_PESDOWN", "_JER", "_PER"};

  TString Era = "Summer16";

  const int n = 8;
  Double_t Mass[n] = {500, 1000, 2000, 3000, 4000, 5000, 6000, 7000};//, 8000, 9000};

  TString InPath = "PA_Results/PToverM/NoDeta_NoPToverM/MC/";
  //TString InPath = "PA_Results/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/ResShapes_n_Syst/";

  //TString OutPath = "MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/";
  TString OutPath = "MassXDist/PToverM/NoDeta_NoPToverM/";
  //Making dir to store output files
  //boost::filesystem::create_directories("MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/");
  boost::filesystem::create_directories("MassXDist/PToverM/NoDeta_NoPToverM/");

  ostringstream OutFile;
  ostringstream InHist;

  for(unsigned int is = 0; is < n_syst; is++){

    InHist.str("");  OutFile.str("");
    InHist.clear();  OutFile.clear();
    InHist << "h_mass_X_bin1" << SystHist[is];
    OutFile << OutPath << model << "_XDists_" << coupling << SystFile[is] << ".root";

    TFile *fout = new TFile(OutFile.str().c_str(), "RECREATE");

    for(unsigned int im = 0; im < n; im++){

      ostringstream InFile;
      InFile.str("");
      InFile.clear();
      InFile << InPath << model << "/" << model << "ToGJ_M" << Mass[im]   << "_" << coupling << ".root";

      TFile *f_sig = new TFile(InFile.str().c_str(), "READ");
      TH1F *h_sig = (TH1F*)f_sig->Get(InHist.str().c_str());

      h_sig->Rebin(10);

      ostringstream OutHist;
      OutHist << "h_" << model << "ToGJ_" << Era << "_M" << Mass[im] << "_XDists";

      h_sig->GetXaxis()->SetTitle("");
      h_sig->GetYaxis()->SetTitle("");
      h_sig->SetName(OutHist.str().c_str());

      fout->cd();
      h_sig->Write(OutHist.str().c_str(), TObject::kWriteDelete);
    }
    fout->Close();
  }

}


void BstarXdists(){

  TString model = "Bstar";
  TString coupling = "f0p1";
  TString tag = "1bTag";

  const int n_syst = 9;
  TString SystHist[n_syst] = {"", "_JES_up_1bTag", "_JES_down_1bTag", "_PES_up_1bTag", "_PES_down_1bTag", "_JER_1bTag", "_PER_1bTag", "_BSF_1bTag_up", "_BSF_1bTag_down"};
  TString SystFile[n_syst] = {"", "_JESUP", "_JESDOWN", "_PESUP", "_PESDOWN", "_JER", "_PER", "_BSFUP", "_BSFDOWN"};

  TString Era = "Summer16";

  const int n = 10;
  Double_t Mass[n] = {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};

  TString InPath = "PA_Results/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/ResShapes_n_Syst/";

  TString OutPath = "MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/";
  //Making dir to store output files
  boost::filesystem::create_directories("MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/");

  ostringstream OutFile;
  ostringstream InHist;

  for(unsigned int is = 0; is < n_syst; is++){

    InHist.str("");  OutFile.str("");
    InHist.clear();  OutFile.clear();
    InHist << "h_mass_X_bin1" << SystHist[is];
    OutFile << OutPath << model << "_" << tag << "_XDists_" << coupling << SystFile[is] << ".root";

    TFile *fout = new TFile(OutFile.str().c_str(), "RECREATE");

    for(unsigned int im = 0; im < n; im++){

      ostringstream InFile;
      InFile.str("");
      InFile.clear();
      InFile << InPath << model << "/" << model << "ToGJ_M" << Mass[im]   << "_" << coupling << ".root";

      TFile *f_sig = new TFile(InFile.str().c_str(), "READ");
      TH1F *h_sig = (TH1F*)f_sig->Get(InHist.str().c_str());

      h_sig->Rebin(10);

      ostringstream OutHist;
      OutHist << "h_" << model << "ToGJ_" << Era << "_M" << Mass[im] << "_XDists";

      h_sig->GetXaxis()->SetTitle("");
      h_sig->GetYaxis()->SetTitle("");
      h_sig->SetName(OutHist.str().c_str());

      fout->cd();
      h_sig->Write(OutHist.str().c_str(), TObject::kWriteDelete);
    }
    fout->Close();
  }
}

void BstarXdists_0tag(){

  TString model = "Bstar";
  TString coupling = "f0p1";
  TString tag = "0bTag";

  const int n_syst = 9;
  TString SystHist[n_syst] = {"", "_JES_up_0bTag", "_JES_down_0bTag", "_PES_up_0bTag", "_PES_down_0bTag", "_JER_0bTag", "_PER_0bTag", "_BSF_0bTag_up", "_BSF_0bTag_down"};
  TString SystFile[n_syst] = {"", "_JESUP", "_JESDOWN", "_PESUP", "_PESDOWN", "_JER", "_PER", "_BSFUP", "_BSFDOWN"};

  TString Era = "Summer16";

  const int n = 10;
  Double_t Mass[n] = {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};

  TString InPath = "PA_Results/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/ResShapes_n_Syst/";

  TString OutPath = "MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/";
  //Making dir to store output files
  boost::filesystem::create_directories("MassXDist/PhMID_JetTID_Pt200_170_DEta1p5_NoDPhi_CSVL_Summer16_35866pb/");

  ostringstream OutFile;
  ostringstream InHist;

  for(unsigned int is = 0; is < n_syst; is++){

    InHist.str("");  OutFile.str("");
    InHist.clear();  OutFile.clear();
    InHist << "h_mass_X_bin1" << SystHist[is];
    OutFile << OutPath << model << "_" << tag << "_XDists_" << coupling << SystFile[is] << ".root";

    TFile *fout = new TFile(OutFile.str().c_str(), "RECREATE");

    for(unsigned int im = 0; im < n; im++){

      ostringstream InFile;
      InFile.str("");
      InFile.clear();
      InFile << InPath << model << "/" << model << "ToGJ_M" << Mass[im]   << "_" << coupling << ".root";

      TFile *f_sig = new TFile(InFile.str().c_str(), "READ");
      TH1F *h_sig = (TH1F*)f_sig->Get(InHist.str().c_str());

      h_sig->Rebin(10);

      ostringstream OutHist;
      OutHist << "h_" << model << "ToGJ_" << Era << "_M" << Mass[im] << "_XDists";

      h_sig->GetXaxis()->SetTitle("");
      h_sig->GetYaxis()->SetTitle("");
      h_sig->SetName(OutHist.str().c_str());

      fout->cd();
      h_sig->Write(OutHist.str().c_str(), TObject::kWriteDelete);
    }
    fout->Close();
  }
}


void QstarXdists_forOpti(){

  TString model = "Qstar";
  TString coupling = "f-1p0";

  const int n_opti = 8;
  TString Opti[n_opti] = {"DEta_1.0", "DEta_1.2", "DEta_1.5", "DEta_1.8", "DEta_2.0", "DEta_2.2", "DEta_2.5", "NoDEta"};
 
  const int n_syst = 7;
  TString SystHist[n_syst] = {"", "_JES_up", "_JES_down", "_PES_up", "_PES_down", "_JER", "_PER"};
  TString SystFile[n_syst] = {"", "_JESUP", "_JESDOWN", "_PESUP", "_PESDOWN", "_JER", "_PER"};

  TString Era = "Spring16";

  const int n = 10;
  Double_t Mass[n] = {500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000};

  for(unsigned int ip = 0; ip < n_opti; ip++){

    cout << "Running for " << Opti[ip] << endl;
    ostringstream InPath, OutPath;
    InPath.str("");  OutPath.str("");
    InPath.clear();  OutPath.clear();
    InPath << "PA_Results/Spring16_36813pb_80X/Optimization/" << Opti[ip] << "/MC/";
    OutPath << "MassXDist/Spring16_36813pb_80X/Optimization/" << Opti[ip] << "/";
    //Making dir to store output files
    boost::filesystem::create_directories(OutPath.str().c_str());

    ostringstream OutFile;
    ostringstream InHist;

    for(unsigned int is = 0; is < n_syst; is++){

      cout << "Running for " << SystFile[is] << endl;
      InHist.str("");  OutFile.str("");
      InHist.clear();  OutFile.clear();
      InHist << "h_mass_X_bin1" << SystHist[is];
      OutFile << OutPath.str() << model << "_XDists_" << coupling << SystFile[is] << ".root";

      TFile *fout = new TFile(OutFile.str().c_str(), "RECREATE");

      for(unsigned int im = 0; im < n; im++){

	ostringstream InFile;
	InFile.str("");
	InFile.clear();
	InFile << InPath.str() << model << "/" << model << "ToGJ_M-" << Mass[im]   << "_" << coupling << ".root";

	TFile *f_sig = new TFile(InFile.str().c_str(), "READ");
	TH1F *h_sig = (TH1F*)f_sig->Get(InHist.str().c_str());

	h_sig->Rebin(10);

	ostringstream OutHist;
	OutHist << "h_" << model << "ToGJ_" << Era << "_M" << Mass[im] << "_XDists";

	h_sig->GetXaxis()->SetTitle("");
	h_sig->GetYaxis()->SetTitle("");
	h_sig->SetName(OutHist.str().c_str());

	fout->cd();
	h_sig->Write(OutHist.str().c_str(), TObject::kWriteDelete);
      }    
      fout->Close();
    }
  }
}

void QstarXdists_CheckingInterpoltionQuality(){

  TString model = "Qstar";
  TString coupling = "f-1p0";

  TString Era = "Spring16";

  const int n = 2;
  Double_t Mass[n] = {5000, 7000};//, 8000, 9000};

  TString InPath = "PA_Results/Spring16_36813pb_80X/ResShapes_n_Syst/";

  TString OutPath = "MassXDist/Spring16_36813pb_80X/QualityCheck/";
  //Making dir to store output files
  boost::filesystem::create_directories("MassXDist/Spring16_36813pb_80X/QualityCheck/");

  ostringstream OutFile;
  ostringstream InHist;

  InHist.str("");  OutFile.str("");
  InHist.clear();  OutFile.clear();
  InHist << "h_mass_X_bin1";
  OutFile << OutPath << model << "_XDists_" << coupling << ".root";

  TFile *fout = new TFile(OutFile.str().c_str(), "RECREATE");

  for(unsigned int im = 0; im < n; im++){

    ostringstream InFile;
    InFile.str("");
    InFile.clear();
    InFile << InPath << model << "/" << model << "ToGJ_M-" << Mass[im]   << "_" << coupling << ".root";

    TFile *f_sig = new TFile(InFile.str().c_str(), "READ");
    TH1F *h_sig = (TH1F*)f_sig->Get(InHist.str().c_str());

    h_sig->Rebin(10);

    ostringstream OutHist;
    OutHist << "h_" << model << "ToGJ_" << Era << "_M" << Mass[im] << "_XDists";

    h_sig->GetXaxis()->SetTitle("");
    h_sig->GetYaxis()->SetTitle("");
    h_sig->SetName(OutHist.str().c_str());

    fout->cd();
    h_sig->Write(OutHist.str().c_str(), TObject::kWriteDelete);
  }
  fout->Close();
}




