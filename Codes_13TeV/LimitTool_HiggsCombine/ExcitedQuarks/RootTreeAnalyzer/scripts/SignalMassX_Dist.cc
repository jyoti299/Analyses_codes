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

#include <boost/filesystem.hpp>

using namespace std;
using namespace ROOT;


///TRYING TO MAKE THE MASS_X DIST DIRECTLY FROM UNITBIN INVT MASS DISTRIBUTION. HOWEVER, THE DISTRIBUTION FROM THIS IS DIFFERENT FROM THE ONE THAT WE
//GET BY FILLING THE M(GJ)/Mq* IN POSTANALYZER. THE REASON IS: HERE IN THIS, THE BIN SIZE OF INVT MASS HISTOGARM IS 1 BIN (0, 14000) WHILE THAT OF 
//X DIST HISTOGRAM IS 5/14000 = 0.00035. WHICH IS VERY SMALL THAN 1. SO WHILE USING THIS SCRIPT HERE, SUPPOSE TWO EVENTS WITH MASSES 1200.0 AND 
//1200.5 WILL BE IN THE SAME BIN FOR INVT MASS HIST, BUT THIS WILL BE IN DIFFERENT BINS FOR X HISTOGRAM AND HERE IN THE SCRIPT, SINCE WE CAN EITHER 
//TAKE THE BINCENTER OR THE BINUPEDGE, SO BOTH 1200.0 AND 1200.5 WILL BE CONSIDERED SAME AND WILL GET SAME X VALUE AND WILL GO IN SAME HIST. SO THIS
//MAKE THE DISTRIBUTION LITTLE DIFFERENT. THIS CAN BE USED IN EMERGENCEY BUT BEST IS TO PRODUCE THE X DISTRIBUTION FROM RUNNING POSTANALYZER AS IT 
//GIVES THE ACTUAL SHAPES. 
void MassX(){

  const Int_t n = 10;
  Double_t Mass[n] = {500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000};

  TSting InPath = "../PostAna_Results/";
  TString InDir = "PA_Qstar_13TeV_76X_2015/";
  TStrng SignalDir = "Signal_Qstar/";
  TString coupling = "f1p0";

  TString OutPath = "Interpolation/MassXDist/";
  //Making dir to store output files
  boost::filesystem::create_directories("Interpolation/MassXDist/Signal_Qstar/f1p0/");

  for(Int_t i = 0; i < n; i++){
    TString filename = "QstarToGJ_M" + Mass[i] + "_" + coupling + "_1.root";

    TFile *fIn = new TFile(InPath+InDir+SignalDir+coupling+filename, "READ");
    TFile *fOut = new TFile(OutPath+SignalDir+coupling+filename, "RECREATE");

    TH1F *hIn = (TH1F*)fIn->Get("h_mass_bin1_noMassCut");
    TF1F *hOut = new TH1F("h_X", "h_X", 14000, 0.0, 5.0);

    for(Int_t i = 0; i < hIn->GetNBinsX(); i++){

      Double_t N = hIn->GetBinContent(i+1);
      Double_t mbin = hIn->GetXaxis()->GetBinCenter(i+1);
      Double_t x = mbin/1000.0;
      Int_t xbin = hOut->GetXaxis()->FindBin(x);
      hOut->SetBinContent(xbin, N);

    }
  }
}
