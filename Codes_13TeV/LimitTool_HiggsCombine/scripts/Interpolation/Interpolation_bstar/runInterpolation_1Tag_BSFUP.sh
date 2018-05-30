#!/bin/bash

pwd=$PWD

#Only Change Input Dir, coupling, tag and mass points required (defined inside the runInterpolation() function).
#Put inputDir inside Bstar_MassXDist dir and structure inside inputDir should be like ${inputDir}/MC/Signal_Bstar/${coupling}/files here.

Model=Bstar ##Bstar   ##TO_CHANGE

tag=0bTag

InPath=${pwd}/MassXDist
InSigDir=${Model}f1p0_80X_ResShapes_latestSyst_full2016lumi/ ##TO_CHANGE

outputDir=${pwd}/InterpolatedFiles/80X/

if [ ! -d ${outputDir} ]; then
echo "--------- Making Directory ${outputDir} -----------"
mkdir -pv ${outputDir}
chmod 775 ${outputDir}
fi

#for coupling in f1p0 f0p5 f0p1
for coupling in f-1p0 
do

echo "++++++++++++"
echo "f = ${coupling}"
echo "++++++++++++"

###Name of output shape file: ResonanceShapes_ModelName_coupling_Tag(IfAny)_13TeV_EvtSelection_CMSSWversion_MCEra_SystUP/DOWN(IfAny).root
OutFileName=ResonanceShapes_${Model}_${tag}_${coupling}_13TeV_PhLID-JetTID-Pt190-nodeta-nodphi-CSVM_mass695_80X_Spring16_BSFUP  ##TO_CHANGE

InDir=${InPath}/${InSigDir}/

cat>Interpolation.C<<EOF
#define Interpolation_cxx
#include "Interpolation.h"
#include "TROOT.h"

int main(){
    gROOT->ProcessLine("#include <vector>");
    Interpolation m;
    m.runInterpolation();
    return 0; 
}

//Method to run interpolation
void Interpolation::runInterpolation(){

  TString OutFile = "${outputDir}/${OutFileName}.root";
  TString signalPath = "${InDir}";

  //output root file
  TFile *fnew = new TFile(OutFile,"RECREATE");
  fnew->cd();

  //Input X histoname
  TString Mass_X = "h_mass_X_bin1_BSF_${tag}_up"; //TO_CHANGE

  //Bin array for mass histogram
  const int nMassBins = 119;
  double MassBins[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};

  //Array of mass points over which interpolation is being done
  Int_t nMassPts = 10;
  
  //Int_t Signal_Mass[nMassPts] = {500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000};  //TO_CHANGE
  Int_t Signal_Mass[nMassPts] = {500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000}; 

  for(Int_t i = 0 ; i < nMassPts-1 ; ++i){

    Int_t Mass = 0;

    std::ostringstream mass1;
    mass1 << Signal_Mass[i];
    std::string mass1_str = mass1.str();
    std::ostringstream mass2;
    mass2 << Signal_Mass[i+1];
    std::string mass2_str = mass2.str();

    TString SignalFile1 = "${Model}ToGJ_M-"+mass1_str+"_${coupling}.root";
    TString SignalFile2 = "${Model}ToGJ_M-"+mass2_str+"_${coupling}.root";

    Mass = Signal_Mass[i];

    //loop from M1 to M2
    while(Mass < Signal_Mass[i+1]){
	
      std::vector<TH1D*> histOut;
      //histOut = getInterpolatedMass(signalPath+SignalFile1, Mass_X, signalPath+SignalFile2, Mass_X, Signal_Mass[i], Signal_Mass[i+1], Mass, nMassBins, MassBins);
      histOut = getInterpolatedMass_UnitBin(signalPath+SignalFile1, Mass_X, signalPath+SignalFile2, Mass_X, Signal_Mass[i], Signal_Mass[i+1], Mass);

      char OutHisto_ProbM[1000], OutHisto_ProbX[1000];

      sprintf(OutHisto_ProbM,"h_%s_%d","${Model}",Mass); //hist in the form of h_model_mass required for Higgs Combine 
      //sprintf(OutHisto_ProbM,"h_ProbM_Bstar%d%s",Mass,"_massVarBin");
      //sprintf(OutHisto_ProbX,"h_ProbX_Bstar%d%s",Mass,"_X");

      //mass histo normalized
      TH1D* finalHisto_ProbM = (TH1D*)histOut[0]->Clone(OutHisto_ProbM);
      //X distribution
      //TH1D* finalHisto_ProbX = (TH1D*)histOut[1]->Clone(OutHisto_ProbX);

      fnew->cd();
      //write histos to root file
      finalHisto_ProbM->Write();
      //finalHisto_ProbX->Write();

      //Mass Increment
      Mass += 50;

    } //while loop
  }//for loop

  fnew->Write();
  fnew->Close();

}


//Function to get Interpolated mass in variable binning
std::vector<TH1D*> Interpolation::getInterpolatedMass(TString filename1, TString histname1, TString filename2, TString histname2, int M1, int M2, int M, int nbins, double* bins){

  cout << "INTERPOLATED MASS = " << M << endl;
  cout << "GENERATED MASS USED = " << M1 << " - " << M2 << endl;
  if(M1 == M) cout << "INTERPOLATED MASS EXACTLY SAME AS GENERATED ONE" << endl;

  //Output Histogram vector	     
  std::vector<TH1D*> HistVector;
  HistVector.clear();

  //Get the X = Mgj/Mb* files
  TFile *file1 = new TFile(filename1, "READ"); 
  TFile *file2 = new TFile(filename2, "READ");

  //get the distributions 
  TH1D *hist1 = (TH1D*)file1->Get(histname1);
  TH1D *hist2 = (TH1D*)file2->Get(histname2);

  //output histogram for X distribution
  int Nbins = 14000; double lowbin = 0.0; double hibin = 5.0;
  TH1D *histout_X = new TH1D("histout_X", "X distribution for Mass M", Nbins, lowbin, hibin); //same as input X distribution

  //Output M distribution - one used in limits
  TH1D *histout_M = new TH1D("histout_M", "M distribution of Mass M", nbins, bins); //variable mass distribution

  double prob_X_M = 0;

  for(Int_t i = 1; i <= hist2->GetXaxis()->GetNbins(); ++i){
    
    prob_X_M = 0;

    //First get the probablity 
    hist1->Scale(1./hist1->Integral());
    hist2->Scale(1./hist2->Integral());

    double prob_X_M1 = hist1->GetBinContent(i);
    double prob_X_M2 = hist2->GetBinContent(i);

    prob_X_M = prob_X_M1 + ((prob_X_M2 - prob_X_M1)*((M-M1)/(M2-M1)));

    histout_X->SetBinContent(i,prob_X_M);

  }//for loop for output X distribution

  cout << "Integral of X = " << (double)histout_X->Integral() << endl;

  //Rebinning the X distribution to get M distribution back by doing X*M 
  //This is the final histogram required for study
  for(Int_t i = 1; i <= histout_M->GetXaxis()->GetNbins(); ++i){                   

    double bin_Min = histout_M->GetXaxis()->GetBinLowEdge(i);
    double bin_Max = histout_M->GetXaxis()->GetBinUpEdge(i);

    double val = 0.0;
    double thisMass_low = 0.0;
    double thisMass_hi = 0.0;

    //Get all y values between min and max
    for(Int_t bx = 1;  bx <= histout_X->GetXaxis()->GetNbins(); ++bx){
       
      thisMass_low = (histout_X->GetXaxis()->GetBinLowEdge(bx)) * ((double)M);
      thisMass_hi =  (histout_X->GetXaxis()->GetBinUpEdge(bx)) * ((double)M);

      if(thisMass_hi <= bin_Max && thisMass_low >= bin_Min){
	val += (histout_X->GetBinContent(bx));
      }else{
        double thisMass_mean = (thisMass_low + thisMass_hi)/2.0;
	if(thisMass_mean <= bin_Max && thisMass_mean > bin_Min) val += (histout_X->GetBinContent(bx));
      }
    }

    histout_M->SetBinContent(i, val);

  }

  cout << "Integral of M = " << (double)histout_M->Integral() << endl;

  HistVector.push_back(histout_M);
  HistVector.push_back(histout_X);

  if(HistVector.size() < 2 || HistVector.size() > 2) cout << " WARNING: SOMETHING IS WRONG WITH OUTPUT" << endl;

  cout << "------------------------------------------------" << endl;   
  cout<<""<<endl;
  cout<<""<<endl;

  return HistVector;

}

//Function to get Interpolated mass in Unit Binning
std::vector<TH1D*> Interpolation::getInterpolatedMass_UnitBin(TString filename1, TString histname1, TString filename2, TString histname2, int M1, int M2, int M){

  cout << "INTERPOLATED MASS = " << M << endl;
  cout << "GENERATED MASS USED = " << M1 << " - " << M2 << endl;
  if(M1 == M) cout << "INTERPOLATED MASS EXACTLY SAME AS GENERATED ONE" << endl;

  //Output Histogram vector	     
  std::vector<TH1D*> HistVector;
  HistVector.clear();

  //Get the X = Mgj/Mq* files
  TFile *file1 = new TFile(filename1, "READ"); 
  TFile *file2 = new TFile(filename2, "READ");

  //get the distributions 
  TH1D *hist1 = (TH1D*)file1->Get(histname1);
  TH1D *hist2 = (TH1D*)file2->Get(histname2);

  //Histogram for X distribution
  int Nbins = 14000; double lowbin = 0.0; double hibin = 5.0;
  TH1D *histout_X = new TH1D("histout_X", "X distribution for Mass M", Nbins, lowbin, hibin); //same as input X distribution

  //Output M distribution in unit bins - one used in Higgs combine
  TH1D *histout_M = new TH1D("histout_M", "M distribution of Mass M", 14000, 0.0, 14000.0); //UnitBin mass distribution

  double prob_X_M = 0;

  for(Int_t i = 1; i <= hist2->GetXaxis()->GetNbins(); ++i){
    
    prob_X_M = 0;

    //First get the probablity 
    hist1->Scale(1./hist1->Integral());
    hist2->Scale(1./hist2->Integral());

    double prob_X_M1 = hist1->GetBinContent(i);
    double prob_X_M2 = hist2->GetBinContent(i);

    prob_X_M = prob_X_M1 + ((prob_X_M2 - prob_X_M1)*((M-M1)/(M2-M1)));

    histout_X->SetBinContent(i,prob_X_M);

  }//for loop for output X distribution

  cout << "Integral of X = " << (double)histout_X->Integral() << endl;

  //Rebinning the X distribution to get M distribution back by doing X*M 
  //This is the final histogram required for study
  for(Int_t i = 1; i <= histout_M->GetXaxis()->GetNbins(); ++i){                   

    double bin_Min = histout_M->GetXaxis()->GetBinLowEdge(i);
    double bin_Max = histout_M->GetXaxis()->GetBinUpEdge(i);

    double val = 0.0;
    double thisMass_low = 0.0;
    double thisMass_hi = 0.0;

    //Get all y values between min and max
    for(Int_t bx = 1;  bx <= histout_X->GetXaxis()->GetNbins(); ++bx){
       
      thisMass_low = (histout_X->GetXaxis()->GetBinLowEdge(bx)) * ((double)M);
      thisMass_hi =  (histout_X->GetXaxis()->GetBinUpEdge(bx)) * ((double)M);

      if(thisMass_hi <= bin_Max && thisMass_low >= bin_Min){
	val += (histout_X->GetBinContent(bx));
      }else{
        double thisMass_mean = (thisMass_low + thisMass_hi)/2.0;
	if(thisMass_mean <= bin_Max && thisMass_mean > bin_Min) val += (histout_X->GetBinContent(bx));
      }
    }

    histout_M->SetBinContent(i, val);

  }

  cout << "Integral of M = " << (double)histout_M->Integral() << endl;

  HistVector.push_back(histout_M);
  HistVector.push_back(histout_X);

  if(HistVector.size() < 2 || HistVector.size() > 2) cout << " WARNING: SOMETHING IS WRONG WITH OUTPUT" << endl;

  cout << "------------------------------------------------" << endl;   
  cout<<""<<endl;
  cout<<""<<endl;

  return HistVector;

}

EOF

cat>Interpolation.h<<EOF
#ifndef Interpolation_h
#define Interpolation_h

#include<TH1.h>
#include<TFile.h>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TVectorD.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#ifdef __MAKECINT__                           
#pragma link C++ class vector<TH1D*>+;
#pragma link C++ class vector<double>+;
#endif

using namespace std;
using namespace ROOT;
 
class Interpolation {
public :
   Interpolation();
   virtual ~Interpolation();

   //My functions
   void runInterpolation(); 
   std::vector<TH1D*> getInterpolatedMass(TString filename1, TString histname1, TString filename2, TString histname2, int M1, int M2, int M, int nbins, double* bins);
   std::vector<TH1D*> getInterpolatedMass_UnitBin(TString filename1, TString histname1, TString filename2, TString histname2, int M1, int M2, int M);

};
#endif
 
#ifdef Interpolation_cxx
Interpolation::Interpolation()
{            
}
 
 
Interpolation::~Interpolation()
{
}
#endif 
EOF


g++ -Wno-deprecated Interpolation.C -o runInterpolation_${OutFileName}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`

./runInterpolation_${OutFileName}.exe

echo "INTERPOLATION DONE for ${coupling}"

rm runInterpolation_${OutFileName}.exe

rm Interpolation.C 
rm Interpolation.h

done
