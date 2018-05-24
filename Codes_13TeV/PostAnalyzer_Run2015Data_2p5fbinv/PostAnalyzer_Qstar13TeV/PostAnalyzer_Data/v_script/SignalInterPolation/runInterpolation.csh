#!/bin/tcsh

###################################################
# To interpolate signal shapes for intermediate mass points
# q* -> q gamma study
###################################################

setenv pwd $PWD
setenv inputDir Pt190jetEta2p4dphi2p5dEta2p0M560TID_2p2fb
setenv outputDir $PWD/InterpolatedFiles/${inputDir}/

if( ! -d ${outputDir}) then
echo "--------- Making Directory ${outputDir} -----------"
mkdir -pv ${outputDir}
chmod 775 ${outputDir}
endif


### Not to be touched below this in most cases -----------------------------------------------------------
##--------------------------------------------------------------------------------------------------------
foreach coupling (f1p0  f0p5  f0p1)

#############################################################
cat>QstarInterpolation.C<<EOF
#define QstarInterpolation_cxx
#include "QstarInterpolation.h"
#include "TROOT.h"


int main(){

    gROOT->ProcessLine("#include <vector>");
    QstarInterpolation m;
    m.runInterpolation();
    return 0; 
}


//Method to run interpolation
void QstarInterpolation::runInterpolation()
{


       TString OutFile = "${outputDir}/IP_${inputDir}_Qstar${coupling}.root";
       TString signalPath = "/eos/uscms/store/user/varun/work/13TeV/QstarGJ/PostAnalyzer/${inputDir}/Interpolation_Qstar${coupling}/";
       
       const int nSigFiles = 7;
       TString SignalFiles[nSigFiles] = {
       "QstarToGJ_M500_${coupling}_1.root",
       "QstarToGJ_M1000_${coupling}_1.root",
       "QstarToGJ_M2000_${coupling}_1.root",
       "QstarToGJ_M3000_${coupling}_1.root",
       "QstarToGJ_M4000_${coupling}_1.root",
       "QstarToGJ_M5000_${coupling}_1.root",
       "QstarToGJ_M6000_${coupling}_1.root"
//       "QstarToGJ_M7000_${coupling}_1.root",
//       "QstarToGJ_M8000_${coupling}_1.root",
//       "QstarToGJ_M9000_${coupling}_1.root"
       };

       int SignalMass[nSigFiles] = {
       500,
       1000,
       2000,
       3000,
       4000,
       5000,
       6000
//       7000,
//       8000,
//       9000
       };


    //Input X histoname
    TString Mass_X   = "h_mass_X_bin1";

    //bins we want for mass
   const int nMassBins = 119;
   double MassBins[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};


    //output root file with all histo
    TFile* f1=new TFile(OutFile,"RECREATE");
    f1->cd();

    //loop over all mass files
    for(Int_t i = 0 ; i < nSigFiles-1 ; ++i ){

	int Mass;
	Mass=0;

	//write the output histo of X and M with correct name
	Mass = SignalMass[i];

	//loop from M1 to M2
	while(Mass < SignalMass[i+1])
	{

	    char OutHistoM[1000], OutHistoProbM[1000], OutHistoXM[1000];

	    sprintf(OutHistoM,"h_Qstar%d%s",Mass,"_massVarBin");
	    sprintf(OutHistoProbM,"h_Prob_Qstar%d%s",Mass,"_massVarBin");
	    sprintf(OutHistoXM,"h_X_Qstar%d%s",Mass,"_X");

	    std::vector<TH1D*> histOut;
	    histOut = getInterpolatedMass(signalPath+SignalFiles[i], Mass_X, signalPath+SignalFiles[i+1], Mass_X, SignalMass[i], SignalMass[i+1], Mass, nMassBins, MassBins);

	    //mass histo normalized
	    TH1D* finalHisto_ProbM = (TH1D*)histOut[0]->Clone(OutHistoProbM);
	    //X distribution
	    TH1D* finalHisto_ProbX = (TH1D*)histOut[1]->Clone(OutHistoXM);

	    f1->cd();
	    //write these to root file
	    finalHisto_ProbM->Write();
	    finalHisto_ProbX->Write();

	    //write the output histo of X and M with correct name
	    if( Mass < 1000)  Mass += 100;
	    if( Mass >= 1000) Mass += 200;

	}//while loop


    }//for loop

    f1->Write();
    f1->Close();

}


//funcitn to get interpolated mass
std::vector<TH1D*> QstarInterpolation::getInterpolatedMass(TString filename1, TString histname1, 
	TString filename2, TString histname2, 
	int M1, int M2, int M,
	int nbins, double* bins)
{

    cout<<"Interpolated Mass will be = "<<M<<endl;
    cout<<"Generate Mass used for this = "<<M1<<" - "<<M2<<endl;
    if(M1 == M)cout<<"Interpolated mass SHOULD be EXACTLY same as generated one"<<endl;

    //Output Vector of Histo	     
    std::vector<TH1D*> HistVector;
    HistVector.clear();

    //Get the X= Mgj/Mq* files
    TFile* file1 = new TFile(filename1,"READ"); 
    TFile* file2 = new TFile(filename2,"READ");

    //get the distributions 
    TH1D* hist1=(TH1D*)file1->Get(histname1);
    TH1D* hist2=(TH1D*)file2->Get(histname2);


    //output hist for X
    int Nbins = 14000; 
    double lowbin=0.;
    double hibin=5.0;

    double prob_X_M=0;

    TH1D* histout_X = new TH1D("histout_X","X distribution for Mass M",Nbins,lowbin,hibin); //same as input X distributina nd same for all mass points
    TH1D* histout_M = new TH1D("histout_M","M distribution of Mass M",nbins,bins); //variable mass bins array

    for(Int_t i=1 ; i<= hist2->GetXaxis()->GetNbins(); ++i)
    {
	prob_X_M=0;
	//First get the probablity 
	hist1->Scale(1./hist1->Integral());
	hist2->Scale(1./hist2->Integral());

	double prob_X_M1 = hist1->GetBinContent(i);
	double prob_X_M2 = hist2->GetBinContent(i);

	prob_X_M = prob_X_M1 + ((prob_X_M2 - prob_X_M1)*((M-M1)/(M2-M1)));
	histout_X->SetBinContent(i,prob_X_M);

    }//loop over X distribution bins

    cout<<"Integral of X = "<<(double)histout_X->Integral()<<endl;

    //rebin the X distribution so that when we do X*M we get minimum of 1 GeV 
    //This will help to rebin the histogram as our binning for final mass histo
    for(int i=1 ; i<= histout_M->GetXaxis()->GetNbins(); ++i)                   
    {  

	double bin_Min = histout_M->GetXaxis()->GetBinLowEdge(i);
	double bin_Max = histout_M->GetXaxis()->GetBinUpEdge(i);

	double  val =0.;
	double thisMass_low=0.;
	double thisMass_hi =0.;

	//get all those y value between min and max
	for(int bx=1;  bx<=histout_X->GetXaxis()->GetNbins(); ++bx)
	{    
	    thisMass_low =  (histout_X->GetXaxis()->GetBinLowEdge(bx))* ((double)M);
	    thisMass_hi =  (histout_X->GetXaxis()->GetBinUpEdge(bx))* ((double)M);

	    if(thisMass_hi <= bin_Max && thisMass_low >=bin_Min)
	    {
		val +=(histout_X->GetBinContent(bx));
		//if(i>20 && i<23)cout<<" mass bin = "<<i<<" , bx ="<<bx<<" binMin-binMax "<<bin_Min<<" - "<<bin_Max<<"  mass:"<<thisMass_low<<"-"<<thisMass_hi<<endl;
	    }else{  double thisMass_mean= (thisMass_low+thisMass_hi)/2.0;
		if(thisMass_mean <= bin_Max && thisMass_mean >bin_Min)val +=(histout_X->GetBinContent(bx));

	    }
	}

	histout_M->SetBinContent(i, val);

    }


    cout<<"Integral of M = "<<(double)histout_M->Integral()<<endl;

    HistVector.push_back(histout_M);
    HistVector.push_back(histout_X);


    if(HistVector.size()<2 || HistVector.size()>2)std::cout<<" WARNING: Something wrong with output"<<std::endl;

    cout<<"---------------------------------------------"<<endl;   
    cout<<""<<endl;
    cout<<""<<endl;

    return HistVector;

}
EOF
######################################################

cat>QstarInterpolation.h<<EOF
#ifndef QstarInterpolation_h
#define QstarInterpolation_h

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



#ifdef __MAKECINT__                           
#pragma link C++ class vector<TH1D*>+;
#pragma link C++ class vector<double>+;
#endif

using namespace std;
using namespace ROOT;

 
class QstarInterpolation {
public :

   QstarInterpolation();

   virtual ~QstarInterpolation();

//My functions
void runInterpolation(); 

//std::vector<TH1D*> getInterpolatedMass(const char* filename1, const char* histname1, 
//                                	const char* filename2, const char* histname2, 
//	                                   int M1, int M2, int M, int nbins, double* bins);

std::vector<TH1D*> getInterpolatedMass(TString filename1, TString histname1, TString filename2, TString histname2, int M1, int M2, int M, int nbins, double* bins);

};
#endif
 
#ifdef QstarInterpolation_cxx
QstarInterpolation::QstarInterpolation()
{            
}
 
 
QstarInterpolation::~QstarInterpolation()
{
}
#endif 
EOF


g++ -Wno-deprecated QstarInterpolation.C -o runInterpolation_${coupling}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`
./runInterpolation_${coupling}.exe
echo "done for interpolation "
rm runInterpolation_${coupling}.exe

rm QstarInterpolation.C 
rm QstarInterpolation.h

end
