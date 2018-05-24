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
