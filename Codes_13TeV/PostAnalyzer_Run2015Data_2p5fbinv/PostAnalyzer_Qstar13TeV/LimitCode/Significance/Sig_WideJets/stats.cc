#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>

#include <TGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <BAT/BCLog.h>

#include "binneddata.hh"
#include "fit.hh"
#include "statistics.hh"

using namespace std;

//##################################################################################################################################
// User Section 1
//
// (Change the code outside the User Sections only if you know what you are doing)
//

////////////////////////////////////////////////////////////////////////////////
// magic numbers
////////////////////////////////////////////////////////////////////////////////

// use Markov chain Monte Carlo (MCMC) to marginalize nuisance parameters
const int useMCMC = 0;
// IMPORTANT: With useMCMC = 1, the systematic uncertanties are included in the limit calculation by default. Use the PAR_NUIS[] array below to control what uncertainties are included

// number of samples of nuisance parameters for Bayesian MC integration
const int NSAMPLES=0; // 10000 (larger value is better but it also slows down the code. 10000 is a reasonable compromise between the speed and precision)
// IMPORTANT: With useMCMC = 0, the systematic uncertanties are included in the limit calculation only when NSAMPLES is greater than 0. Use the PAR_NUIS[] array below to control what uncertainties are included

// number of pseudoexperiments (when greater than 0, expected limit with +/- 1 and 2 sigma bands is calculated)
int NPES=0; // 200 (the more pseudo-experiments, the better. However, 200 is a reasonable choice)

// calculate significance estimator Sig = sgn(S)*sqrt{-2ln[L(B)/L(S+B)]}
const int calcSig = 1; // needs to be set to 0 for limit calculation

// use 6-parameter background fit function
const bool use6ParFit = 0;

// set the factor that defines the upper bound for the signal xs used by the MCMC as xsUpperBoundFactor*stat-only_limit
const double xsUpperBoundFactor=3.0;

// alpha (1-alpha=confidence interval)
const double ALPHA=0.05;

// left side tail
const double LEFTSIDETAIL=0.0;

// output file name
const string OUTPUTFILE="stats.root";

// center-of-mass energy
const double sqrtS = 13000.;

// histogram binning
//q* binning

const int NBINS=77;
double BOUNDARIES[NBINS+1] = {/*560, 593, */626, 660, 695, 731, 768,  806, 846, 887, 929, 972, 1017, 
                              1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609,   
                              1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 
                              2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592,
                              3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 
                              5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461,  
                              7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612,  
                              10954};
/*
//for dijet
const int NBINS=53;
double BOUNDARIES[NBINS+1] = {       1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
                               1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546,
                               2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704,
                               3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253,
                               5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320,
                               7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072,
                               10430 };
*/

// parameters
double SIGMASS=0;
const int NPARS=16;
const int NBKGPARS=(use6ParFit ? 6 : 4);
const int POIINDEX=0; // which parameter is "of interest"

string PAR_NAMES[NPARS]   = { "xs",  "lumi",  "jes", "jer",        "p0",        "p1",        "p2",         "p3",        "p4",         "p5", "n0", "n1", "n2", "n3", "n4", "n5" };

/*
//for dijet
double PAR_GUESSES[NPARS] = { 1E-3,    2445,    1.0,   1.0, 4.35374e-03, 7.46697e+00, 6.37033e+00,  2.14709e-01,          0.,           0.,    0,    0,    0,    0,    0,    0 };
double PAR_MIN[NPARS]     = {    0,     0.0,    0.0,   0.0,        -1E4,       -9999,       -9999,        -9999,       -9999,        -9999, -100, -100, -100, -100, -100, -100 };
double PAR_MAX[NPARS]     = {  1E3,     5E3,    2.0,   2.0,         1E4,        9999,        9999,         9999,        9999,         9999,  100,  100,  100,  100,  100,  100 };
double PAR_ERR[NPARS]     = { 1E-3,  122.25,   0.02,  0.10,       1e-04,       1e-01,       1e-01,        1e-03,       1e-02,        1e-03,    1,    1,    1,    1,    1,    1 };
*/


//for q*
//double PAR_GUESSES[NPARS] = { 1E-3,    2502.816,    1.0,   1.0,     2.54274e-10,   -5.27705e+00,  1.30363e+01,   1.40433e+00,       0.,           0.,   0,    0,    0,    0,    0,    0 };
double PAR_GUESSES[NPARS] = { 1E-3,    2502.816,    1.0,   1.0,     2.56158e-10,   -5.26605e+00, 1.30324e+01, 1.40376e+00,       0.,           0.,    0,    0,    0,    0,    0,    0 };
double PAR_MIN[NPARS]     = {    0,     0.0,    0.0,   0.0,        1E-20,       -9999,       -9999,        -9999,       -9999,        -9999, -100, -100, -100, -100, -100, -100 };
double PAR_MAX[NPARS]     = {  1E3,     5E3,    2.0,   2.0,         1E-2,        9999,        9999,         9999,        9999,         9999,  100,  100,  100,  100,  100,  100 };
double PAR_ERR[NPARS]     = { 1E-3,  150.169,   0.0125,  0.10,       1e-6,       1e-01,       1e-02,        1e-02,       1e-02,        1e-03,    1,    1,    1,    1,    1,    1 };



int PAR_TYPE[NPARS]       = {    1,       2,      2,     2,           0,           0,           0,            0,           0,            0,    3,    3,    3,    3,    3,    3 }; // // 1,2 = signal (2 not used in the fit); 0,3 = background (3 not used in the fit)
int PAR_NUIS[NPARS]       = {    0,       1,      1,     1,           0,           0,           0,            0,           0,            0,    4,    4,    4,    4,    4,    4 }; // 0 = not varied, >=1 = nuisance parameters with different priors (1 = Lognormal, 2 = Gaussian, 3 = Gamma, >=4 = Uniform)

//int PAR_NUIS[NPARS]       = {    0,       1,      1,     1,          0,            0,           0,            0,           0,            0,    4,    4,    4,    4,    4,    4 }; // all (same as above)
//int PAR_NUIS[NPARS]       = {    0,       0,      0,     0,          0,            0,           0,            0,           0,            0,    0,    0,    0,    0,    0,    0 }; // none
//int PAR_NUIS[NPARS]       = {    0,       1,      0,     0,          0,            0,           0,            0,           0,            0,    0,    0,    0,    0,    0,    0 }; // lumi only
//int PAR_NUIS[NPARS]       = {    0,       0,      1,     0,          0,            0,           0,            0,           0,            0,    0,    0,    0,    0,    0,    0 }; // jes only
//int PAR_NUIS[NPARS]       = {    0,       0,      0,     1,          0,            0,           0,            0,           0,            0,    0,    0,    0,    0,    0,    0 }; // jer only
//int PAR_NUIS[NPARS]       = {    0,       1,      1,     1,          0,            0,           0,            0,           0,            0,    0,    0,    0,    0,    0,    0 }; // all except background
//int PAR_NUIS[NPARS]       = {    0,       0,      0,     0,          0,            0,           0,            0,           0,            0,    4,    4,    4,    4,    4,    4 }; // background only

//
// End of User Section 1
//##################################################################################################################################

// input files vector
vector<string> INPUTFILES;

// covariance matrix
double COV_MATRIX[NPARS][NPARS];
TMatrixDSym covMatrix = TMatrixDSym(NBKGPARS);
TVectorD eigenValues = TVectorD(NBKGPARS);
TMatrixD eigenVectors = TMatrixD(NBKGPARS,NBKGPARS);

// constrain S to be positive in the S+B fit
const bool posS = 0;

// use B-only fit with fixed but non-zero signal when calculating the covariance matrix used for background systematics
const bool BonlyFitForSyst = 1;

// shift in the counter used to extract the covariance matrix
int shift = 1;

TH1D* HISTCDF=0; // signal CDF

////////////////////////////////////////////////////////////////////////////////
// printout function
////////////////////////////////////////////////////////////////////////////////
void printSignalXS(const Fitter* fit)
{
  cout << "Fitted signal xs: " << fit->getParameter(POIINDEX) << " +/- " << fit->getParError(POIINDEX) << endl;
}

////////////////////////////////////////////////////////////////////////////////
// fit functions
////////////////////////////////////////////////////////////////////////////////
Double_t fitQCD4Par(Double_t *m, Double_t *p)
{
    double x=m[0]/sqrtS;
    double logx=log(x);
    return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*logx);
}

Double_t fitQCD6Par(Double_t *m, Double_t *p)
{
    double x=m[0]/sqrtS;
    double logx=log(x);
    return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*logx+p[4]*pow(logx,2)+p[5]*pow(logx,3));
}

TF1 fit4par("fit4par",fitQCD4Par,BOUNDARIES[0],BOUNDARIES[NBINS],4);
TF1 fit6par("fit6par",fitQCD6Par,BOUNDARIES[0],BOUNDARIES[NBINS],6);

////////////////////////////////////////////////////////////////////////////////
// function integral -- 4-parameter background fit function
////////////////////////////////////////////////////////////////////////////////
double INTEGRAL_4PAR(double *x0, double *xf, double *par)
{
  double xs=par[0];
  double lumi=par[1];
  double jes=par[2];
  double jer=par[3];
  double p0=par[4];
  double p1=par[5];
  double p2=par[6];
  double p3=par[7];
  double n[4] = {0.};
  n[0]=par[10];
  n[1]=par[11];
  n[2]=par[12];
  n[3]=par[13];

  if( COV_MATRIX[0+shift][0+shift]>0. && (n[0]!=0. || n[1]!=0. || n[2]!=0. || n[3]!=0.) )
  {
    double g[4] = {0.};
    for(int v=0; v<4; ++v)
    {
      for(int k=0; k<4; ++k) g[k]=n[v]*eigenValues(v)*eigenVectors[k][v];
      p0 += g[0];
      p1 += g[1];
      p2 += g[2];
      p3 += g[3];
    }
  }

  fit4par.SetParameter(0,p0);
  fit4par.SetParameter(1,p1);
  fit4par.SetParameter(2,p2);
  fit4par.SetParameter(3,p3);

  // uses Simpson's 3/8th rule to compute the background integral over a short interval
  //double h=(xf[0]-x0[0])/3.;
  double a=x0[0];
  double b=xf[0];
  //double f1=fit4par.Eval(a);
  //double f2=fit4par.Eval(a+h);
  //double f3=fit4par.Eval(b-h);
  //double f4=fit4par.Eval(b);

  //double bkg=0.375*h*(f1 + 3*(f2 + f3) + f4);
  double bkg=fit4par.Integral(a,b);
  if(bkg<0.) bkg=1e-7;

  if(xs==0.0) return bkg;

  double xprimef=jes*(jer*(xf[0]-SIGMASS)+SIGMASS);
  double xprime0=jes*(jer*(x0[0]-SIGMASS)+SIGMASS);
  int bin1=HISTCDF->GetXaxis()->FindBin(xprimef);
  int bin2=HISTCDF->GetXaxis()->FindBin(xprime0);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF->GetNbinsX()) bin1=HISTCDF->GetNbinsX();
  if(bin2<1) bin2=1;
  if(bin2>HISTCDF->GetNbinsX()) bin2=HISTCDF->GetNbinsX();
  double sig=xs*lumi*(HISTCDF->GetBinContent(bin1)-HISTCDF->GetBinContent(bin2));

  return bkg+sig;
}

////////////////////////////////////////////////////////////////////////////////
// function integral -- 6-parameter background fit function
////////////////////////////////////////////////////////////////////////////////
double INTEGRAL_6PAR(double *x0, double *xf, double *par)
{
  double xs=par[0];
  double lumi=par[1];
  double jes=par[2];
  double jer=par[3];
  double p0=par[4];
  double p1=par[5];
  double p2=par[6];
  double p3=par[7];
  double p4=par[8];
  double p5=par[9];
  double n[6] = {0.};
  n[0]=par[10];
  n[1]=par[11];
  n[2]=par[12];
  n[3]=par[13];
  n[4]=par[14];
  n[5]=par[15];

  if( COV_MATRIX[0+shift][0+shift]>0. && (n[0]!=0. || n[1]!=0. || n[2]!=0. || n[3]!=0. || n[4]!=0. || n[5]!=0.) )
  {
    double g[6] = {0.};
    for(int v=0; v<6; ++v)
    {
      for(int k=0; k<6; ++k) g[k]=n[v]*eigenValues(v)*eigenVectors[k][v];
      p0 += g[0];
      p1 += g[1];
      p2 += g[2];
      p3 += g[3];
      p4 += g[4];
      p5 += g[5];
    }
  }

  fit6par.SetParameter(0,p0);
  fit6par.SetParameter(1,p1);
  fit6par.SetParameter(2,p2);
  fit6par.SetParameter(3,p3);
  fit6par.SetParameter(4,p4);
  fit6par.SetParameter(5,p5);

  // uses Simpson's 3/8th rule to compute the background integral over a short interval
  //double h=(xf[0]-x0[0])/3.;
  double a=x0[0];
  double b=xf[0];
  //double f1=fit6par.Eval(a);
  //double f2=fit6par.Eval(a+h);
  //double f3=fit6par.Eval(b-h);
  //double f4=fit6par.Eval(b);

  //double bkg=0.375*h*(f1 + 3*(f2 + f3) + f4);
  double bkg=fit6par.Integral(a,b);
  if(bkg<0.) bkg=1e-7;

  if(xs==0.0) return bkg;

  double xprimef=jes*(jer*(xf[0]-SIGMASS)+SIGMASS);
  double xprime0=jes*(jer*(x0[0]-SIGMASS)+SIGMASS);
  int bin1=HISTCDF->GetXaxis()->FindBin(xprimef);
  int bin2=HISTCDF->GetXaxis()->FindBin(xprime0);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF->GetNbinsX()) bin1=HISTCDF->GetNbinsX();
  if(bin2<1) bin2=1;
  if(bin2>HISTCDF->GetNbinsX()) bin2=HISTCDF->GetNbinsX();
  double sig=xs*lumi*(HISTCDF->GetBinContent(bin1)-HISTCDF->GetBinContent(bin2));

  return bkg+sig;
}

////////////////////////////////////////////////////////////////////////////////
// main function integral
////////////////////////////////////////////////////////////////////////////////
double INTEGRAL(double *x0, double *xf, double *par)
{
  if(use6ParFit) return INTEGRAL_6PAR(x0, xf, par);
  else           return INTEGRAL_4PAR(x0, xf, par);
}

////////////////////////////////////////////////////////////////////////////////
// main function
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  if(argc<=2) {
    cout << "Usage: stats MASS FINAL_STATE" << endl;
    return 0;
  }

  SIGMASS = atof(argv[1]);
  string masspoint = argv[1];

  string final_state = "qq";
  if(argc>2) final_state = argv[2];
  if(argc>3) NPES = atoi(argv[3]);
  int jobID = 0;
  if(argc>4) jobID = atoi(argv[4]);

  //##################################################################################################################################
  // User Section 2
  //
  // (Change the code outside the User Sections only if you know what you are doing)
  //

  /*
  //for dijet
  // input data file
  INPUTFILES.push_back("Data_and_ResonanceShapes/histo_data_mjj_Run2015D_silverJSON_2445invpb_JEC_Summer15_25nsV6_clean.root");

  // data histogram name
  string datahistname = "h_dat";

  // input signal files with resonance shapes
  string filename = "Data_and_ResonanceShapes/ResonanceShapes_qq_13TeV_Phys14Spring15Mix.root";
  if(final_state=="qg") filename = "Data_and_ResonanceShapes/ResonanceShapes_qg_13TeV_PU30_Spring15.root";
  if(final_state=="gg") filename = "Data_and_ResonanceShapes/ResonanceShapes_gg_13TeV_PU30_Spring15.root";

  // signal histogram name
  string histname = "h_" + final_state + "_" + masspoint;
  */
    
  //q*
  // input data file  
  INPUTFILES.push_back("Info_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_WideJets.root");

  // data histogram name
  string datahistname = "h_Data_massBin1";

  // input signal files with resonance shapes                                                                                      
  string filename = "IP_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_WideJets_Qstarf1p0.root";

 // signal histogram name               
  string histname = "h_Prob_Qstar" + masspoint + "_massVarBin";
  


  //
  // End of User Section 2
  //##################################################################################################################################

  if(BonlyFitForSyst) shift = 0;

  if(!posS) PAR_MIN[0] = -PAR_MAX[0];

  if(!use6ParFit) { PAR_TYPE[8]=3; PAR_TYPE[9]=3; PAR_NUIS[14]=0; PAR_NUIS[15]=0; }

  // initialize the covariance matrix
  for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

  // enable more detailed printout from the BAT MCMC
  BCLog::SetLogLevel(BCLog::detail);

  HISTCDF=getSignalCDF(filename.c_str(), histname.c_str(), filename.c_str(), histname.c_str(), 1., 1., 1.);

  assert(HISTCDF && SIGMASS>0);

  // get the data
  TH1D* data=getData(INPUTFILES, datahistname.c_str(), NBINS, BOUNDARIES);

  // create the output file
  string outputfile = OUTPUTFILE.substr(0,OUTPUTFILE.find(".root")) + "_" + masspoint + "_" + final_state + (argc>4 ? "_" : "") + (argc>4 ? to_string(jobID) : "") + ".root";
  TFile* rootfile=new TFile(outputfile.c_str(), "RECREATE");  rootfile->cd();

  // setup an initial fitter just to get pseudo-data, not to perform any actual fits
  //Fitter* initfit = new Fitter(data, INTEGRAL);
  //initfit->setRandomSeed(3157);
  //TH1D* pseudodata = initfit->makePseudoDataFromMC("data_0");
  //pseudodata->SaveAs("pseudodata.root");

  // setup the limit values
  double observedLowerBound, observedUpperBound;
  vector<double> expectedLowerBounds;
  vector<double> expectedUpperBounds;

  cout << "********************** pe=0 (data) **********************" << endl;

  // setup the fitter with the input from the signal+background fit
  Fitter* fit_data = new Fitter(data, INTEGRAL); // replace 'data' with 'pseudodata' if interested in fitting to pseudo-data derived from QCD MC
  fit_data->setRandomSeed(31415+jobID*100);
  fit_data->setPOIIndex(POIINDEX);
  //fit_data->setPrintLevel(0);
  for(int i=0; i<NPARS; i++) fit_data->defineParameter(i, PAR_NAMES[i].c_str(), PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

  // perform a signal+background fit possibly followed by a background-only fit with a fixed but non-zero signal
  for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) fit_data->fixParameter(i);
  if(BonlyFitForSyst) { fit_data->doFit(); cout << "Fit status: " << fit_data->getFitStatus() << endl; if(fit_data->getFitStatus().find("CONVERGED")==string::npos) { fit_data->fixParameter(0); fit_data->setParameter(0, 0.0); } else { printSignalXS(fit_data); fit_data->fixParameter(0); } }
  //fit_data->fixParameter(0); fit_data->setParameter(0, 0.0); // for MC studies with expected limits, fixing the signal xs to 0 in the fit
  fit_data->doFit(&COV_MATRIX[0][0], NPARS);
  if(!BonlyFitForSyst) printSignalXS(fit_data);
  cout << "Data fit status: " << fit_data->getFitStatus() << endl;
  double nll_SpB_data = fit_data->evalNLL();
  //cout << "NLL(S+B) = " << nll_SpB_data << endl;
  double sign_data = fit_data->getParameter(0)/fabs(fit_data->getParameter(0));
  fit_data->fixParameter(0); // a parameter needs to be fixed before its value can be changed
  fit_data->setParameter(0, 0.0); // set the xs value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
  if(calcSig) fit_data->doFit();
  double nll_B_data = fit_data->evalNLL();
  //cout << "NLL(B) = " << nll_B_data << endl;
  fit_data->setPrintLevel(0);
  if(jobID==0) fit_data->calcPull("pull_bkg_0")->Write();
  if(jobID==0) fit_data->calcDiff("diff_bkg_0")->Write();
  if(jobID==0) fit_data->write("fit_bkg_0");

  // Significance estimator: Sig = sgn(S)*sqrt{-2ln[L(B)/L(S+B)]}
  double nll_Diff_data = nll_B_data-nll_SpB_data;
  if(calcSig) cout << "Significance(data) = " << ( nll_Diff_data>0 ? sign_data*sqrt(2*nll_Diff_data) : 0. ) << endl;

  // calculate eigenvalues and eigenvectors
  for(int i = 0; i<NBKGPARS; ++i) { for(int j = 0; j<NBKGPARS; ++j) { covMatrix(i,j)=COV_MATRIX[i+shift][j+shift]; } }
  //covMatrix.Print();
  const TMatrixDSymEigen eigen_data(covMatrix);
  eigenValues = eigen_data.GetEigenValues();
  eigenValues.Sqrt();
  //eigenValues.Print();
  eigenVectors = eigen_data.GetEigenVectors();
  //eigenVectors.Print();

  fit_data->setParLimits(0, 0.0, PAR_MAX[0]); // for the posterior calculation, the signal xs has to be positive
  TGraph* post_data = 0;
  if(useMCMC==0)
  {
    post_data=fit_data->calculatePosterior(NSAMPLES);
    if(jobID==0) post_data->Write("post_0");
    cout << "Call limit reached: " << (fit_data->callLimitReached() ? "True" : "False") << endl;
  }
  else
  {
    post_data=fit_data->calculatePosterior(0);
    pair<double, double> statonly_bounds=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
    fit_data->setParLimits(0, 0.0, xsUpperBoundFactor*(statonly_bounds.second));
    post_data=fit_data->calculatePosterior((useMCMC ? 1 : NSAMPLES), useMCMC);
    //fit_data->PrintAllMarginalized("plots.ps");
    //fit_data->PrintResults("results.txt");
    if(jobID==0) post_data->Write("post_0");
  }

  // evaluate the limit
  pair<double, double> bounds_data=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
  observedLowerBound=bounds_data.first;
  observedUpperBound=bounds_data.second;

  // reset the covariance matrix
  for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

  // perform the PEs
  for(int pe=(jobID*NPES+1); pe<=(jobID*NPES+NPES); ++pe) {

    cout << "********************** pe=" << pe << " **********************" << endl;
    ostringstream pestr;
    pestr << "_" << pe;

    fit_data->fixParameter(0); // a parameter needs to be fixed before its value can be changed
    // setup the fitter with the input from the signal+background fit
    fit_data->setParameter(0, 0.0); // set the xs value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
    TH1D* hist = fit_data->makePseudoData((string("data")+pestr.str()).c_str()); // makes pseudo-data from the background fit
    fit_data->setParameter(0, PAR_GUESSES[0]);

    Fitter* fit = new Fitter(hist, INTEGRAL);
    fit->setPOIIndex(POIINDEX);
    fit->setPrintLevel(0);
    for(int i=0; i<NPARS; i++) fit->defineParameter(i, PAR_NAMES[i].c_str(), PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

    // perform a signal+background fit possibly followed by a background-only fit with a fixed but non-zero signal
    for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) fit->fixParameter(i);
    if(BonlyFitForSyst) { fit->doFit(); if(fit->getFitStatus().find("CONVERGED")==string::npos) { fit->fixParameter(0); fit->setParameter(0, 0.0); } else fit->fixParameter(0); }
    fit->doFit(&COV_MATRIX[0][0], NPARS);
    if(fit->getFitStatus().find("CONVERGED")==string::npos) continue; // skip the PE if the fit did not converge
    double nll_SpB = fit->evalNLL();
    //cout << "NLL(S+B) = " << nll_SpB << endl;
    double sign = fit->getParameter(0)/fabs(fit->getParameter(0));
    fit->fixParameter(0); // a parameter needs to be fixed before its value can be changed
    fit->setParameter(0, 0.0); // set the xs value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
    if(calcSig) fit->doFit();
    double nll_B = fit->evalNLL();
    //cout << "NLL(B) = " << nll_B << endl;
    fit->calcPull((string("pull_bkg")+pestr.str()).c_str())->Write();
    fit->calcDiff((string("diff_bkg")+pestr.str()).c_str())->Write();
    fit->write((string("fit_bkg")+pestr.str()).c_str());

    // Significance estimator: Sig = sgn(S)*sqrt{-2ln[L(B)/L(S+B)]}
    double nll_Diff = nll_B-nll_SpB;
    //if(nll_Diff<0.) cout << "-----> Negative NLL difference: " << nll_Diff << endl;
    if(calcSig) cout << "Significance(" << pe << ") = " << ( nll_Diff>0 ? sign*sqrt(2*nll_Diff) : 0. ) << endl;

    // calculate eigenvalues and eigenvectors
    for(int i = 0; i<NBKGPARS; ++i) { for(int j = 0; j<NBKGPARS; ++j) { covMatrix(i,j)=COV_MATRIX[i+shift][j+shift]; } }
    const TMatrixDSymEigen eigen(covMatrix);
    eigenValues = eigen.GetEigenValues();
    bool hasNegativeElement = false;
    for(int i = 0; i<NBKGPARS; ++i) { if(eigenValues(i)<0.) hasNegativeElement = true; }
    if(hasNegativeElement) continue; // this is principle should never happen. However, if it does, skip the PE
    eigenValues.Sqrt();
    eigenVectors = eigen.GetEigenVectors();

    fit->setParLimits(0, 0.0, PAR_MAX[0]); // for the posterior calculation, the signal xs has to be positive
    TGraph* post = 0;
    if(useMCMC==0)
    {
      post=fit->calculatePosterior(NSAMPLES);
      post->Write((string("post")+pestr.str()).c_str());
      cout << "Call limit reached in pe=" << pe << ": " << (fit->callLimitReached() ? "True" : "False") << endl;
    }
    else
    {
      post=fit->calculatePosterior(0);
      pair<double, double> statonly_bounds=evaluateInterval(post, ALPHA, LEFTSIDETAIL);
      fit->setParLimits(0, 0.0, xsUpperBoundFactor*(statonly_bounds.second));
      post=fit->calculatePosterior((useMCMC ? 1 : NSAMPLES), useMCMC);
      post->Write((string("post")+pestr.str()).c_str());
    }

    // evaluate the limit
    pair<double, double> bounds=evaluateInterval(post, ALPHA, LEFTSIDETAIL);
    if(bounds.first==0. && bounds.second>0.)
    {
      expectedLowerBounds.push_back(bounds.first);
      expectedUpperBounds.push_back(bounds.second);
    }

    // reset the covariance matrix
    for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

    delete fit;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // print the results
  ////////////////////////////////////////////////////////////////////////////////

  cout << "**********************************************************************" << endl;
  for(unsigned int i=0; i<expectedLowerBounds.size(); i++)
    cout << "expected bound(" << (jobID*NPES+i+1) << ") = [ " << expectedLowerBounds[i] << " , " << expectedUpperBounds[i] << " ]" << endl;

  cout << "\nobserved bound = [ " << observedLowerBound << " , " << observedUpperBound << " ]" << endl;

  if(LEFTSIDETAIL>0.0 && NPES>0) {
    cout << "\n***** expected lower bounds *****" << endl;
    double median;
    pair<double, double> onesigma;
    pair<double, double> twosigma;
    getQuantiles(expectedLowerBounds, median, onesigma, twosigma);
    cout << "median: " << median << endl;
    cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
    cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;
  }

  if(LEFTSIDETAIL<1.0 && NPES>0) {
    cout << "\n***** expected upper bounds *****" << endl;
    double median;
    pair<double, double> onesigma;
    pair<double, double> twosigma;
    getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
    cout << "median: " << median << endl;
    cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
    cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;
  }

  // close the output file
  rootfile->Close();

  return 0;
}
