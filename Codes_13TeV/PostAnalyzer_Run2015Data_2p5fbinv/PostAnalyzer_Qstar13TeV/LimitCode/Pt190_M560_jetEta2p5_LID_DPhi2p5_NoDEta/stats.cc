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
const int useMCMC = 1;
// IMPORTANT: With useMCMC = 1, the systematic uncertanties are included in the limit calculation by default. Use the PAR_NUIS[] array below to control what uncertainties are included

// number of samples of nuisance parameters for Bayesian MC integration
const int NSAMPLES=5000; // 10000 (larger value is better but it also slows down the code. 10000 is a reasonable compromise between the speed and precision)
// IMPORTANT: With useMCMC = 0, the systematic uncertanties are included in the limit calculation only when NSAMPLES is greater than 0. Use the PAR_NUIS[] array below to control what uncertainties are included

// number of pseudoexperiments (when greater than 0, expected limit with +/- 1 and 2 sigma bands is calculated)
const int NPES=200; // 200 (the more pseudo-experiments, the better. However, 200 is a reasonable choice)

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
const int NBINS=75;
double BOUNDARIES[NBINS+1] = {/*560, 593, 626, 660, */695, 731, 768, 806, 846, 887, 929, 972, 1017, 
                              1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609,   
                              1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 
                              2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592,
                              3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 
                              5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461,  
                              7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612,  
                              10954};


// parameters
double SIGMASS=0;
const int NPARS=12;
const int NBKGPARS=4;//(use6ParFit ? 6 : 4);
const int POIINDEX=0; // which parameter is "of interest"
string PAR_NAMES[NPARS]   = { "xs",  "lumi",  "jes", "jer",        "p0",        "p1",        "p2",         "p3",    "n0", "n1", "n2", "n3"};
double PAR_GUESSES[NPARS] = { 1E-3,    2502.816,    1.0,   1.0,   5.21917e-06,  -5.13468e+00,  1.13301e+01,  1.40498e+00,  0,    0,    0,    0 };
//double PAR_GUESSES[NPARS] = { 1E-3,    2197.327,    1.0,   1.0,   3.11523e-10,   -2.90812e+00,   1.27400e+01,   1.38250e+00,  0,    0,    0,    0 };
double PAR_MIN[NPARS]     = {    0,     0.0,    0.0,   0.0,       1.0E-20,       -9999,       -9999,        -9999,   -100, -100, -100, -100 };
double PAR_MAX[NPARS]     = {  1E3,     5E3,    2.0,   2.0,        1.0E-02,        9999,        9999,         9999,   100,  100,  100,  100 };
double PAR_ERR[NPARS]     = { 1E-3,     150.169,   0.0125,  0.10,       1e-11,       1e-01,       1e-02,        1e-03,  1,    1,    1,    1 };
int PAR_TYPE[NPARS]       = {    1,       2,      2,     2,           0,           0,           0,            0,    3,    3,    3,    3 }; // // 1,2 = signal (2 not used in the fit); 0,3 = background (3 not used in the fit)
int PAR_NUIS[NPARS]       = {    0,       1,      1,     1,           0,           0,           0,            0,    4,    4,    4,    4 }; // 0 = not varied, >=1 = nuisance parameters with different priors (1 = Lognormal, 2 = Gaussian, 3 = Gamma, >=4 = Uniform)

//int PAR_NUIS[NPARS]       = {    0,       1,      1,     1,          0,            0,           0,            0,     4,    4,    4,    4 }; // all (same as above)
//int PAR_NUIS[NPARS]       = {    0,       1,      0,     0,          0,            0,           0,            0,     0,    0,    0,    0 }; // lumi only
//int PAR_NUIS[NPARS]       = {    0,       0,      1,     0,          0,            0,           0,            0,     0,    0,    0,    0 }; // jes only
//int PAR_NUIS[NPARS]       = {    0,       0,      0,     1,          0,            0,           0,            0,     0,    0,    0,    0 }; // jer only
//int PAR_NUIS[NPARS]       = {    0,       1,      1,     1,          0,            0,           0,            0,     0,    0,    0,    0 }; // all except background
//int PAR_NUIS[NPARS]       = {    0,       0,      0,     0,          0,            0,           0,            0,     4,    4,    4,    4 }; // background only

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
const bool BonlyFitForSyst = 1; // disabled to speed things up since not including systematics in the limits

// shift in the counter used to extract the covariance matrix
int shift = 1;

TH1D* HISTCDF=0; // signal CDF

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
  double n[NBKGPARS] = {0.};
  n[0]=par[8];
  n[1]=par[9];
  n[2]=par[10];
  n[3]=par[11];

  if( COV_MATRIX[0+shift][0+shift]>0. && (n[0]!=0. || n[1]!=0. || n[2]!=0. || n[3]!=0.) )
  {
    double g[NBKGPARS] = {0.};
    for(int v=0; v<NBKGPARS; ++v)
    {
      for(int k=0; k<NBKGPARS; ++k) g[k]=n[v]*eigenValues(v)*eigenVectors[k][v];
      p0 += g[0];
      p1 += g[1];
      p2 += g[2];
      p3 += g[3];
    }
  }

  // uses Simpson's 3/8th rule to compute the background integral over a short interval
  // also use a power series expansion to determine the intermediate intervals since the pow() call is expensive

  double dx=(xf[0]-x0[0])/3./sqrtS;
  double x=x0[0]/sqrtS;
  double logx=log(x);

  double a=pow(1-x,p1)/pow(x,p2+p3*logx);
  double b=dx*a/x/(x-1)*(p2+p1*x-p2*x-2*p3*(x-1)*logx);
  double c=0.5*dx*dx*a*( (p1-1)*p1/(x-1)/(x-1) - 2*p1*(p2+2*p3*logx)/(x-1)/x + (p2+p2*p2-2*p3+2*p3*logx*(1+2*p2+2*p3*logx))/x/x );
  double d=0.166666667*dx*dx*dx*a*( (p1-2)*(p1-1)*p1/(x-1)/(x-1)/(x-1) - 3*(p1-1)*p1*(p2+2*p3*logx)/(x-1)/(x-1)/x - (1+p2+2*p3*logx)*(p2*(2+p2) - 6*p3 + 4*p3*logx*(1+p2*p3*logx))/x/x/x + 3*p1*(p2+p2*p2-2*p3+2*p3*logx*(1+2*p2+2*p3*logx))/(x-1)/x/x );

  double bkg=(xf[0]-x0[0])*p0*(a+0.375*(b+c+d)+0.375*(2*b+4*c+8*d)+0.125*(3*b+9*c+27*d));
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
  return INTEGRAL_4PAR(x0, xf, par);
}

////////////////////////////////////////////////////////////////////////////////
// main function
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  if(argc<=1) {
    cout << "Usage: stats MASS" << endl;
    return 0;
  }

  SIGMASS = atof(argv[1]);
  string masspoint = argv[1];

  //##################################################################################################################################
  // User Section 2
  //
  // (Change the code outside the User Sections only if you know what you are doing)
  //

  // input data file
  INPUTFILES.push_back("LimitInput_Data_Pt190_M560_jetEta2p5_LID_DPhi2p5_NoDEta.root");

  // data histogram name
  string datahistname = "hGJetInvtMass_data_fine";

  // input signal files with resonance shapes
  string filename = "ResonanceShapes_Pt190_M560_jetEta2p5_LID_DPhi2p5_NoDEta_Qstarf1p0.root";

  // signal histogram name
  string histname = "h_Prob_Qstar" + masspoint + "_massVarBin";

  //
  // End of User Section 2
  //##################################################################################################################################

  if(BonlyFitForSyst) shift = 0;

  if(!posS) PAR_MIN[0] = -PAR_MAX[0];

  // initialize the covariance matrix
  for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

  // enable more detailed printout from the BAT MCMC
  BCLog::SetLogLevel(BCLog::detail);

  HISTCDF=getSignalCDF(filename.c_str(), histname.c_str(), filename.c_str(), histname.c_str(), 1., 1., 1.);

  assert(HISTCDF && SIGMASS>0);

  // get the data
  TH1D* data=getData(INPUTFILES, datahistname.c_str(), NBINS, BOUNDARIES);

  // create the output file
  string outputfile = OUTPUTFILE.substr(0,OUTPUTFILE.find(".root")) + "_" + masspoint + ".root";
  TFile* rootfile=new TFile(outputfile.c_str(), "RECREATE");  rootfile->cd();

  // setup the limit values
  double observedLowerBound, observedUpperBound;
  vector<double> expectedLowerBounds;
  vector<double> expectedUpperBounds;

  cout << "********************** pe=0 (data) **********************" << endl;

  // setup the fitter with the input from the signal+background fit
  Fitter* fit_data = new Fitter(data, INTEGRAL);
  fit_data->setPOIIndex(POIINDEX);
  //fit_data->setPrintLevel(0);
  for(int i=0; i<NPARS; i++) fit_data->defineParameter(i, PAR_NAMES[i].c_str(), PAR_GUESSES[i], PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

  // perform a signal+background fit possibly followed by a background-only fit with a fixed but non-zero signal
  for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) fit_data->fixParameter(i);
  if(BonlyFitForSyst) { fit_data->doFit(); if(fit_data->getFitStatus().find("CONVERGED")==string::npos) { fit_data->fixParameter(0); fit_data->setParameter(0, 0.0); } else fit_data->fixParameter(0); }
  fit_data->doFit(&COV_MATRIX[0][0], NPARS);
  cout << "Data fit status: " << fit_data->getFitStatus() << endl;
  fit_data->fixParameter(0); // a parameter needs to be fixed before its value can be changed
  fit_data->setParameter(0, 0.0); // set the xs value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
  fit_data->setPrintLevel(0);
  fit_data->calcPull("pull_bkg_0")->Write();
  fit_data->calcDiff("diff_bkg_0")->Write();
  fit_data->write("fit_bkg_0");

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
    post_data->Write("post_0");
    cout << "Call limit reached: " << (fit_data->callLimitReached() ? "True" : "False") << endl;
  }
  else
  {
    post_data=fit_data->calculatePosterior(0);
    pair<double, double> statonly_bounds=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
    fit_data->setParLimits(0, 0.0, xsUpperBoundFactor*(statonly_bounds.second));
    post_data=fit_data->calculatePosterior(NSAMPLES, useMCMC);
    //fit_data->PrintAllMarginalized("plots.ps");
    //fit_data->PrintResults("results.txt");
    post_data->Write("post_0");
  }

  // evaluate the limit
  pair<double, double> bounds_data=evaluateInterval(post_data, ALPHA, LEFTSIDETAIL);
  observedLowerBound=bounds_data.first;
  observedUpperBound=bounds_data.second;

  // reset the covariance matrix
  for(int i = 0; i<NPARS; ++i) { for(int j = 0; j<NPARS; ++j) COV_MATRIX[i][j]=0.; }

  // perform the PEs
  for(int pe=1; pe<=NPES; ++pe) {

    cout << "********************** pe=" << pe << " **********************" << endl;
    ostringstream pestr;
    pestr << "_" << pe;

    fit_data->fixParameter(0); // a parameter needs to be fixed before its value can be changed
    // setup the fitter with the input from the signal+background fit
    fit_data->setParameter(0, 0.0); // set the xs value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
    TH1D* hist = fit_data->makePseudoData((string("data")+pestr.str()).c_str());
    fit_data->setParameter(0, PAR_GUESSES[0]);

    Fitter* fit = new Fitter(hist, INTEGRAL);
    fit->setPOIIndex(POIINDEX);
    fit->setPrintLevel(0);
    for(int i=0; i<NPARS; i++) fit->defineParameter(i, PAR_NAMES[i].c_str(), fit_data->getParameter(i), PAR_ERR[i], PAR_MIN[i], PAR_MAX[i], PAR_NUIS[i]);

    // perform a signal+background fit possibly followed by a background-only fit with a fixed but non-zero signal
    for(int i=0; i<NPARS; i++) if(PAR_TYPE[i]>=2 || PAR_MIN[i]==PAR_MAX[i]) fit->fixParameter(i);
    if(BonlyFitForSyst) { fit->doFit(); if(fit->getFitStatus().find("CONVERGED")==string::npos) { fit->fixParameter(0); fit->setParameter(0, 0.0); } else fit->fixParameter(0); }
    fit->doFit(&COV_MATRIX[0][0], NPARS);
    if(fit->getFitStatus().find("CONVERGED")==string::npos) continue; // skip the PE if the fit did not converge
    fit->fixParameter(0); // a parameter needs to be fixed before its value can be changed
    fit->setParameter(0, 0.0); // set the xs value to 0 to get the B component of the S+B fit (for calculating pulls and generating pseudo-data)
    fit->calcPull((string("pull_bkg")+pestr.str()).c_str())->Write();
    fit->calcDiff((string("diff_bkg")+pestr.str()).c_str())->Write();
    fit->write((string("fit_bkg")+pestr.str()).c_str());

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
      post=fit->calculatePosterior(NSAMPLES, useMCMC);
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
    cout << "expected bound(" << (i+1) << ") = [ " << expectedLowerBounds[i] << " , " << expectedUpperBounds[i] << " ]" << endl;

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
