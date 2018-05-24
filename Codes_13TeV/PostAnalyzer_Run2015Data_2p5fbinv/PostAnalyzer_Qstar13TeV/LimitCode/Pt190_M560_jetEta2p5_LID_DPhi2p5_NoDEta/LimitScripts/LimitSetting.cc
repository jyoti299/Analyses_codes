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

void LimitSettingFiles_Data(){

  TFile *f = new TFile("/eos/uscms/store/user/rocky86/13TeV/PostAnalyzer_Results/Pt190_M560_jetEta2p5_LID_DPhi2p5_NoDEta/Data/Run2015D_PromptReco.root");
  TFile *fnew = new TFile("../LimitInput_Data_Pt190_M560_jetEta2p5_LID_DPhi2p5_NoDEta.root", "RECREATE");

  TH1F *h = (TH1F*)f->Get("h_GJetInvtMass_UnitBin_MassCut");
  h->GetXaxis()->SetTitle(""); 
  h->GetYaxis()->SetTitle("");

  h->SetName("hGJetInvtMass_data_fine");
  h->Write("hGJetInvtMass_data_fine", TObject::kWriteDelete);

  fnew->Close();

}

void LimitSettingFiles_qstarSignal_f1p0(){

  string coupling = "f1p0";
  int mass[6] = {1000, 2000, 3000, 4000, 5000, 7000};

  ostringstream OutFile;
  OutFile << "../Resonance_Shapes_qstar_" << coupling << ".root";
  TFile *fnew = new TFile(OutFile.str().c_str(), "RECREATE");

  int nbins = 119;
  double xbins[120] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000 };
    
  for(int i = 0; i < 6; i++){

    ostringstream InputFile;
    InputFile << "/eos/uscms/store/user/lpcqstar/13TeV/PostAnalyzer_Results/Qstar13TeV/MC/Signal_Qstar/QstarToGJ_M" << mass[i] << "_" << coupling << ".root";
    TFile *f = new TFile(InputFile.str().c_str());
    fnew->cd();

    ostringstream hist;
    hist << "h_GJetInvtMass_VarBin_MassCut";

    TH1F *h = (TH1F*)f->Get(hist.str().c_str());

    ostringstream histname;
    histname << "h_GJet_" << mass[i];

    TH1F *hnew = new TH1F(histname.str().c_str(), histname.str().c_str(), nbins, xbins);
    hnew->GetYaxis()->SetTitle("Probability");
    hnew->GetXaxis()->SetTitle("M_#gamma+Jet");

    for(int i = 0; i < nbins; i++){

      double val = h->GetBinContent(i+1);
      hnew->SetBinContent(i+1, val);
  }

    hnew->Scale(1.0/hnew->Integral());

    double max = 0;
    for(int j = 0; j < nbins; j++){
      if(hnew->GetBinContent(j+1) > 0){
	if(hnew->GetBinContent(j+1) > max){
	  max = hnew->GetBinContent(j+1);
	}
      }
    }

    hnew->GetYaxis()->SetRangeUser(0, max+0.02);

    hnew->Write();
  }
  fnew->Close();

}

void LimitSetting_signal_ProbWithHistError_example(){

  string coupling = "f0p1";
  string tag = "1bTag";
  int mass[10] = {500, 700, 1000, 1200, 1500, 1700, 2000, 2500, 3000, 4000};

  ostringstream OutFile;
  OutFile << "Resonance_Shapes_GbJet_" << coupling << "_" << tag << ".root";
  TFile *fnew = new TFile(OutFile.str().c_str(), "RECREATE");

  for(int i = 0; i < 10; i++){
    ostringstream InputFile;
    InputFile << "BStarSignal/BstarToGJ_M" << mass[i] << "_" << coupling << ".root";
    TFile *f = new TFile(InputFile.str().c_str());
    fnew->cd();

    ostringstream hist;
    hist << "h_CSVL_GBJetInvtMass_VariableBinning_noErr_" << tag;
    TH1F *h = (TH1F*)f->Get(hist.str().c_str());

    h->Scale(1.0/h->Integral());      
    h->GetXaxis()->SetTitle("");
    h->GetYaxis()->SetTitle("Probability"); h->GetYaxis()->CenterTitle();
    h->SetLineColor(1);
    h->SetFillColor(kWhite);
    ostringstream histname;
    histname << "h_GbJet_" << mass[i];
    h->SetName(histname.str().c_str());
   
    h->Write(histname.str().c_str(), TObject::kWriteDelete);
  }
  fnew->Close();
}


double histError(double val)
{
  const double alpha = 1 - 0.6827;

  if(val<25. && val>0.) return (0.5*TMath::ChisquareQuantile(1-alpha/2,2*(val+1))-0.5*TMath::ChisquareQuantile(alpha/2,2*val))/2.0; 
  else if(val==0.) return 0.5*TMath::ChisquareQuantile(1-alpha/2,2*(val+1)); // special case of 0 events with one-sided error bar
  else return sqrt(val); // for val>25 error bars with correct coverage are essentially symmetric
}


TH1D* getSignalCDF(const char* filename1, const char* histname1, const char* filename2, const char* histname2, const double BR, const double eff_h, const double eff_l, const std::string& postfix)
{
  TH1D* h_pdf=new TH1D("h_pdf", "Resonance Shape", 14000, 0, 14000);
  TH1D* h_cdf=new TH1D(("h_cdf"+postfix).c_str(), "Resonance Shape CDF", 14000, 0, 14000);

  TFile* histfile1=new TFile(filename1);
  TFile* histfile2=new TFile(filename2);

  TH1D* hist1=(TH1D*)histfile1->Get(histname1);
  TH1D* hist2=(TH1D*)histfile2->Get(histname2);

  hist1->Scale(eff_h*BR);
  hist2->Scale(eff_l*(1-BR));
  hist1->Add(hist2);

  for(int i=1; i<=hist1->GetNbinsX(); i++){

    int bin_min = h_pdf->GetXaxis()->FindBin(hist1->GetXaxis()->GetBinLowEdge(i)+0.5);
    int bin_max = h_pdf->GetXaxis()->FindBin(hist1->GetXaxis()->GetBinUpEdge(i)-0.5);
    double bin_content = hist1->GetBinContent(i)/double(bin_max-bin_min+1);
    if(bin_content<0.) bin_content=0.; // some extrapolated resonance shapes can have a small negative bin content in some of the bins. this protects against such cases
    for(int b=bin_min; b<=bin_max; b++){
       h_pdf->SetBinContent(b, bin_content);
    }
  }

  h_pdf->Scale(1./h_pdf->Integral());

  for(int i=1; i<=h_cdf->GetNbinsX(); i++){

    int bin_min = h_pdf->GetXaxis()->FindBin(h_cdf->GetXaxis()->GetBinLowEdge(i)+0.5);
    int bin_max = h_pdf->GetXaxis()->FindBin(h_cdf->GetXaxis()->GetBinUpEdge(i)-0.5);

    double curr = 0;
    for(int b=bin_min; b<=bin_max; b++){
       curr+=h_pdf->GetBinContent(b);
    }

    double prev=h_cdf->GetBinContent(i-1);

    h_cdf->SetBinContent(i, prev+curr);
  }

  histfile1->Close();
  histfile2->Close();
  delete h_pdf;

  return h_cdf;
}


TH1D* SignalCDF(){
  string filename_0bTag = "Resonance_Shapes_GbJet_f1p0_0bTag.root";
  string histname = "h_GbJet_1000";
  TH1D *Hist_cdf = getSignalCDF(filename_0bTag.c_str(), histname.c_str(), filename_0bTag.c_str(), histname.c_str(), 1, 1, 1, "_0Tag");

  return Hist_cdf;
}



double Integral(double x0, double xf){

  double xs = 1;
  double lumi = 19711.4;
  double jes= 1.0;
  double jer= 1.0;
  double eff= 0.406169;
  double norm= 1.02536e-07;
  double p1= -1.46933;
  double p2= 11.9965;
  double p3= 1.33157;

  double dx=(xf-x0)/3./8000.0;
  double x=x0/8000.0;
  double logx=log(x);

  double a=pow(1-x,p1)/pow(x,p2+p3*logx);
  double b=dx*a/x/(x-1)*(p2+p1*x-p2*x-2*p3*(x-1)*logx);
  double c=0.5*dx*dx*a*( (p1-1)*p1/(x-1)/(x-1) - 2*p1*(p2+2*p3*logx)/(x-1)/x + (p2+p2*p2-2*p3+2*p3*logx*(1+2*p2+2*p3*logx))/x/x );
  double d=0.166666667*dx*dx*dx*a*( (p1-2)*(p1-1)*p1/(x-1)/(x-1)/(x-1) - 3*(p1-1)*p1*(p2+2*p3*logx)/(x-1)/(x-1)/x - (1+p2+2*p3*logx)*(p2*(2+p2) - 6*p3 + 4*p3*logx*(1+p2*p3*logx))/x/x/x + 3*p1*(p2+p2*p2-2*p3+2*p3*logx*(1+2*p2+2*p3*logx))/(x-1)/x/x );

  double bkg=(xf-x0)*norm*(a+0.375*(b+c+d)+0.375*(2*b+4*c+8*d)+0.125*(3*b+9*c+27*d));
  if(bkg<0.) bkg=0.0000001;

  if(xs==0.0) return bkg;

  TH1D* HISTCDF_1Tag = SignalCDF();
  double SIGMASS = 1000;

  double xprimef=jes*(jer*(xf-SIGMASS)+SIGMASS);
  double xprime0=jes*(jer*(x0-SIGMASS)+SIGMASS);
  int bin1=HISTCDF_1Tag->GetXaxis()->FindBin(xprimef);
  int bin2=HISTCDF_1Tag->GetXaxis()->FindBin(xprime0);
  if(bin1<1) bin1=1;
  if(bin1>HISTCDF_1Tag->GetNbinsX()) bin1=HISTCDF_1Tag->GetNbinsX();
  if(bin2<1) bin1=1;
  if(bin2>HISTCDF_1Tag->GetNbinsX()) bin2=HISTCDF_1Tag->GetNbinsX();
  double sig=xs*eff*lumi*(HISTCDF_1Tag->GetBinContent(bin1)-HISTCDF_1Tag->GetBinContent(bin2));

  return bkg+sig;
}


TH1D* getData(){

  const int nbins=38;
  double bins[nbins+1] = {   532,   576,   622,   669,   718,   768,   819,   872,   926,   982,  1040,  1099,  1160,  1223,
				 1288,  1355,  1424,  1495,  1568,  1643,  1720,  1799,  1881,  1965,  2052,  2141,  2233,  2328,
			     2425,  2525,  2628,  2734,  2843,  2955,  3071,  3190,  3313,  3439,  3569};

  TH1D *hist = new TH1D("normalized_data", "normalized data", nbins, bins);

  double categoryoffset = 10000;

  std::vector<std::string> filenames;
  filenames.push_back("Histos_data_GbJetInvtMass_0bTag.root");
  //  filenames.push_back("Histos_data_GbJetInvtMass_1bTag.root");

  unsigned nFiles = filenames.size();

  double minvalue = bins[0];
  double maxvalue = bins[nbins-1];

  assert( bins[0] < categoryoffset && maxvalue <= ( nFiles*categoryoffset ) );

  TFile *files[nFiles];
  TH1D *hists[nFiles];

  // open input files and get input histograms
  for(unsigned i=0; i<nFiles; ++i)
  {
    files[i] = new TFile(filenames[i].c_str());
    hists[i] = (TH1D*)files[i]->Get("hGbJetInvtMass_data_fine");
  }

  // normalize the data in the histogram
  for(int i=1; i<=hist->GetNbinsX(); ++i) {
    if(hist->GetXaxis()->GetBinUpEdge(i)<=(minvalue+categoryoffset))
    {
      double offset = 0.;
      int bin_min = hists[0]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)-offset+0.5);
      int bin_max = hists[0]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-offset-0.5);
      double val=hists[0]->Integral(bin_min,bin_max);
      double err=histError(val);
      double width=hist->GetBinWidth(i);
      hist->SetBinContent(i, val/width);
      hist->SetBinError(i, err/width);
    }
    else if(hist->GetXaxis()->GetBinUpEdge(i)>(minvalue+categoryoffset) && hist->GetXaxis()->GetBinUpEdge(i)<=(minvalue+2*categoryoffset) && nFiles>1)
    {
      double offset = categoryoffset;
      int bin_min = hists[1]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)-offset+0.5);
      int bin_max = hists[1]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-offset-0.5);
      double val=hists[1]->Integral(bin_min,bin_max);
      double err=histError(val);
      double width=hist->GetBinWidth(i);
      hist->SetBinContent(i, val/width);
      hist->SetBinError(i, err/width);
    }
    else if(hist->GetXaxis()->GetBinUpEdge(i)>(minvalue+2*categoryoffset) && nFiles>2)
    {
      double offset = 2*categoryoffset;
      int bin_min = hists[2]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinLowEdge(i)-offset+0.5);
      int bin_max = hists[2]->GetXaxis()->FindBin(hist->GetXaxis()->GetBinUpEdge(i)-offset-0.5);
      double val=hists[2]->Integral(bin_min,bin_max);
      double err=histError(val);
      double width=hist->GetBinWidth(i);
      hist->SetBinContent(i, val/width);
      hist->SetBinError(i, err/width);
      }
    }
  for(unsigned i=0; i<nFiles; ++i)
  {
    files[i]->Close();
  }


  return hist;
}


void nll(){

  TH1D* data = getData();

  double f=0.0;
  for(int bin=1; bin<=data->GetNbinsX(); bin++) {
    double binwidth=data->GetBinWidth(bin);
    double N=data->GetBinContent(bin)*binwidth;

    double x0=data->GetBinLowEdge(bin);
    double xf=data->GetBinLowEdge(bin+1);
    double mu=Integral(x0, xf);

    if(N==0.0) f += mu;
    else f -= (N*TMath::Log(mu) - TMath::LnGamma(N+1) - mu);

  }
  cout << "f = " << f << endl;
}


void Eff_f1p0(){
    std::string files[6] = {"BStarSignal/BstarToGJ_M1500_f1p0.root", "BStarSignal/BstarToGJ_M1700_f1p0.root", "BStarSignal/BstarToGJ_M2000_f1p0.root", "BStarSignal/BstarToGJ_M2500_f1p0.root", "BStarSignal/BstarToGJ_M3000_f1p0.root", "BStarSignal/BstarToGJ_M4000_f1p0.root"};


  for(int i = 0; i < 6; i++){
    cout << "-------------------------------------------------------------" << endl;
    cout << "for the file " << files[i].c_str() << " :::::::::: " << endl;
    cout << "-------------------------------------------------------------" << endl;

    TFile *f = new TFile(files[i].c_str(), "READ");
    TH1F *h = (TH1F*)f->Get("h_CutFlowTableWithWeights");

    Double_t binValue_total = h->GetBinContent(1);
    std::string binName = h->GetXaxis()->GetBinLabel(16);
    Double_t binValue = h->GetBinContent(16);

    Double_t div = binValue/binValue_total;
    
    cout << " eff for binName = " << binName << " is " << div << endl;
  
  }
} 
    

void Eff_f0p5(){
    std::string files[6] = {"BStarSignal/BstarToGJ_M1500_f0p5.root", "BStarSignal/BstarToGJ_M1700_f0p5.root", "BStarSignal/BstarToGJ_M2000_f0p5.root", "BStarSignal/BstarToGJ_M2500_f0p5.root", "BStarSignal/BstarToGJ_M3000_f0p5.root", "BStarSignal/BstarToGJ_M4000_f0p5.root"};


  for(int i = 0; i < 6; i++){
    cout << "-------------------------------------------------------------" << endl;
    cout << "for the file " << files[i].c_str() << " :::::::::: " << endl;
    cout << "-------------------------------------------------------------" << endl;

    TFile *f = new TFile(files[i].c_str(), "READ");
    TH1F *h = (TH1F*)f->Get("h_CutFlowTableWithWeights");

    Double_t binValue_total = h->GetBinContent(1);
    std::string binName = h->GetXaxis()->GetBinLabel(16);
    Double_t binValue = h->GetBinContent(16);

    Double_t div = binValue/binValue_total;
    
    cout << " eff for binName = " << binName << " is " << div << endl;
  
  }
} 
    

void Eff_f0p1(){
    std::string files[6] = {"BStarSignal/BstarToGJ_M1500_f0p1.root", "BStarSignal/BstarToGJ_M1700_f0p1.root", "BStarSignal/BstarToGJ_M2000_f0p1.root", "BStarSignal/BstarToGJ_M2500_f0p1.root", "BStarSignal/BstarToGJ_M3000_f0p1.root", "BStarSignal/BstarToGJ_M4000_f0p1.root"};


  for(int i = 0; i < 6; i++){
    cout << "-------------------------------------------------------------" << endl;
    cout << "for the file " << files[i].c_str() << " :::::::::: " << endl;
    cout << "-------------------------------------------------------------" << endl;

    TFile *f = new TFile(files[i].c_str(), "READ");
    TH1F *h = (TH1F*)f->Get("h_CutFlowTableWithWeights");

    Double_t binValue_total = h->GetBinContent(1);
    std::string binName = h->GetXaxis()->GetBinLabel(16);
    Double_t binValue = h->GetBinContent(16);

    Double_t div = binValue/binValue_total;
    
    cout << " eff for binName = " << binName << " is " << div << endl;
  
  }
} 
    

void getQuantiles(std::vector<double>& limits, double &median_, std::pair<double, double>& onesigma_, std::pair<double, double>& twosigma_) {
  unsigned int nit=limits.size();
  if(nit==0) return;

  // sort the vector with limits
  std::sort(limits.begin(), limits.end());

  // median for the expected limit
  median_ = TMath::Median(nit, &limits[0]);

  // quantiles for the expected limit bands
  double prob[4]; // array with quantile boundaries
  prob[0] = 0.021;
  prob[1] = 0.159;
  prob[2] = 0.841;
  prob[3] = 0.979;

  // array for the results
  double quantiles[4];
  TMath::Quantiles(nit, 4, &limits[0], quantiles, prob); // evaluate quantiles
  onesigma_.first=quantiles[1];
  onesigma_.second=quantiles[2];
  twosigma_.first=quantiles[0];
  twosigma_.second=quantiles[3];

  return;
}



void vec_700_f1p0(){

  vector<double> expectedUpperBounds;

  //for 700
  expectedUpperBounds.push_back( 0.023412  );
  expectedUpperBounds.push_back( 0.021796  );
  expectedUpperBounds.push_back( 0.0266338 );
  expectedUpperBounds.push_back( 0.0337399 ); 
  expectedUpperBounds.push_back( 0.0205289 );
  expectedUpperBounds.push_back( 0.0262506 );
  expectedUpperBounds.push_back( 0.0280901 );
  expectedUpperBounds.push_back( 0.0237389 );
  expectedUpperBounds.push_back( 0.040767  );
  expectedUpperBounds.push_back( 0.0404194 );
  expectedUpperBounds.push_back( 0.0217809 );
  expectedUpperBounds.push_back( 0.0131376 );
  expectedUpperBounds.push_back( 0.00987642);
  expectedUpperBounds.push_back( 0.0186185 );
  expectedUpperBounds.push_back( 0.0238894 );
  expectedUpperBounds.push_back( 0.0313396 );
  expectedUpperBounds.push_back( 0.0350311 );
  expectedUpperBounds.push_back( 0.0191511 );
  expectedUpperBounds.push_back( 0.0225638 );
  expectedUpperBounds.push_back( 0.0120219 );
  expectedUpperBounds.push_back( 0.0353878 );
  expectedUpperBounds.push_back( 0.0144611 );
  expectedUpperBounds.push_back( 0.0173779 );
  expectedUpperBounds.push_back( 0.0209759 );
  expectedUpperBounds.push_back( 0.0246181 );
  expectedUpperBounds.push_back( 0.0386755 );
  expectedUpperBounds.push_back( 0.0231747 );
  expectedUpperBounds.push_back( 0.0124035 );
  expectedUpperBounds.push_back( 0.0151178 );
  expectedUpperBounds.push_back( 0.0207849 );
  expectedUpperBounds.push_back( 0.0157627 );
  expectedUpperBounds.push_back( 0.0251036 );
  expectedUpperBounds.push_back( 0.0233862 );
  expectedUpperBounds.push_back( 0.0187476 );
  expectedUpperBounds.push_back( 0.0136487 );
  expectedUpperBounds.push_back( 0.0260121 );
  expectedUpperBounds.push_back( 0.0145157 );
  expectedUpperBounds.push_back( 0.0279511 );
  expectedUpperBounds.push_back( 0.0303302 );
  expectedUpperBounds.push_back( 0.0142907 );
  expectedUpperBounds.push_back( 0.0269211 );
  expectedUpperBounds.push_back( 0.0233173 );
  expectedUpperBounds.push_back( 0.040729  );
  expectedUpperBounds.push_back( 0.0104686 );
  expectedUpperBounds.push_back( 0.0124778 );
  expectedUpperBounds.push_back( 0.0456042 );
  expectedUpperBounds.push_back( 0.0347336 );
  expectedUpperBounds.push_back( 0.0196015 );
  expectedUpperBounds.push_back( 0.0164539 );
  expectedUpperBounds.push_back( 0.0309728 );
  expectedUpperBounds.push_back( 0.0218864 );
  expectedUpperBounds.push_back( 0.0384213 );
  expectedUpperBounds.push_back( 0.0206651 );
  expectedUpperBounds.push_back( 0.0200246 );
  expectedUpperBounds.push_back( 0.0320686 );
  expectedUpperBounds.push_back( 0.0144755 );
  expectedUpperBounds.push_back( 0.0598274 );
  expectedUpperBounds.push_back( 0.0179839 );
  expectedUpperBounds.push_back( 0.0315774 );
  expectedUpperBounds.push_back( 0.0212732 );
  expectedUpperBounds.push_back( 0.016481  );
  expectedUpperBounds.push_back( 0.0152657 );
  expectedUpperBounds.push_back( 0.0166026 );
  expectedUpperBounds.push_back( 0.0263277 );
  expectedUpperBounds.push_back( 0.0139103 );
  expectedUpperBounds.push_back( 0.0215133 );
  expectedUpperBounds.push_back( 0.0252056 );
  expectedUpperBounds.push_back( 0.0247376 );
  expectedUpperBounds.push_back( 0.023102  );
  expectedUpperBounds.push_back( 0.0177394 );
  expectedUpperBounds.push_back( 0.0124719 );
  expectedUpperBounds.push_back( 0.0149319 );
  expectedUpperBounds.push_back( 0.0185052 );
  expectedUpperBounds.push_back( 0.0396905 );
  expectedUpperBounds.push_back( 0.0178023 );
  expectedUpperBounds.push_back( 0.0248403 );
  expectedUpperBounds.push_back( 0.021955  );
  expectedUpperBounds.push_back( 0.0283869 );
  expectedUpperBounds.push_back( 0.0229585 );
  expectedUpperBounds.push_back( 0.0186887 );
  expectedUpperBounds.push_back( 0.0129286 );
  expectedUpperBounds.push_back( 0.0226899 );
  expectedUpperBounds.push_back( 0.0167671 );
  expectedUpperBounds.push_back( 0.00996496);
  expectedUpperBounds.push_back( 0.032557  );
  expectedUpperBounds.push_back( 0.0206787 );
  expectedUpperBounds.push_back( 0.0353958 );
  expectedUpperBounds.push_back( 0.0162372 );
  expectedUpperBounds.push_back( 0.0327624 );
  expectedUpperBounds.push_back( 0.0125892 );
  expectedUpperBounds.push_back( 0.0167509 );
  expectedUpperBounds.push_back( 0.0270854 );
  expectedUpperBounds.push_back( 0.0113644 );
  expectedUpperBounds.push_back( 0.027818  );
  expectedUpperBounds.push_back( 0.0162949 );
  expectedUpperBounds.push_back( 0.0435977 );
  expectedUpperBounds.push_back( 0.0162266 );
  expectedUpperBounds.push_back( 0.0174848 );
  expectedUpperBounds.push_back( 0.0252666 );
  expectedUpperBounds.push_back( 0.0193136 );
  expectedUpperBounds.push_back( 0.0224101 );
  expectedUpperBounds.push_back( 0.0125578 );
  expectedUpperBounds.push_back( 0.00469286);
  expectedUpperBounds.push_back( 0.0390504 );
  expectedUpperBounds.push_back( 0.0327268 );
  expectedUpperBounds.push_back( 0.0295909 );
  expectedUpperBounds.push_back( 0.0172767 );
  expectedUpperBounds.push_back( 0.021008  );
  expectedUpperBounds.push_back( 0.0312994 );
  expectedUpperBounds.push_back( 0.0112383 );
  expectedUpperBounds.push_back( 0.0184999 );
  expectedUpperBounds.push_back( 0.0123253 );
  expectedUpperBounds.push_back( 0.0382533 );
  expectedUpperBounds.push_back( 0.0104235 );
  expectedUpperBounds.push_back( 0.0126215 );
  expectedUpperBounds.push_back( 0.0183315 );
  expectedUpperBounds.push_back( 0.0220425 );
  expectedUpperBounds.push_back( 0.0142935 );
  expectedUpperBounds.push_back( 0.030207  );
  expectedUpperBounds.push_back( 0.0146657 );
  expectedUpperBounds.push_back( 0.041469  );
  expectedUpperBounds.push_back( 0.0228197 );
  expectedUpperBounds.push_back( 0.0109923 );
  expectedUpperBounds.push_back( 0.0196594 );
  expectedUpperBounds.push_back( 0.0440536 );
  expectedUpperBounds.push_back( 0.01131   );
  expectedUpperBounds.push_back( 0.02332   );
  expectedUpperBounds.push_back( 0.00914728);
  expectedUpperBounds.push_back( 0.0219825 );
  expectedUpperBounds.push_back( 0.0267611 );
  expectedUpperBounds.push_back( 0.0358867 );
  expectedUpperBounds.push_back( 0.0148645 );
  expectedUpperBounds.push_back( 0.0141732 );
  expectedUpperBounds.push_back( 0.025472  );
  expectedUpperBounds.push_back( 0.0148444 );
  expectedUpperBounds.push_back( 0.0180475 );
  expectedUpperBounds.push_back( 0.0234396 );
  expectedUpperBounds.push_back( 0.0326619 );
  expectedUpperBounds.push_back( 0.033486  );
  expectedUpperBounds.push_back( 0.0278049 );
  expectedUpperBounds.push_back( 0.0140947 );
  expectedUpperBounds.push_back( 0.0315105 );
  expectedUpperBounds.push_back( 0.0177444 );
  expectedUpperBounds.push_back( 0.0182161 );
  expectedUpperBounds.push_back( 0.0154919 );
  expectedUpperBounds.push_back( 0.0139297 );
  expectedUpperBounds.push_back( 0.0311908 );
  expectedUpperBounds.push_back( 0.0253498 );
  expectedUpperBounds.push_back( 0.0346194 );
  expectedUpperBounds.push_back( 0.0338241 );
  expectedUpperBounds.push_back( 0.031374  );
  expectedUpperBounds.push_back( 0.0127639 );
  expectedUpperBounds.push_back( 0.0210772 );
  expectedUpperBounds.push_back( 0.00560978);
  expectedUpperBounds.push_back( 0.0370263 );
  expectedUpperBounds.push_back( 0.01094   );
  expectedUpperBounds.push_back( 0.019216  );
  expectedUpperBounds.push_back( 0.0148489 );
  expectedUpperBounds.push_back( 0.0435854 );
  expectedUpperBounds.push_back( 0.0231567 );
  expectedUpperBounds.push_back( 0.0179433 );
  expectedUpperBounds.push_back( 0.047454  );
  expectedUpperBounds.push_back( 0.0276876 );
  expectedUpperBounds.push_back( 0.033979  );
  expectedUpperBounds.push_back( 0.0158326 );
  expectedUpperBounds.push_back( 0.0247411 );
  expectedUpperBounds.push_back( 0.0329395 );
  expectedUpperBounds.push_back( 0.0234387 );
  expectedUpperBounds.push_back( 0.0183245 );
  expectedUpperBounds.push_back( 0.0203523 );
  expectedUpperBounds.push_back( 0.021777  );
  expectedUpperBounds.push_back( 0.024875  );
  expectedUpperBounds.push_back( 0.0117607 );
  expectedUpperBounds.push_back( 0.0508472 );
  expectedUpperBounds.push_back( 0.022899  );
  expectedUpperBounds.push_back( 0.0179808 );
  expectedUpperBounds.push_back( 0.0062458 );
  expectedUpperBounds.push_back( 0.0106948 );
  expectedUpperBounds.push_back( 0.0204405 );
  expectedUpperBounds.push_back( 0.0132593 );
  expectedUpperBounds.push_back( 0.0129525 );
  expectedUpperBounds.push_back( 0.0133798 );
  expectedUpperBounds.push_back( 0.0149584 );
  expectedUpperBounds.push_back( 0.017462  );
  expectedUpperBounds.push_back( 0.017393  );
  expectedUpperBounds.push_back( 0.0154586 );
  expectedUpperBounds.push_back( 0.0263634 );
  expectedUpperBounds.push_back( 0.0199106 );
  expectedUpperBounds.push_back( 0.0175692 );
  expectedUpperBounds.push_back( 0.0178456 );
  expectedUpperBounds.push_back( 0.0263718 );
  expectedUpperBounds.push_back( 0.0135008 );
  expectedUpperBounds.push_back( 0.0274007 );
  expectedUpperBounds.push_back( 0.0274383 );
  expectedUpperBounds.push_back( 0.028699  );
  expectedUpperBounds.push_back( 0.0288704 );
  expectedUpperBounds.push_back( 0.0173995 );
  expectedUpperBounds.push_back( 0.0437453 );
  expectedUpperBounds.push_back( 0.0190598 );
  expectedUpperBounds.push_back( 0.0229029 );
  expectedUpperBounds.push_back( 0.0201335 );
  expectedUpperBounds.push_back( 0.0144457 );
  expectedUpperBounds.push_back( 0.0211218 );
  expectedUpperBounds.push_back( 0.0115863 );
  expectedUpperBounds.push_back( 0.0278765 );
  expectedUpperBounds.push_back( 0.0109632 );
  expectedUpperBounds.push_back( 0.0276103 );
  expectedUpperBounds.push_back( 0.0291187 );
  expectedUpperBounds.push_back( 0.0201573 );
  expectedUpperBounds.push_back( 0.0145157 );
  expectedUpperBounds.push_back( 0.0224448 );
  expectedUpperBounds.push_back( 0.0284639 );
  expectedUpperBounds.push_back( 0.0138714 );
  expectedUpperBounds.push_back( 0.0363775 );
  expectedUpperBounds.push_back( 0.0518937 );
  expectedUpperBounds.push_back( 0.0480965 );
  expectedUpperBounds.push_back( 0.0167523 );
  expectedUpperBounds.push_back( 0.0181142 );
  expectedUpperBounds.push_back( 0.0243036 );
  expectedUpperBounds.push_back( 0.0259994 );
  expectedUpperBounds.push_back( 0.0311689 );
  expectedUpperBounds.push_back( 0.0253329 );
  expectedUpperBounds.push_back( 0.0260988 );
  expectedUpperBounds.push_back( 0.00921062);
  expectedUpperBounds.push_back( 0.0241768 );
  expectedUpperBounds.push_back( 0.0234787 );
  expectedUpperBounds.push_back( 0.0178839 );
  expectedUpperBounds.push_back( 0.0196004 );
  expectedUpperBounds.push_back( 0.0313058 );
  expectedUpperBounds.push_back( 0.0102912 );
  expectedUpperBounds.push_back( 0.0252483 );
  expectedUpperBounds.push_back( 0.0367811 );
  expectedUpperBounds.push_back( 0.0169587 );
  
  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}

void vec_1200_f1p0(){

  vector<double> expectedUpperBounds;

  //for 1200
  expectedUpperBounds.push_back( 0.00295826 );
  expectedUpperBounds.push_back( 0.00630312 );
  expectedUpperBounds.push_back( 0.00367944 );
  expectedUpperBounds.push_back( 0.00388822 );
  expectedUpperBounds.push_back( 0.00434596 );
  expectedUpperBounds.push_back( 0.00545273 );
  expectedUpperBounds.push_back( 0.00538862 );
  expectedUpperBounds.push_back( 0.00726297 );
  expectedUpperBounds.push_back( 0.00862754 );
  expectedUpperBounds.push_back( 0.00539716 );
  expectedUpperBounds.push_back( 0.0111781  );
  expectedUpperBounds.push_back( 0.00623898 );
  expectedUpperBounds.push_back( 0.00432841 );
  expectedUpperBounds.push_back( 0.00811718 );
  expectedUpperBounds.push_back( 0.00565287 );
  expectedUpperBounds.push_back( 0.00951123 );
  expectedUpperBounds.push_back( 0.00734528 );
  expectedUpperBounds.push_back( 0.00731713 );
  expectedUpperBounds.push_back( 0.00635077 );
  expectedUpperBounds.push_back( 0.00595004 );
  expectedUpperBounds.push_back( 0.00356809 );
  expectedUpperBounds.push_back( 0.00656074 );
  expectedUpperBounds.push_back( 0.00468025 );
  expectedUpperBounds.push_back( 0.00506596 );
  expectedUpperBounds.push_back( 0.00490704 );
  expectedUpperBounds.push_back( 0.00943956 );
  expectedUpperBounds.push_back( 0.00505405 );
  expectedUpperBounds.push_back( 0.00621376 );
  expectedUpperBounds.push_back( 0.00333808 );
  expectedUpperBounds.push_back( 0.00723838 );
  expectedUpperBounds.push_back( 0.00942348 );
  expectedUpperBounds.push_back( 0.006292   );
  expectedUpperBounds.push_back( 0.00283945 );
  expectedUpperBounds.push_back( 0.00864033 );
  expectedUpperBounds.push_back( 0.00684278 );
  expectedUpperBounds.push_back( 0.0043875  );
  expectedUpperBounds.push_back( 0.00966373 );
  expectedUpperBounds.push_back( 0.00385247 );
  expectedUpperBounds.push_back( 0.0107606  );
  expectedUpperBounds.push_back( 0.00541071 );
  expectedUpperBounds.push_back( 0.00439581 );
  expectedUpperBounds.push_back( 0.00533984 );
  expectedUpperBounds.push_back( 0.0102658  );
  expectedUpperBounds.push_back( 0.00781691 );
  expectedUpperBounds.push_back( 0.00660997 );
  expectedUpperBounds.push_back( 0.00420674 );
  expectedUpperBounds.push_back( 0.00816178 );
  expectedUpperBounds.push_back( 0.00563417 );
  expectedUpperBounds.push_back( 0.0138786  );
  expectedUpperBounds.push_back( 0.00630876 );
  expectedUpperBounds.push_back( 0.00429129 );
  expectedUpperBounds.push_back( 0.00503815 );
  expectedUpperBounds.push_back( 0.00989202 );
  expectedUpperBounds.push_back( 0.00795169 );
  expectedUpperBounds.push_back( 0.00519256 );
  expectedUpperBounds.push_back( 0.00579568 );
  expectedUpperBounds.push_back( 0.00747241 );
  expectedUpperBounds.push_back( 0.00495549 );
  expectedUpperBounds.push_back( 0.0047214  );
  expectedUpperBounds.push_back( 0.00611978 );
  expectedUpperBounds.push_back( 0.00650526 );
  expectedUpperBounds.push_back( 0.0052604  );
  expectedUpperBounds.push_back( 0.00740148 );
  expectedUpperBounds.push_back( 0.00444231 );
  expectedUpperBounds.push_back( 0.00563374 );
  expectedUpperBounds.push_back( 0.0159705  );
  expectedUpperBounds.push_back( 0.00338417 );
  expectedUpperBounds.push_back( 0.0137241  );
  expectedUpperBounds.push_back( 0.00740456 );
  expectedUpperBounds.push_back( 0.00705698 );
  expectedUpperBounds.push_back( 0.00614694 );
  expectedUpperBounds.push_back( 0.00611207 );
  expectedUpperBounds.push_back( 0.0112424  );
  expectedUpperBounds.push_back( 0.0104042  );
  expectedUpperBounds.push_back( 0.00320295 );
  expectedUpperBounds.push_back( 0.00718351 );
  expectedUpperBounds.push_back( 0.00505542 );
  expectedUpperBounds.push_back( 0.00965458 );
  expectedUpperBounds.push_back( 0.00281405 );
  expectedUpperBounds.push_back( 0.00474659 );
  expectedUpperBounds.push_back( 0.0109049  );
  expectedUpperBounds.push_back( 0.00942392 );
  expectedUpperBounds.push_back( 0.00912152 );
  expectedUpperBounds.push_back( 0.0058679  );
  expectedUpperBounds.push_back( 0.00540342 );
  expectedUpperBounds.push_back( 0.0088658  );
  expectedUpperBounds.push_back( 0.00822416 );
  expectedUpperBounds.push_back( 0.00463496 );
  expectedUpperBounds.push_back( 0.00623127 );
  expectedUpperBounds.push_back( 0.00406823 );
  expectedUpperBounds.push_back( 0.00951913 );
  expectedUpperBounds.push_back( 0.00438157 );
  expectedUpperBounds.push_back( 0.00553212 );
  expectedUpperBounds.push_back( 0.0057864  );
  expectedUpperBounds.push_back( 0.00756349 );
  expectedUpperBounds.push_back( 0.00199818 );
  expectedUpperBounds.push_back( 0.00555939 );
  expectedUpperBounds.push_back( 0.00925514 );
  expectedUpperBounds.push_back( 0.00570553 );
  expectedUpperBounds.push_back( 0.00942681 );
  expectedUpperBounds.push_back( 0.00411137 );
  expectedUpperBounds.push_back( 0.00526604 );
  expectedUpperBounds.push_back( 0.00834102 );
  expectedUpperBounds.push_back( 0.00695762 );
  expectedUpperBounds.push_back( 0.0082437  );
  expectedUpperBounds.push_back( 0.00626601 );
  expectedUpperBounds.push_back( 0.00528922 );
  expectedUpperBounds.push_back( 0.00857217 );
  expectedUpperBounds.push_back( 0.00352909 );
  expectedUpperBounds.push_back( 0.00714064 );
  expectedUpperBounds.push_back( 0.00503811 );
  expectedUpperBounds.push_back( 0.00659258 );
  expectedUpperBounds.push_back( 0.00515337 );
  expectedUpperBounds.push_back( 0.0030673  );
  expectedUpperBounds.push_back( 0.006115   );
  expectedUpperBounds.push_back( 0.00716373 );
  expectedUpperBounds.push_back( 0.00540908 );
  expectedUpperBounds.push_back( 0.0028771  );
  expectedUpperBounds.push_back( 0.00828893 );
  expectedUpperBounds.push_back( 0.00517586 );
  expectedUpperBounds.push_back( 0.00699363 );
  expectedUpperBounds.push_back( 0.0045     );	
  expectedUpperBounds.push_back( 0.00722694 );
  expectedUpperBounds.push_back( 0.00644056 );
  expectedUpperBounds.push_back( 0.00212108 );
  expectedUpperBounds.push_back( 0.0109689  );
  expectedUpperBounds.push_back( 0.00671766 );
  expectedUpperBounds.push_back( 0.00591788 );
  expectedUpperBounds.push_back( 0.00390228 );
  expectedUpperBounds.push_back( 0.00423233 );
  expectedUpperBounds.push_back( 0.00326898 );
  expectedUpperBounds.push_back( 0.00528671 );
  expectedUpperBounds.push_back( 0.0047278  );
  expectedUpperBounds.push_back( 0.00617612 );
  expectedUpperBounds.push_back( 0.00761021 );
  expectedUpperBounds.push_back( 0.0101708  );
  expectedUpperBounds.push_back( 0.00676264 );
  expectedUpperBounds.push_back( 0.00725677 );
  expectedUpperBounds.push_back( 0.00535242 );
  expectedUpperBounds.push_back( 0.00601505 );
  expectedUpperBounds.push_back( 0.00384988 );
  expectedUpperBounds.push_back( 0.011548   );
  expectedUpperBounds.push_back( 0.00769075 );
  expectedUpperBounds.push_back( 0.00534926 );
  expectedUpperBounds.push_back( 0.0111321  );
  expectedUpperBounds.push_back( 0.00335123 );
  expectedUpperBounds.push_back( 0.0066607  );
  expectedUpperBounds.push_back( 0.00725862 );
  expectedUpperBounds.push_back( 0.0045536  );
  expectedUpperBounds.push_back( 0.00413822 );
  expectedUpperBounds.push_back( 0.00993615 );
  expectedUpperBounds.push_back( 0.00247493 );
  expectedUpperBounds.push_back( 0.00739815 );
  expectedUpperBounds.push_back( 0.00938521 );
  expectedUpperBounds.push_back( 0.00363771 );
  expectedUpperBounds.push_back( 0.00534872 );
  expectedUpperBounds.push_back( 0.00525532 );
  expectedUpperBounds.push_back( 0.00628658 );
  expectedUpperBounds.push_back( 0.00789648 );
  expectedUpperBounds.push_back( 0.0039738  );
  expectedUpperBounds.push_back( 0.0042727  );
  expectedUpperBounds.push_back( 0.0053241  );
  expectedUpperBounds.push_back( 0.00538278 );
  expectedUpperBounds.push_back( 0.00774493 );
  expectedUpperBounds.push_back( 0.00444589 );
  expectedUpperBounds.push_back( 0.00734073 );
  expectedUpperBounds.push_back( 0.00547255 );
  expectedUpperBounds.push_back( 0.00546401 );
  expectedUpperBounds.push_back( 0.0053265  );
  expectedUpperBounds.push_back( 0.00641435 );
  expectedUpperBounds.push_back( 0.00857782 );
  expectedUpperBounds.push_back( 0.00494668 );
  expectedUpperBounds.push_back( 0.003816   );
  expectedUpperBounds.push_back( 0.003498   );
  expectedUpperBounds.push_back( 0.00487124 );
  expectedUpperBounds.push_back( 0.00781141 );
  expectedUpperBounds.push_back( 0.00717258 );
  expectedUpperBounds.push_back( 0.00997402 );
  expectedUpperBounds.push_back( 0.00453176 );
  expectedUpperBounds.push_back( 0.00604386 );
  expectedUpperBounds.push_back( 0.00548304 );
  expectedUpperBounds.push_back( 0.00697721 );
  expectedUpperBounds.push_back( 0.0115032  );
  expectedUpperBounds.push_back( 0.00555574 );
  expectedUpperBounds.push_back( 0.00962137 );
  expectedUpperBounds.push_back( 0.00628537 );
  expectedUpperBounds.push_back( 0.00287674 );
  expectedUpperBounds.push_back( 0.00317512 );
  expectedUpperBounds.push_back( 0.00492202 );
  expectedUpperBounds.push_back( 0.0044599  );
  expectedUpperBounds.push_back( 0.00538598 );
  expectedUpperBounds.push_back( 0.0035113  );
  expectedUpperBounds.push_back( 0.0067288  );
  expectedUpperBounds.push_back( 0.00644449 );
  expectedUpperBounds.push_back( 0.0054461  );
  expectedUpperBounds.push_back( 0.00734243 );
  expectedUpperBounds.push_back( 0.00733895 );
  expectedUpperBounds.push_back( 0.00710546 );
  expectedUpperBounds.push_back( 0.00483633 );
  expectedUpperBounds.push_back( 0.00593397 );
  expectedUpperBounds.push_back( 0.00962707 );
  expectedUpperBounds.push_back( 0.00314353 );
  expectedUpperBounds.push_back( 0.00601045 );
  expectedUpperBounds.push_back( 0.00407202 );
  expectedUpperBounds.push_back( 0.00473658 );
  expectedUpperBounds.push_back( 0.00724046 );
  expectedUpperBounds.push_back( 0.00293361 );
  expectedUpperBounds.push_back( 0.0063421  );
  expectedUpperBounds.push_back( 0.00413248 );
  expectedUpperBounds.push_back( 0.00620015 );
  expectedUpperBounds.push_back( 0.00249873 );
  expectedUpperBounds.push_back( 0.00449556 );
  expectedUpperBounds.push_back( 0.00438219 );
  expectedUpperBounds.push_back( 0.00809956 );
  expectedUpperBounds.push_back( 0.00912373 );
  expectedUpperBounds.push_back( 0.00749076 );
  expectedUpperBounds.push_back( 0.00636498 );
  expectedUpperBounds.push_back( 0.00679774 );
  expectedUpperBounds.push_back( 0.00603712 );
  expectedUpperBounds.push_back( 0.00944087 );
  expectedUpperBounds.push_back( 0.00764302 );
  expectedUpperBounds.push_back( 0.00610264 );
  expectedUpperBounds.push_back( 0.00909835 );
  expectedUpperBounds.push_back( 0.00638015 );
  expectedUpperBounds.push_back( 0.00586824 );
  expectedUpperBounds.push_back( 0.00689675 );
  expectedUpperBounds.push_back( 0.00326294 );
  expectedUpperBounds.push_back( 0.00333373 );
  expectedUpperBounds.push_back( 0.00461297 );
  expectedUpperBounds.push_back( 0.008943   );
  expectedUpperBounds.push_back( 0.00483751 );
  expectedUpperBounds.push_back( 0.00907461 );
  expectedUpperBounds.push_back( 0.012548   );
  expectedUpperBounds.push_back( 0.00552243 );
  expectedUpperBounds.push_back( 0.00912074 );
  expectedUpperBounds.push_back( 0.00941578 );


  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}

void vec_3000_f1p0(){

  vector<double> expectedUpperBounds;

  //for 3000
  expectedUpperBounds.push_back( 0.000509995 );
  expectedUpperBounds.push_back( 0.000385683 );
  expectedUpperBounds.push_back( 0.000516278 );
  expectedUpperBounds.push_back( 0.000530295 );
  expectedUpperBounds.push_back( 0.000591133 );
  expectedUpperBounds.push_back( 0.000505589 );
  expectedUpperBounds.push_back( 0.000548643 );
  expectedUpperBounds.push_back( 0.000632181 );
  expectedUpperBounds.push_back( 0.000515813 );
  expectedUpperBounds.push_back( 0.000467964 );
  expectedUpperBounds.push_back( 0.000638667 );
  expectedUpperBounds.push_back( 0.000861837 );
  expectedUpperBounds.push_back( 0.000280555 );
  expectedUpperBounds.push_back( 0.000904176 );
  expectedUpperBounds.push_back( 0.00104824  );
  expectedUpperBounds.push_back( 0.000570718 );
  expectedUpperBounds.push_back( 0.000630994 );
  expectedUpperBounds.push_back( 0.000793137 );
  expectedUpperBounds.push_back( 0.000354215 );
  expectedUpperBounds.push_back( 0.00111414  );
  expectedUpperBounds.push_back( 0.000869649 );
  expectedUpperBounds.push_back( 0.000881929 );
  expectedUpperBounds.push_back( 0.000385222 );
  expectedUpperBounds.push_back( 0.00085481  );
  expectedUpperBounds.push_back( 0.000538348 );
  expectedUpperBounds.push_back( 0.000468887 );
  expectedUpperBounds.push_back( 0.000504794 );
  expectedUpperBounds.push_back( 0.000585425 );
  expectedUpperBounds.push_back( 0.000641127 );
  expectedUpperBounds.push_back( 0.000692037 );
  expectedUpperBounds.push_back( 0.00034995  );
  expectedUpperBounds.push_back( 0.000438798 );
  expectedUpperBounds.push_back( 0.000473284 );
  expectedUpperBounds.push_back( 0.00086658  );
  expectedUpperBounds.push_back( 0.000310442 );
  expectedUpperBounds.push_back( 0.000530372 );
  expectedUpperBounds.push_back( 0.000876681 );
  expectedUpperBounds.push_back( 0.000422246 );
  expectedUpperBounds.push_back( 0.000534517 );
  expectedUpperBounds.push_back( 0.000536165 );
  expectedUpperBounds.push_back( 0.000510264 );
  expectedUpperBounds.push_back( 0.000428096 );
  expectedUpperBounds.push_back( 0.000531173 );
  expectedUpperBounds.push_back( 0.000373337 );
  expectedUpperBounds.push_back( 0.000685651 );
  expectedUpperBounds.push_back( 0.000891929 );
  expectedUpperBounds.push_back( 0.000400768 );
  expectedUpperBounds.push_back( 0.000506942 );
  expectedUpperBounds.push_back( 0.000798447 );
  expectedUpperBounds.push_back( 0.000392092 );
  expectedUpperBounds.push_back( 0.000559053 );
  expectedUpperBounds.push_back( 0.000360086 );
  expectedUpperBounds.push_back( 0.000367961 );
  expectedUpperBounds.push_back( 0.000700531 );
  expectedUpperBounds.push_back( 0.00063966  );
  expectedUpperBounds.push_back( 0.000660095 );
  expectedUpperBounds.push_back( 0.000454756 );
  expectedUpperBounds.push_back( 0.000321433 );
  expectedUpperBounds.push_back( 0.000501883 );
  expectedUpperBounds.push_back( 0.000453682 );
  expectedUpperBounds.push_back( 0.000615641 );
  expectedUpperBounds.push_back( 0.00108412  );
  expectedUpperBounds.push_back( 0.000457688 );
  expectedUpperBounds.push_back( 0.000800954 );
  expectedUpperBounds.push_back( 0.00058839  );
  expectedUpperBounds.push_back( 0.00059773  );
  expectedUpperBounds.push_back( 0.000332055 );
  expectedUpperBounds.push_back( 0.000479548 );
  expectedUpperBounds.push_back( 0.000466172 );
  expectedUpperBounds.push_back( 0.000291761 );
  expectedUpperBounds.push_back( 0.00068389  );
  expectedUpperBounds.push_back( 0.000619756 );
  expectedUpperBounds.push_back( 0.000496065 );
  expectedUpperBounds.push_back( 0.000498372 );
  expectedUpperBounds.push_back( 0.00073567  );
  expectedUpperBounds.push_back( 0.000375689 );
  expectedUpperBounds.push_back( 0.000886455 );
  expectedUpperBounds.push_back( 0.000513442 );
  expectedUpperBounds.push_back( 0.000408562 );
  expectedUpperBounds.push_back( 0.000858832 );
  expectedUpperBounds.push_back( 0.000701447 );
  expectedUpperBounds.push_back( 0.00048229  );
  expectedUpperBounds.push_back( 0.000344477 );
  expectedUpperBounds.push_back( 0.000390722 );
  expectedUpperBounds.push_back( 0.000505789 );
  expectedUpperBounds.push_back( 0.000730323 );
  expectedUpperBounds.push_back( 0.000938888 );
  expectedUpperBounds.push_back( 0.00064487  );
  expectedUpperBounds.push_back( 0.000429202 );
  expectedUpperBounds.push_back( 0.000597202 );
  expectedUpperBounds.push_back( 0.000220266 );
  expectedUpperBounds.push_back( 0.000769021 );
  expectedUpperBounds.push_back( 0.000537855 );
  expectedUpperBounds.push_back( 0.00064488  );
  expectedUpperBounds.push_back( 0.000505378 );
  expectedUpperBounds.push_back( 0.000457449 );
  expectedUpperBounds.push_back( 0.000466554 );
  expectedUpperBounds.push_back( 0.000554225 );
  expectedUpperBounds.push_back( 0.000691167 );
  expectedUpperBounds.push_back( 0.000474157 );
  expectedUpperBounds.push_back( 0.000578515 );
  expectedUpperBounds.push_back( 0.000403801 );
  expectedUpperBounds.push_back( 0.000500124 );
  expectedUpperBounds.push_back( 0.0010566   );
  expectedUpperBounds.push_back( 0.000936766 );
  expectedUpperBounds.push_back( 0.000800487 );
  expectedUpperBounds.push_back( 0.000421873 );
  expectedUpperBounds.push_back( 0.000575153 );
  expectedUpperBounds.push_back( 0.00048848  );
  expectedUpperBounds.push_back( 0.000391331 );
  expectedUpperBounds.push_back( 0.000619846 );
  expectedUpperBounds.push_back( 0.000976432 );
  expectedUpperBounds.push_back( 0.000665387 );
  expectedUpperBounds.push_back( 0.000548494 );
  expectedUpperBounds.push_back( 0.000447988 );
  expectedUpperBounds.push_back( 0.000407646 );
  expectedUpperBounds.push_back( 0.000953943 );
  expectedUpperBounds.push_back( 0.000630741 );
  expectedUpperBounds.push_back( 0.000847271 );
  expectedUpperBounds.push_back( 0.000512333 );
  expectedUpperBounds.push_back( 0.000517054 );
  expectedUpperBounds.push_back( 0.000964687 );
  expectedUpperBounds.push_back( 0.000750879 );
  expectedUpperBounds.push_back( 0.000525082 );
  expectedUpperBounds.push_back( 0.000868687 );
  expectedUpperBounds.push_back( 0.000517839 );
  expectedUpperBounds.push_back( 0.000368791 );
  expectedUpperBounds.push_back( 0.000437316 );
  expectedUpperBounds.push_back( 0.000516622 );
  expectedUpperBounds.push_back( 0.000542514 );
  expectedUpperBounds.push_back( 0.000735314 );
  expectedUpperBounds.push_back( 0.000533846 );
  expectedUpperBounds.push_back( 0.000638052 );
  expectedUpperBounds.push_back( 0.000334841 );
  expectedUpperBounds.push_back( 0.000393953 );
  expectedUpperBounds.push_back( 0.000280974 );
  expectedUpperBounds.push_back( 0.000860275 );
  expectedUpperBounds.push_back( 0.000382168 );
  expectedUpperBounds.push_back( 0.000544895 );
  expectedUpperBounds.push_back( 0.000417322 );
  expectedUpperBounds.push_back( 0.000374809 );
  expectedUpperBounds.push_back( 0.00033246  );
  expectedUpperBounds.push_back( 0.000584972 );
  expectedUpperBounds.push_back( 0.000387199 );
  expectedUpperBounds.push_back( 0.000566477 );
  expectedUpperBounds.push_back( 0.000351404 );
  expectedUpperBounds.push_back( 0.00065211  );
  expectedUpperBounds.push_back( 0.000462423 );
  expectedUpperBounds.push_back( 0.000426822 );
  expectedUpperBounds.push_back( 0.00108849  );
  expectedUpperBounds.push_back( 0.000639355 );
  expectedUpperBounds.push_back( 0.000441454 );
  expectedUpperBounds.push_back( 0.000451294 );
  expectedUpperBounds.push_back( 0.000651022 );
  expectedUpperBounds.push_back( 0.000492952 );
  expectedUpperBounds.push_back( 0.000564205 );
  expectedUpperBounds.push_back( 0.000474429 );
  expectedUpperBounds.push_back( 0.000382291 );
  expectedUpperBounds.push_back( 0.000473644 );
  expectedUpperBounds.push_back( 0.000731222 );
  expectedUpperBounds.push_back( 0.00080686  );
  expectedUpperBounds.push_back( 0.000522959 );
  expectedUpperBounds.push_back( 0.000780012 );
  expectedUpperBounds.push_back( 0.000287214 );
  expectedUpperBounds.push_back( 0.000871049 );
  expectedUpperBounds.push_back( 0.000212181 );
  expectedUpperBounds.push_back( 0.000376087 );
  expectedUpperBounds.push_back( 0.000570962 );
  expectedUpperBounds.push_back( 0.000682877 );
  expectedUpperBounds.push_back( 0.000574201 );
  expectedUpperBounds.push_back( 0.000396704 );
  expectedUpperBounds.push_back( 0.000377122 );
  expectedUpperBounds.push_back( 0.000821005 );
  expectedUpperBounds.push_back( 0.00036972  );
  expectedUpperBounds.push_back( 0.000697686 );
  expectedUpperBounds.push_back( 0.000561645 );
  expectedUpperBounds.push_back( 0.000754717 );
  expectedUpperBounds.push_back( 0.000936805 );
  expectedUpperBounds.push_back( 0.000385278 );
  expectedUpperBounds.push_back( 0.000510133 );
  expectedUpperBounds.push_back( 0.000706377 );
  expectedUpperBounds.push_back( 0.0012726   );
  expectedUpperBounds.push_back( 0.000671869 );
  expectedUpperBounds.push_back( 0.000982412 );
  expectedUpperBounds.push_back( 0.000797179 );
  expectedUpperBounds.push_back( 0.000837539 );
  expectedUpperBounds.push_back( 0.000382836 );
  expectedUpperBounds.push_back( 0.000538244 );
  expectedUpperBounds.push_back( 0.000656867 );
  expectedUpperBounds.push_back( 0.000280834 );
  expectedUpperBounds.push_back( 0.000386611 );
  expectedUpperBounds.push_back( 0.000364294 );
  expectedUpperBounds.push_back( 0.00128854  );
  expectedUpperBounds.push_back( 0.000346362 );
  expectedUpperBounds.push_back( 0.000912378 );
  expectedUpperBounds.push_back( 0.000428437 );
  expectedUpperBounds.push_back( 0.00068881  );
  expectedUpperBounds.push_back( 0.000610683 );
  expectedUpperBounds.push_back( 0.000892912 );
  expectedUpperBounds.push_back( 0.000562486 );
  expectedUpperBounds.push_back( 0.000394644 );
  expectedUpperBounds.push_back( 0.00047455  );
  expectedUpperBounds.push_back( 0.000528234 );
  expectedUpperBounds.push_back( 0.000676502 );
  expectedUpperBounds.push_back( 0.000618298 );
  expectedUpperBounds.push_back( 0.000913659 );
  expectedUpperBounds.push_back( 0.000684474 );
  expectedUpperBounds.push_back( 0.000708822 );
  expectedUpperBounds.push_back( 0.00106556  );
  expectedUpperBounds.push_back( 0.000382132 );
  expectedUpperBounds.push_back( 0.000666357 );
  expectedUpperBounds.push_back( 0.000352691 );
  expectedUpperBounds.push_back( 0.0011981   );
  expectedUpperBounds.push_back( 0.000903553 );
  expectedUpperBounds.push_back( 0.000451579 );
  expectedUpperBounds.push_back( 0.000725631 );
  expectedUpperBounds.push_back( 0.000442866 );
  expectedUpperBounds.push_back( 0.000444186 );
  expectedUpperBounds.push_back( 0.000349733 );
  expectedUpperBounds.push_back( 0.000412408 );
  expectedUpperBounds.push_back( 0.000964238 );
  expectedUpperBounds.push_back( 0.000671693 );
  expectedUpperBounds.push_back( 0.000623584 );
  expectedUpperBounds.push_back( 0.000521258 );
  expectedUpperBounds.push_back( 0.000402494 );
  expectedUpperBounds.push_back( 0.000721996 );
  expectedUpperBounds.push_back( 0.000605963 );
  expectedUpperBounds.push_back( 0.000486161 );

  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}

void vec_4000_f1p0(){

  vector<double> expectedUpperBounds;

  //for 4000
  expectedUpperBounds.push_back( 0.000428594 );
  expectedUpperBounds.push_back( 0.000626615 );
  expectedUpperBounds.push_back( 0.000428776 );
  expectedUpperBounds.push_back( 0.000573807 );
  expectedUpperBounds.push_back( 0.000678109 );
  expectedUpperBounds.push_back( 0.000612673 );
  expectedUpperBounds.push_back( 0.000730578 );
  expectedUpperBounds.push_back( 0.000818145 );
  expectedUpperBounds.push_back( 0.000572129 );
  expectedUpperBounds.push_back( 0.000371621 );
  expectedUpperBounds.push_back( 0.000370114 );
  expectedUpperBounds.push_back( 0.000761789 );
  expectedUpperBounds.push_back( 0.000463118 );
  expectedUpperBounds.push_back( 0.000556552 );
  expectedUpperBounds.push_back( 0.000465295 );
  expectedUpperBounds.push_back( 0.000473486 );
  expectedUpperBounds.push_back( 0.000408154 );
  expectedUpperBounds.push_back( 0.000528419 );
  expectedUpperBounds.push_back( 0.000416503 );
  expectedUpperBounds.push_back( 0.000511934 );
  expectedUpperBounds.push_back( 0.000455277 );
  expectedUpperBounds.push_back( 0.00052637  );
  expectedUpperBounds.push_back( 0.000636506 );
  expectedUpperBounds.push_back( 0.000436845 );
  expectedUpperBounds.push_back( 0.000897239 );
  expectedUpperBounds.push_back( 0.000416287 );
  expectedUpperBounds.push_back( 0.000571086 );
  expectedUpperBounds.push_back( 0.001434    );
  expectedUpperBounds.push_back( 0.000524059 );
  expectedUpperBounds.push_back( 0.00072818  );
  expectedUpperBounds.push_back( 0.000516407 );
  expectedUpperBounds.push_back( 0.000430884 );
  expectedUpperBounds.push_back( 0.000737917 );
  expectedUpperBounds.push_back( 0.000516262 );
  expectedUpperBounds.push_back( 0.000497772 );
  expectedUpperBounds.push_back( 0.000499187 );
  expectedUpperBounds.push_back( 0.000469617 );
  expectedUpperBounds.push_back( 0.000577076 );
  expectedUpperBounds.push_back( 0.000618101 );
  expectedUpperBounds.push_back( 0.000711313 );
  expectedUpperBounds.push_back( 0.000654982 );
  expectedUpperBounds.push_back( 0.00055505  );
  expectedUpperBounds.push_back( 0.000483352 );
  expectedUpperBounds.push_back( 0.000371006 );
  expectedUpperBounds.push_back( 0.000342676 );
  expectedUpperBounds.push_back( 0.000616953 );
  expectedUpperBounds.push_back( 0.000797422 );
  expectedUpperBounds.push_back( 0.000399494 );
  expectedUpperBounds.push_back( 0.000548835 );
  expectedUpperBounds.push_back( 0.000506763 );
  expectedUpperBounds.push_back( 0.000625614 );
  expectedUpperBounds.push_back( 0.000676586 );
  expectedUpperBounds.push_back( 0.00043752  );
  expectedUpperBounds.push_back( 0.000472404 );
  expectedUpperBounds.push_back( 0.00093672  );
  expectedUpperBounds.push_back( 0.000636638 );
  expectedUpperBounds.push_back( 0.000471864 );
  expectedUpperBounds.push_back( 0.0005345   );
  expectedUpperBounds.push_back( 0.000414484 );
  expectedUpperBounds.push_back( 0.00045184  );
  expectedUpperBounds.push_back( 0.000450883 );
  expectedUpperBounds.push_back( 0.000559106 );
  expectedUpperBounds.push_back( 0.000455022 );
  expectedUpperBounds.push_back( 0.00031356  );
  expectedUpperBounds.push_back( 0.00056667  );
  expectedUpperBounds.push_back( 0.000804911 );
  expectedUpperBounds.push_back( 0.000695465 );
  expectedUpperBounds.push_back( 0.000524728 );
  expectedUpperBounds.push_back( 0.000510722 );
  expectedUpperBounds.push_back( 0.000600732 );
  expectedUpperBounds.push_back( 0.000273227 );
  expectedUpperBounds.push_back( 0.000566175 );
  expectedUpperBounds.push_back( 0.000467599 );
  expectedUpperBounds.push_back( 0.00049265  );
  expectedUpperBounds.push_back( 0.000542731 );
  expectedUpperBounds.push_back( 0.000671525 );
  expectedUpperBounds.push_back( 0.000537856 );
  expectedUpperBounds.push_back( 0.000495357 );
  expectedUpperBounds.push_back( 0.000542526 );
  expectedUpperBounds.push_back( 0.000366398 );
  expectedUpperBounds.push_back( 0.000475446 );
  expectedUpperBounds.push_back( 0.00046659  );
  expectedUpperBounds.push_back( 0.000463479 );
  expectedUpperBounds.push_back( 0.000650867 );
  expectedUpperBounds.push_back( 0.000672163 );
  expectedUpperBounds.push_back( 0.000457534 );
  expectedUpperBounds.push_back( 0.000416598 );
  expectedUpperBounds.push_back( 0.000526966 );
  expectedUpperBounds.push_back( 0.000514285 );
  expectedUpperBounds.push_back( 0.000473462 );
  expectedUpperBounds.push_back( 0.00104674  );
  expectedUpperBounds.push_back( 0.00037506  );
  expectedUpperBounds.push_back( 0.00054965  );
  expectedUpperBounds.push_back( 0.000373716 );
  expectedUpperBounds.push_back( 0.000516154 );
  expectedUpperBounds.push_back( 0.000438945 );
  expectedUpperBounds.push_back( 0.000507497 );
  expectedUpperBounds.push_back( 0.000768304 );
  expectedUpperBounds.push_back( 0.000449509 );
  expectedUpperBounds.push_back( 0.000542804 );
  expectedUpperBounds.push_back( 0.000517118 );
  expectedUpperBounds.push_back( 0.000649398 );
  expectedUpperBounds.push_back( 0.000462917 );
  expectedUpperBounds.push_back( 0.000611521 );
  expectedUpperBounds.push_back( 0.000415409 );
  expectedUpperBounds.push_back( 0.000494847 );
  expectedUpperBounds.push_back( 0.000340034 );
  expectedUpperBounds.push_back( 0.000559695 );
  expectedUpperBounds.push_back( 0.000364658 );
  expectedUpperBounds.push_back( 0.000545315 );
  expectedUpperBounds.push_back( 0.000569141 );
  expectedUpperBounds.push_back( 0.000550054 );
  expectedUpperBounds.push_back( 0.000605736 );
  expectedUpperBounds.push_back( 0.00058196  );
  expectedUpperBounds.push_back( 0.000556875 );
  expectedUpperBounds.push_back( 0.000737528 );
  expectedUpperBounds.push_back( 0.000865728 );
  expectedUpperBounds.push_back( 0.0005194   );
  expectedUpperBounds.push_back( 0.000536016 );
  expectedUpperBounds.push_back( 0.000551367 );
  expectedUpperBounds.push_back( 0.000496103 );
  expectedUpperBounds.push_back( 0.000469437 );
  expectedUpperBounds.push_back( 0.000448063 );
  expectedUpperBounds.push_back( 0.000418085 );
  expectedUpperBounds.push_back( 0.000332565 );
  expectedUpperBounds.push_back( 0.000367322 );
  expectedUpperBounds.push_back( 0.000464483 );
  expectedUpperBounds.push_back( 0.000547208 );
  expectedUpperBounds.push_back( 0.000382271 );
  expectedUpperBounds.push_back( 0.000739992 );
  expectedUpperBounds.push_back( 0.000388673 );
  expectedUpperBounds.push_back( 0.000396455 );
  expectedUpperBounds.push_back( 0.0006455   );
  expectedUpperBounds.push_back( 0.000518417 );
  expectedUpperBounds.push_back( 0.000573714 );
  expectedUpperBounds.push_back( 0.000701155 );
  expectedUpperBounds.push_back( 0.000521257 );
  expectedUpperBounds.push_back( 0.000528149 );
  expectedUpperBounds.push_back( 0.000330324 );
  expectedUpperBounds.push_back( 0.000534584 );
  expectedUpperBounds.push_back( 0.000663002 );
  expectedUpperBounds.push_back( 0.000549276 );
  expectedUpperBounds.push_back( 0.000595821 );
  expectedUpperBounds.push_back( 0.000385471 );
  expectedUpperBounds.push_back( 0.00064678  );
  expectedUpperBounds.push_back( 0.000382754 );
  expectedUpperBounds.push_back( 0.000478921 );
  expectedUpperBounds.push_back( 0.000501364 );
  expectedUpperBounds.push_back( 0.000504307 );
  expectedUpperBounds.push_back( 0.000761121 );
  expectedUpperBounds.push_back( 0.000665075 );
  expectedUpperBounds.push_back( 0.000405416 );
  expectedUpperBounds.push_back( 0.000637322 );
  expectedUpperBounds.push_back( 0.000372863 );
  expectedUpperBounds.push_back( 0.000358371 );
  expectedUpperBounds.push_back( 0.000604228 );
  expectedUpperBounds.push_back( 0.000546353 );
  expectedUpperBounds.push_back( 0.000400892 );
  expectedUpperBounds.push_back( 0.000858137 );
  expectedUpperBounds.push_back( 0.000513872 );
  expectedUpperBounds.push_back( 0.000669237 );
  expectedUpperBounds.push_back( 0.000485437 );
  expectedUpperBounds.push_back( 0.000460078 );
  expectedUpperBounds.push_back( 0.000530049 );
  expectedUpperBounds.push_back( 0.000536738 );
  expectedUpperBounds.push_back( 0.00104538  );
  expectedUpperBounds.push_back( 0.00044486  );
  expectedUpperBounds.push_back( 0.000453735 );
  expectedUpperBounds.push_back( 0.000598486 );
  expectedUpperBounds.push_back( 0.000657443 );
  expectedUpperBounds.push_back( 0.000524635 );
  expectedUpperBounds.push_back( 0.000507789 );
  expectedUpperBounds.push_back( 0.00046838  );
  expectedUpperBounds.push_back( 0.00036962  );
  expectedUpperBounds.push_back( 0.000445208 );
  expectedUpperBounds.push_back( 0.000658868 );
  expectedUpperBounds.push_back( 0.00054941  );
  expectedUpperBounds.push_back( 0.000569687 );
  expectedUpperBounds.push_back( 0.000341035 );
  expectedUpperBounds.push_back( 0.000368318 );
  expectedUpperBounds.push_back( 0.000527912 );
  expectedUpperBounds.push_back( 0.000551576 );
  expectedUpperBounds.push_back( 0.000482072 );
  expectedUpperBounds.push_back( 0.000373473 );
  expectedUpperBounds.push_back( 0.000574942 );
  expectedUpperBounds.push_back( 0.00101674  );
  expectedUpperBounds.push_back( 0.000495577 );
  expectedUpperBounds.push_back( 0.000579633 );
  expectedUpperBounds.push_back( 0.000704538 );
  expectedUpperBounds.push_back( 0.000382014 );
  expectedUpperBounds.push_back( 0.000536671 );


  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}

void vec_700_f0p5(){

  vector<double> expectedUpperBounds;




  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}


void vec_3000_f0p5(){

  vector<double> expectedUpperBounds;

  //for 3000
  expectedUpperBounds.push_back( 0.0002909   );
  expectedUpperBounds.push_back( 0.000692701 );
  expectedUpperBounds.push_back( 0.00047381  );
  expectedUpperBounds.push_back( 0.0007998   );
  expectedUpperBounds.push_back( 0.000340142 );
  expectedUpperBounds.push_back( 0.000490616 );
  expectedUpperBounds.push_back( 0.000510403 );
  expectedUpperBounds.push_back( 0.000274629 );
  expectedUpperBounds.push_back( 0.000538764 );
  expectedUpperBounds.push_back( 0.000369982 );
  expectedUpperBounds.push_back( 0.00033351  );
  expectedUpperBounds.push_back( 0.000582693 );
  expectedUpperBounds.push_back( 0.000721869 );
  expectedUpperBounds.push_back( 0.000490788 );
  expectedUpperBounds.push_back( 0.000633982 );
  expectedUpperBounds.push_back( 0.000334795 );
  expectedUpperBounds.push_back( 0.000533001 );
  expectedUpperBounds.push_back( 0.000235724 );
  expectedUpperBounds.push_back( 0.000377585 );
  expectedUpperBounds.push_back( 0.000493438 );
  expectedUpperBounds.push_back( 0.000739977 );
  expectedUpperBounds.push_back( 0.000215998 );
  expectedUpperBounds.push_back( 0.000286596 );
  expectedUpperBounds.push_back( 0.000590261 );
  expectedUpperBounds.push_back( 0.000172987 );
  expectedUpperBounds.push_back( 0.000359321 );
  expectedUpperBounds.push_back( 0.000454033 );
  expectedUpperBounds.push_back( 0.000381519 );
  expectedUpperBounds.push_back( 0.000460133 );
  expectedUpperBounds.push_back( 0.000405986 );
  expectedUpperBounds.push_back( 0.000592146 );
  expectedUpperBounds.push_back( 0.000637323 );
  expectedUpperBounds.push_back( 0.000611879 );
  expectedUpperBounds.push_back( 0.00047368  );
  expectedUpperBounds.push_back( 0.000294743 );
  expectedUpperBounds.push_back( 0.000474184 );
  expectedUpperBounds.push_back( 0.000459674 );
  expectedUpperBounds.push_back( 0.000417297 );
  expectedUpperBounds.push_back( 0.000525983 );
  expectedUpperBounds.push_back( 0.000720498 );
  expectedUpperBounds.push_back( 0.000404645 );
  expectedUpperBounds.push_back( 0.000359549 );
  expectedUpperBounds.push_back( 0.000433126 );
  expectedUpperBounds.push_back( 0.000365572 );
  expectedUpperBounds.push_back( 0.00032263  );
  expectedUpperBounds.push_back( 0.000359737 );
  expectedUpperBounds.push_back( 0.000474556 );
  expectedUpperBounds.push_back( 0.000541477 );
  expectedUpperBounds.push_back( 0.000226878 );
  expectedUpperBounds.push_back( 0.000305718 );
  expectedUpperBounds.push_back( 0.000483267 );
  expectedUpperBounds.push_back( 0.000529629 );
  expectedUpperBounds.push_back( 0.00049518  );
  expectedUpperBounds.push_back( 0.000486638 );
  expectedUpperBounds.push_back( 0.000708461 );
  expectedUpperBounds.push_back( 0.00038064  );
  expectedUpperBounds.push_back( 0.000301787 );
  expectedUpperBounds.push_back( 0.000517606 );
  expectedUpperBounds.push_back( 0.000347171 );
  expectedUpperBounds.push_back( 0.000498154 );
  expectedUpperBounds.push_back( 0.000515194 );
  expectedUpperBounds.push_back( 0.000337659 );
  expectedUpperBounds.push_back( 0.000391927 );
  expectedUpperBounds.push_back( 0.000546047 );
  expectedUpperBounds.push_back( 0.00031685  );
  expectedUpperBounds.push_back( 0.00115646  );
  expectedUpperBounds.push_back( 0.00034319  );
  expectedUpperBounds.push_back( 0.00061182  );
  expectedUpperBounds.push_back( 0.00046424  );
  expectedUpperBounds.push_back( 0.000609411 );
  expectedUpperBounds.push_back( 0.000323985 );
  expectedUpperBounds.push_back( 0.000389585 );
  expectedUpperBounds.push_back( 0.000372654 );
  expectedUpperBounds.push_back( 0.000273704 );
  expectedUpperBounds.push_back( 0.000345405 );
  expectedUpperBounds.push_back( 0.00026287  );
  expectedUpperBounds.push_back( 0.000348976 );
  expectedUpperBounds.push_back( 0.000295939 );
  expectedUpperBounds.push_back( 0.000480424 );
  expectedUpperBounds.push_back( 0.000511761 );
  expectedUpperBounds.push_back( 0.000400682 );
  expectedUpperBounds.push_back( 0.000604848 );
  expectedUpperBounds.push_back( 0.000677089 );
  expectedUpperBounds.push_back( 0.000330637 );
  expectedUpperBounds.push_back( 0.000345169 );
  expectedUpperBounds.push_back( 0.000473187 );
  expectedUpperBounds.push_back( 0.000267397 );
  expectedUpperBounds.push_back( 0.000398274 );
  expectedUpperBounds.push_back( 0.000349019 );
  expectedUpperBounds.push_back( 0.000469473 );
  expectedUpperBounds.push_back( 0.000305982 );
  expectedUpperBounds.push_back( 0.000561057 );
  expectedUpperBounds.push_back( 0.000348    );
  expectedUpperBounds.push_back( 0.000559543 );
  expectedUpperBounds.push_back( 0.000385555 );
  expectedUpperBounds.push_back( 0.000618257 );
  expectedUpperBounds.push_back( 0.00100821  );
  expectedUpperBounds.push_back( 0.000483795 );
  expectedUpperBounds.push_back( 0.00039928  );
  expectedUpperBounds.push_back( 0.000445292 );
  expectedUpperBounds.push_back( 0.000442534 );
  expectedUpperBounds.push_back( 0.000675543 );
  expectedUpperBounds.push_back( 0.000490782 );
  expectedUpperBounds.push_back( 0.000293086 );
  expectedUpperBounds.push_back( 0.000504434 );
  expectedUpperBounds.push_back( 0.000320921 );
  expectedUpperBounds.push_back( 0.000439786 );
  expectedUpperBounds.push_back( 0.000437965 );
  expectedUpperBounds.push_back( 0.000425457 );
  expectedUpperBounds.push_back( 0.000428955 );
  expectedUpperBounds.push_back( 0.000363161 );
  expectedUpperBounds.push_back( 0.000183505 );
  expectedUpperBounds.push_back( 0.000555362 );
  expectedUpperBounds.push_back( 0.000544465 );
  expectedUpperBounds.push_back( 0.000241042 );
  expectedUpperBounds.push_back( 0.000350912 );
  expectedUpperBounds.push_back( 0.000713296 );
  expectedUpperBounds.push_back( 0.000560947 );
  expectedUpperBounds.push_back( 0.000518284 );
  expectedUpperBounds.push_back( 0.00030613  );
  expectedUpperBounds.push_back( 0.000336544 );
  expectedUpperBounds.push_back( 0.000484145 );
  expectedUpperBounds.push_back( 0.000307758 );
  expectedUpperBounds.push_back( 0.000464694 );
  expectedUpperBounds.push_back( 0.000477361 );
  expectedUpperBounds.push_back( 0.000464815 );
  expectedUpperBounds.push_back( 0.00037395  );
  expectedUpperBounds.push_back( 0.000424114 );
  expectedUpperBounds.push_back( 0.000625015 );
  expectedUpperBounds.push_back( 0.000394403 );
  expectedUpperBounds.push_back( 0.000375118 );
  expectedUpperBounds.push_back( 0.000550081 );
  expectedUpperBounds.push_back( 0.000298561 );
  expectedUpperBounds.push_back( 0.000280377 );
  expectedUpperBounds.push_back( 0.000559448 );
  expectedUpperBounds.push_back( 0.000678539 );
  expectedUpperBounds.push_back( 0.000308218 );
  expectedUpperBounds.push_back( 0.000337838 );
  expectedUpperBounds.push_back( 0.000254715 );
  expectedUpperBounds.push_back( 0.000412298 );
  expectedUpperBounds.push_back( 0.00029021  );
  expectedUpperBounds.push_back( 0.000250353 );
  expectedUpperBounds.push_back( 0.000531797 );
  expectedUpperBounds.push_back( 0.000428232 );
  expectedUpperBounds.push_back( 0.000405131 );
  expectedUpperBounds.push_back( 0.000322766 );
  expectedUpperBounds.push_back( 0.000429539 );
  expectedUpperBounds.push_back( 0.00057305  );
  expectedUpperBounds.push_back( 0.000472901 );
  expectedUpperBounds.push_back( 0.000440471 );
  expectedUpperBounds.push_back( 0.000298691 );
  expectedUpperBounds.push_back( 0.00039167  );
  expectedUpperBounds.push_back( 0.000314666 );
  expectedUpperBounds.push_back( 0.000453419 );
  expectedUpperBounds.push_back( 0.000354814 );
  expectedUpperBounds.push_back( 0.000891242 );
  expectedUpperBounds.push_back( 0.000438324 );
  expectedUpperBounds.push_back( 0.000436385 );
  expectedUpperBounds.push_back( 0.000760707 );
  expectedUpperBounds.push_back( 0.00051714  );
  expectedUpperBounds.push_back( 0.00034409  );
  expectedUpperBounds.push_back( 0.000332182 );
  expectedUpperBounds.push_back( 0.000347926 );
  expectedUpperBounds.push_back( 0.000484112 );
  expectedUpperBounds.push_back( 0.000577952 );
  expectedUpperBounds.push_back( 0.000466006 );
  expectedUpperBounds.push_back( 0.0003483   );
  expectedUpperBounds.push_back( 0.000723571 );
  expectedUpperBounds.push_back( 0.000505744 );
  expectedUpperBounds.push_back( 0.000287338 );
  expectedUpperBounds.push_back( 0.000758352 );
  expectedUpperBounds.push_back( 0.000471693 );
  expectedUpperBounds.push_back( 0.000396031 );
  expectedUpperBounds.push_back( 0.000232923 );
  expectedUpperBounds.push_back( 0.000425345 );
  expectedUpperBounds.push_back( 0.000645054 );
  expectedUpperBounds.push_back( 0.000908162 );
  expectedUpperBounds.push_back( 0.000635896 );
  expectedUpperBounds.push_back( 0.000857452 );
  expectedUpperBounds.push_back( 0.000519862 );
  expectedUpperBounds.push_back( 0.000440864 );
  expectedUpperBounds.push_back( 0.000278815 );
  expectedUpperBounds.push_back( 0.000313473 );
  expectedUpperBounds.push_back( 0.000290782 );
  expectedUpperBounds.push_back( 0.000441306 );
  expectedUpperBounds.push_back( 0.00033972  );
  expectedUpperBounds.push_back( 0.0002316   );
  expectedUpperBounds.push_back( 0.000579158 );
  expectedUpperBounds.push_back( 0.000308858 );
  expectedUpperBounds.push_back( 0.00080281  );
  expectedUpperBounds.push_back( 0.000395244 );
  expectedUpperBounds.push_back( 0.000255841 );
  expectedUpperBounds.push_back( 0.000329891 );
  expectedUpperBounds.push_back( 0.000416828 );
  expectedUpperBounds.push_back( 0.000385774 );
  expectedUpperBounds.push_back( 0.000546755 );
  expectedUpperBounds.push_back( 0.00112617  );
  expectedUpperBounds.push_back( 0.000516261 );
  expectedUpperBounds.push_back( 0.000657434 );
  expectedUpperBounds.push_back( 0.00035786  );
  expectedUpperBounds.push_back( 0.000294151 );
  expectedUpperBounds.push_back( 0.000263981 );
  expectedUpperBounds.push_back( 0.000173656 );
  expectedUpperBounds.push_back( 0.000254115 );
  expectedUpperBounds.push_back( 0.000431183 );
  expectedUpperBounds.push_back( 0.00027686  );
  expectedUpperBounds.push_back( 0.000469623 );
  expectedUpperBounds.push_back( 0.000334625 );
  expectedUpperBounds.push_back( 0.000599438 );
  expectedUpperBounds.push_back( 0.000559344 );
  expectedUpperBounds.push_back( 0.000667384 );
  expectedUpperBounds.push_back( 0.000708757 );
  expectedUpperBounds.push_back( 0.000597977 );
  expectedUpperBounds.push_back( 0.000455907 );
  expectedUpperBounds.push_back( 0.000260565 );
  expectedUpperBounds.push_back( 0.000868137 );
  expectedUpperBounds.push_back( 0.000461611 );
  expectedUpperBounds.push_back( 0.000259795 );
  expectedUpperBounds.push_back( 0.000543469 );
  expectedUpperBounds.push_back( 0.000285339 );
  expectedUpperBounds.push_back( 0.000573069 );
  expectedUpperBounds.push_back( 0.000557327 );
  expectedUpperBounds.push_back( 0.000405287 );
  expectedUpperBounds.push_back( 0.000536636 );
  expectedUpperBounds.push_back( 0.000385867 );
  expectedUpperBounds.push_back( 0.000347458 );
  expectedUpperBounds.push_back( 0.000405027 );
  expectedUpperBounds.push_back( 0.000489591 );
  expectedUpperBounds.push_back( 0.000251465 );
  expectedUpperBounds.push_back( 0.000593678 );


  double median;
  pair<double, double> onesigma;
  pair<double, double> twosigma;

  getQuantiles(expectedUpperBounds, median, onesigma, twosigma);
  cout << "median: " << median << endl;
  cout << "+/-1 sigma band: [ " << onesigma.first << " , " << onesigma.second << " ] " << endl;
  cout << "+/-2 sigma band: [ " << twosigma.first << " , " << twosigma.second << " ] " << endl;

}
