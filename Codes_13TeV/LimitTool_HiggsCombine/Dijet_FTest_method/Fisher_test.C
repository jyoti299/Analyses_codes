#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TArrow.h"
#include "TLine.h"
#include "TPad.h"
#include "TAxis.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
//#include <iomanip>                                                                                                                                   

using namespace std;
using namespace ROOT;


double RSS(TH1D*, double, double, double, double, double, double, double,
           double);
double getN(TH1D*, double, double);

void Fisher_test()
{

    double lumi = 36813.0; // Integrated luminosity [/pb]

    // Fit range [GeV]
    double minX = 695.0;
    double maxX = 5000.0;

    TFile *f_in;
    TH1F *h_data_1GeV;
    f_in = new TFile("../ExcitedQuarks/RootTreeAnalyzer/inputs/Bias_study/TotalMC_Run2016_ReReco-BCDEFG_PromptReco-H_80X_36813pb_Cut-PhLID_JetTID_Nodeta_Nodphi_CSVM_Mass695_QstarInvtMass.root");
    f_in->GetObject("h_mgj_13TeV_1GeVBin", h_data_1GeV);

    const int nMassBins = 119;
    double massBoundaries[nMassBins+1]
        = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 
           252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 
           731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 
           1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 
           2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 
           3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 
           6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 
           9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};
                           
    TH1D *h_data = static_cast<TH1D*>(h_data_1GeV->Rebin(nMassBins, "",
                                                         massBoundaries));

    double N_data = getN(h_data, minX, maxX);

    /*
    double parameters[4][5] = {
        // p0, p1, p2, p3, p4
        { // PF RECO
            {3.16766240835e-07, 0.0, 6.27895318762, 0.0, 0.0}, // f2
            {1.51335086414e-05, 9.05904632837, 5.02571873249, 0.0, 0.0}, // f3
            {3.46397990934e-06, 7.34961109006, 5.96523987147, 0.163614755034,
             0.0}, // f4
            {2.38002753645e-05, 8.87638444646, 4.18327834544, -0.457885218595,
             -0.0784871491669} // f5
        },
        { // Calo Scouting
            {3.3110576299e-06, 0.0, 5.40474642534, 0.0, 0.0}, // f2
            {8.58951893126e-05, 14.3144967628, 4.57510842619, 0.0, 0.0}, // f3
            {9.49994522633e-07, 5.84090837867, 6.83225395106, 0.299718539449,
             0.0}, // f4
            {9.52739808228e-07, 5.84779556315, 6.83125819802, 0.299773743681,
             3.51085801569e-05} // f5
        }
    };
    */
    //parameters for Old F4
    double parameters[3][5] = {
      // p0, p1, p2, p3, p4
      {2.2349e+05, -1.9690e+00, -2.8513e+01, 0.0, 0.0}, // f3
      {2.2349e+05, -9.1968e-01, -4.9389e+01, 4.6074e+01, 0.0}, // f4
      {2.2349e+05, 4.9986e-02, -7.4393e+01, 1.4009e+02, -1.3337e+02} // f5
    };

    //    double RSS2 = RSS(h_data, parameters[dataset][0][0],
    //                 parameters[dataset][0][1], parameters[dataset][0][2],
    //                 parameters[dataset][0][3], parameters[dataset][0][4],
    //                  minX, maxX, lumi);
    
    double RSS3 = RSS(h_data, parameters[1][0],
                      parameters[1][1], parameters[1][2],
                      parameters[1][3], parameters[1][4],
                      minX, maxX, lumi);
    double RSS4 = RSS(h_data, parameters[2][0],
                      parameters[2][1], parameters[2][2],
                      parameters[2][3], parameters[2][4],
                      minX, maxX, lumi);
    double RSS5 = RSS(h_data, parameters[3][0],
                      parameters[3][1], parameters[3][2],
                      parameters[3][3], parameters[3][4],
                      minX, maxX, lumi);
    
    cout << endl;
    //    cout << "RSS2: " << RSS2 << endl;
    cout << "RSS3: " << RSS3 << endl;
    cout << "RSS4: " << RSS4 << endl;
    cout << "RSS5: " << RSS5 << endl;

    //    double F32 = (RSS2 - RSS3)/(RSS3/(N_data - 3.0));
    double F43 = (RSS3 - RSS4)/(RSS4/(N_data - 4.0));
    double F54 = (RSS4 - RSS5)/(RSS5/(N_data - 5.0));

    cout << endl;
    //    cout << "F32: " << F32 << endl;
    cout << "F43: " << F43 << endl;
    cout << "F54: " << F54 << endl;

    //    double CL32 = 1.0 - TMath::FDistI(F32, 1.0, N_data - 3.0);
    double CL43 = 1.0 - TMath::FDistI(F43, 1.0, N_data - 4.0);
    double CL54 = 1.0 - TMath::FDistI(F54, 1.0, N_data - 5.0);

    cout << endl;
    //    cout << "CL32: " << CL32 << endl;
    cout << "CL43: " << CL43 << endl;
    cout << "CL54: " << CL54 << endl;

    return;
}

double RSS(TH1D *hist, double p0, double p1, double p2, double p3, double p4,
           double minX, double maxX, double lumi)
{ 
  //    TF1 *func = new TF1("func",
  //                        "[0]*TMath::Power(1.0 - x/13000.0, [1])/TMath::Power(x/13000.0, [2] + [3]*TMath::Log(x/13000) + [4]*TMath::Log(x/13000)*TMath::Log(x/13000))", minX, maxX);
  //Old F4
  TF1 *func = new TF1("func",
                        "[0]*TMath::Power(x/13000.0, [1])*TMath::Exp(1+[2]*TMath::Power(x/13000.0, 1)+[3]*TMath::Power(x/13000.0, 2)+[4]*TMath::Power(x/13000.0, 3))", minX, maxX);
  //New F4
  //  TF1 *func = new TF1("func",
  //                    "[0]*TMath::Exp(1+[1]*TMath::Power(x/13000.0, 1)+[2]*TMath::Power(x/13000.0, 2)+[3]*TMath::Power(x/13000.0, 3)+[4]*TMath::Power(x/13000.0, 4))", minX, maxX);

    func->SetParameter(0, p0);
    func->SetParameter(1, p1);
    func->SetParameter(2, p2);
    func->SetParameter(3, p3);
    func->SetParameter(4, p4);

    double rss = 0.0;
    for (int i=1; i<=hist->GetXaxis()->GetNbins(); ++i) {
        double value = hist->GetBinContent(i);
        double lowX = hist->GetXaxis()->GetBinLowEdge(i);
        double upX = hist->GetXaxis()->GetBinUpEdge(i);
        if (lowX < minX || upX > maxX)
            continue;
        if (value == 0.0)
            continue;
        //double fit = func->Integral(lowX, upX)*lumi;
        double fit = func->Integral(lowX, upX);
        rss += (value - fit)*(value - fit);
    }

    return rss;
}

double getN(TH1D *hist, double minX, double maxX)
{
    int N = 0;
    int zeros = 0;
    for (int i=1; i<=hist->GetXaxis()->GetNbins(); ++i) {
        if (hist->GetXaxis()->GetBinLowEdge(i) < minX
            || hist->GetXaxis()->GetBinUpEdge(i) > maxX)
            continue;

        if (hist->GetBinContent(i) > 0.0)
            ++N;
        else
            ++zeros;
    }

    cout << endl;
    cout << N << " data points between " << minX << " and " << maxX
         << " GeV. " << zeros << " bins with zero in that range." << endl;

    return static_cast<double>(N);
}
