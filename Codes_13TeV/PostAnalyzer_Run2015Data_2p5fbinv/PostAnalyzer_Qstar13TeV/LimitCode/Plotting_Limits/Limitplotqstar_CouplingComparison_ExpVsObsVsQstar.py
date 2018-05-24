#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.06, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()


masses = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0])

masses_tev = array('d', [ 0.001 * float(b) for b in masses ])

xs_exp_limits_f1p0 = array('d', [0.0753898, 0.0481741, 0.0331191, 0.0228193, 0.0168115, 0.0110654, 0.00836895, 0.00696791, 0.00545557, 0.00469033, 0.00389036, 0.00333766, 0.00296295, 0.00251753, 0.00229563, 0.00213392, 0.00202398, 0.00186637, 0.00176224, 0.00162506, 0.00160681, 0.00155983, 0.00151451])

xs_exp_limits_f0p5 = array('d', [0.0634588, 0.0397915, 0.0274431, 0.0189898, 0.0142143, 0.00917415, 0.00756611, 0.0060682, 0.00511302, 0.00430562, 0.0033858, 0.00299268, 0.00258906, 0.00236565, 0.00217617, 0.00205313, 0.00186991, 0.00169143, 0.00162074, 0.00153837, 0.0014903, 0.00151204, 0.00141419])

xs_exp_limits_f0p1 = array('d', [0.0600064, 0.0369641, 0.0250958, 0.0175956, 0.013085, 0.00878582, 0.007144, 0.00582527, 0.00466438, 0.00409728, 0.00319358, 0.00283142, 0.00243939, 0.00219535, 0.00210048, 0.00193973, 0.00177434, 0.00162138, 0.00155422, 0.00151111, 0.0014473, 0.0014124, 0.00143429])

xs_obs_limits_f1p0 = array('d', [0.0400255, 0.0335698, 0.0366197, 0.0309246, 0.0323951, 0.0275996, 0.0134593, 0.00852663, 0.00386713, 0.00356322, 0.00333334, 0.002743, 0.00223023, 0.00214148, 0.00196558, 0.00171071, 0.00153965, 0.00146057, 0.00140387, 0.00137167, 0.00137552, 0.00136663, 0.0014239])

xs_obs_limits_f0p5 = array('d', [0.0358287, 0.028908, 0.0299882, 0.0228818, 0.0263926, 0.0242519, 0.0111218, 0.00796278, 0.00340258, 0.00323503, 0.00305676, 0.00259461, 0.00194774, 0.00205362, 0.00191006, 0.00160518, 0.00147694, 0.00140558, 0.00137932, 0.00142995, 0.0014136, 0.00139592, 0.00137765])

xs_obs_limits_f0p1 = array('d', [0.0342172, 0.0286537, 0.0284372, 0.0207648, 0.0250083, 0.0233782, 0.0103408, 0.00774566, 0.00322188, 0.00314437, 0.00295733, 0.00252371, 0.00186357, 0.00204071, 0.00189796, 0.00157957, 0.00145849, 0.001389, 0.00136628, 0.00141632, 0.00139346, 0.00137592, 0.00135984])


masses_qstar        = array('d', [ 1.0,       1.1,       1.2,       1.3,       1.4,       1.5,       1.6,       1.7,       1.8,       1.9,
                                   2.0,       2.1,       2.2,       2.3,       2.4,       2.5,       2.6,       2.7,       2.8,       2.9,
                                   3.0,       3.1,       3.2,       3.3,       3.4,       3.5,       3.6,       3.7,       3.8,       3.9,
                                   4.0,       4.1,       4.2,       4.3,       4.4,       4.5,       4.6,       4.7,       4.8,       4.9,
                                   5.0])

xs_qstar_f1p0       = array('d', [ 1.635e+01, 1.068e+01, 7.036e+00, 4.827e+00, 3.363e+00, 2.373e+00, 1.709e+00, 1.243e+00, 9.292e-01, 6.885e-01,
                                 5.244e-01, 3.951e-01, 3.010e-01, 2.322e-01, 1.805e-01, 1.402e-01, 1.099e-01, 8.650e-02, 6.787e-02, 5.372e-02,
                                 4.273e-02, 3.391e-02, 2.720e-02, 2.186e-02, 1.744e-02, 1.417e-02, 1.126e-02, 9.062e-03, 7.276e-03, 5.911e-03,
                                 4.814e-03, 3.870e-03, 3.156e-03, 2.554e-03, 2.057e-03, 1.656e-03, 1.354e-03, 1.089e-03, 8.813e-04, 7.214e-04,
                                 5.836e-04])

xs_qstar_f0p5     = array('d', [ 4.137e+00, 2.642e+00, 1.768e+00, 1.217e+00, 8.445e-01, 6.012e-01, 4.345e-01, 3.179e-01, 2.342e-01, 1.765e-01,
                                 1.328e-01, 1.005e-01, 7.712e-02, 5.922e-02, 4.583e-02, 3.601e-02, 2.799e-02, 2.206e-02, 1.746e-02, 1.378e-02,
                                  1.096e-02, 8.642e-03, 7.002e-03, 5.531e-03, 4.407e-03, 3.554e-03, 2.860e-03, 2.302e-03, 1.851e-03, 1.488e-03,
                                  1.211e-03, 9.753e-04, 7.847e-04, 6.374e-04, 5.156e-04, 4.187e-04, 3.360e-04, 2.728e-04, 2.189e-04, 1.770e-04,
                                  1.437e-04])
                                  
xs_qstar_f0p1     = array('d',  [ 1.655e-01, 1.057e-01, 7.134e-02, 4.932e-02, 3.421e-02, 2.440e-02, 1.750e-02, 1.284e-02, 9.433e-03, 7.075e-03,
                                 5.298e-03, 4.025e-03, 3.107e-03, 2.374e-03, 1.861e-03, 1.431e-03, 1.130e-03, 8.902e-04, 7.051e-04, 5.527e-04,
                                   4.363e-04, 3.511e-04, 2.784e-04, 2.232e-04, 1.791e-04, 1.435e-04, 1.148e-04, 9.267e-05, 7.459e-05, 6.014e-05,
                                   4.852e-05, 3.902e-05, 3.157e-05, 2.536e-05, 2.058e-05, 1.677e-05, 1.344e-05, 1.087e-05, 8.690e-06, 7.102e-06, 
                                   5.739e-06])

#Taken same for all three couplings as not much difference
eff_qstar_LID    = array('d',  [ 0.488302,  0.49915,   0.51,      0.52,     0.53,      0.54,      0.55,      0.56,      0.57,      0.5793, 
                                 0.58861,   0.59115,   0.5937,    0.59634,   0.599,     0.60195,   0.6049,    0.60745,   0.61,      0.5993,
                                 0.615333,  0.61636,   0.6174,    0.61845,   0.6195,    0.62065,   0.6218,    0.62295,   0.6241,    0.62511,
                                 0.626118,  0.62584,   0.625565,  0.625537,  0.62551,   0.62538,   0.62525,   0.625065,  0.62488,   0.62474,
                                 0.624607])

xstimesEff_f1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_LID)])
xstimesEff_f0p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p5, eff_qstar_LID)])
xstimesEff_f0p1 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p1, eff_qstar_LID)])

#print xstimesEff 

graph_exp_f1p0 = TGraph(len(masses_tev), masses_tev, xs_exp_limits_f1p0)
graph_exp_f1p0.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_f1p0.GetYaxis().SetTitle("#sigma #times BR #times A #times #epsilon [pb]")
graph_exp_f1p0.GetYaxis().SetRangeUser(1e-4,0.3)
#graph_exp_f1p0.GetYaxis().SetRangeUser(5e-4,5e-3)
graph_exp_f1p0.GetXaxis().SetNdivisions(510)
graph_exp_f1p0.GetXaxis().SetLimits(0.6,6.0) 
#graph_exp_f1p0.GetXaxis().SetLimits(4.2,4.4) 

graph_exp_f1p0.GetYaxis().CenterTitle()
graph_exp_f1p0.GetYaxis().SetLabelSize(0.05)
graph_exp_f1p0.GetYaxis().SetTitleOffset(1.1)
graph_exp_f1p0.GetXaxis().CenterTitle()
graph_exp_f1p0.GetXaxis().SetLabelSize(0.05)
graph_exp_f1p0.GetXaxis().SetTitleOffset(1.1)

graph_exp_f1p0.SetMarkerStyle(24)
graph_exp_f1p0.SetMarkerColor(9)
graph_exp_f1p0.SetMarkerSize(0.7)
graph_exp_f1p0.SetLineWidth(1)
graph_exp_f1p0.SetLineColor(9)

graph_exp_f0p5 = TGraph(len(masses_tev), masses_tev, xs_exp_limits_f0p5)
graph_exp_f0p5.SetMarkerStyle(26)
graph_exp_f0p5.SetMarkerColor(8)
graph_exp_f0p5.SetMarkerSize(0.7)
graph_exp_f0p5.SetLineWidth(1)
graph_exp_f0p5.SetLineColor(8)

graph_exp_f0p1 = TGraph(len(masses_tev), masses_tev, xs_exp_limits_f0p1)
graph_exp_f0p1.SetMarkerStyle(27)
graph_exp_f0p1.SetMarkerColor(46)
graph_exp_f0p1.SetMarkerSize(0.7)
graph_exp_f0p1.SetLineWidth(1)
graph_exp_f0p1.SetLineColor(46)

graph_obs_f1p0 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_f1p0)
graph_obs_f1p0.SetMarkerStyle(20)
graph_obs_f1p0.SetMarkerColor(9)
graph_obs_f1p0.SetMarkerSize(0.7)
graph_obs_f1p0.SetLineWidth(2)
graph_obs_f1p0.SetLineColor(9)

graph_obs_f0p5 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_f0p5)
graph_obs_f0p5.SetMarkerStyle(22)
graph_obs_f0p5.SetMarkerColor(8)
graph_obs_f0p5.SetMarkerSize(0.7)
graph_obs_f0p5.SetLineWidth(2)
graph_obs_f0p5.SetLineColor(8)

graph_obs_f0p1 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_f0p1)
graph_obs_f0p1.SetMarkerStyle(33)
graph_obs_f0p1.SetMarkerColor(46)
graph_obs_f0p1.SetMarkerSize(0.7)
graph_obs_f0p1.SetLineWidth(2)
graph_obs_f0p1.SetLineColor(46)


graph_qstar_f1p0 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_f1p0)
graph_qstar_f1p0.SetLineWidth(2)
graph_qstar_f1p0.SetLineColor(9)
graph_qstar_f1p0.SetLineStyle(2)

graph_qstar_f0p5 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_f0p5)
graph_qstar_f0p5.SetLineWidth(2)
graph_qstar_f0p5.SetLineColor(8)
graph_qstar_f0p5.SetLineStyle(2)

graph_qstar_f0p1 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_f0p1)
graph_qstar_f0p1.SetLineWidth(2)
graph_qstar_f0p1.SetLineColor(46)
graph_qstar_f0p1.SetLineStyle(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_exp_f1p0.Draw("ALP")
graph_exp_f0p5.Draw("LP")
graph_exp_f0p1.Draw("LP")
graph_obs_f1p0.Draw("LP")
graph_obs_f0p5.Draw("LP")
graph_obs_f0p1.Draw("LP")

graph_qstar_f1p0.Draw("L")
graph_qstar_f0p5.Draw("L")
graph_qstar_f0p1.Draw("L")

legend = TLegend(.55,.69,.85,.92)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.023)
legend.SetHeader('95% CL Upper Limits')
legend.AddEntry(graph_exp_f1p0, "Expected Limit (f = 1.0)","lp" )
legend.AddEntry(graph_exp_f0p5, "Expected Limit (f = 0.5)","lp")
legend.AddEntry(graph_exp_f0p1, "Expected Limit (f = 0.1)","lp")
legend.AddEntry(graph_obs_f1p0, "Observed Limit (f = 1.0)","lp" )
legend.AddEntry(graph_obs_f0p5, "Observed Limit (f = 0.5)","lp")
legend.AddEntry(graph_obs_f0p1, "Observed Limit (f = 0.1)","lp")

legend.AddEntry(graph_qstar_f1p0, "q* (f = 1.0)", "l")
legend.AddEntry(graph_qstar_f0p5, "q* (f = 0.5)", "l")
legend.AddEntry(graph_qstar_f0p1, "q* (f = 0.1)", "l")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.025)
l1.DrawLatex(0.20,0.40, "CMS Preliminary")
#l1.DrawLatex(0.215,0.37, "CMS (unpublished)")
l1.DrawLatex(0.20,0.36, "#sqrt{s} = 13 TeV")
l1.DrawLatex(0.20,0.32, "#int Ldt = 2.5 fb^{-1}")

gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('ExcitedbQuarksToGbJ_CouplingComparison_ExpVsObsVsQstar.pdf')
c.SaveAs('ExcitedbQuarksToGbJ_CouplingComparison_ExpVsObsVsQstar.png')

#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ExpVsQstar_zoomed.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ExpVsQstar_zoomed.png')



