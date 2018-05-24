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

xs_exp_limits_DPhi_1p5 = array('d', [0.0824724, 0.0518581, 0.0368938, 0.0263986, 0.0175979, 0.0119491, 0.00907025, 0.00736071, 0.00586477, 0.00461258, 0.00373162, 0.00318134, 0.0027413, 0.00250953, 0.00226311, 0.00201437, 0.00186699, 0.00162731, 0.00153964, 0.00146261, 0.00147057, 0.00141497, 0.00140955])

xs_exp_limits_DPhi_2p0 = array('d', [0.0806345, 0.0508739, 0.0352324, 0.0241545, 0.016899, 0.0117564, 0.00960138, 0.00749883, 0.00568849, 0.00467033, 0.00381436, 0.0032936, 0.00283435, 0.00243226, 0.00227701, 0.00208425, 0.0019525, 0.0017481, 0.0016516, 0.00156941, 0.00158025, 0.00149348, 0.0014427])

#these numbers for DEta=2p0 from varun (Read README)
xs_exp_limits_DPhi_2p5 = array('d', [0.0753898, 0.0481741, 0.0331191, 0.0228193, 0.0168115, 0.0110654, 0.00836895, 0.00696791, 0.00545557, 0.00469033, 0.00389036, 0.00333766, 0.00296295, 0.00251753, 0.00229563, 0.00213392, 0.00202398, 0.00186637, 0.00176224, 0.00162506, 0.00160681, 0.00155983, 0.00151451])

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

eff_qstar_DPhi_1p5  = array('d', [0.52772,   0.53663,   0.54556,   0.55447,    0.56334,    0.57217,   0.58104,   0.58902,  0.59872,   0.60812,  
                                  0.61753,   0.61887,   0.62021,   0.62155,    0.62289,    0.62423,   0.62557,   0.62692,  0.62826,   0.62960, 
                                  0.63098,   0.63144,   0.63190,   0.63236,    0.63282,    0.63328,   0.63374,   0.63420,  0.63466,   0.63512,
                                  0.63558,   0.63532,   0.63506,   0.63480,    0.63454,    0.63428,   0.63401,   0.63374,  0.63348,   0.63322, 
                                  0.63300])

eff_qstar_DPhi_2p0   = array('d', [0.52143,   0.53068,   0.53993,   0.54918,    0.55843,    0.56768,   0.57693,   0.58618,  0.59543,   0.60468,    
                                   0.61392,   0.61546,   0.61700,   0.61854,    0.62008,    0.62162,   0.62316,   0.62470,  0.62624,   0.62778,
                                   0.62930,   0.62981,   0.63032,   0.63083,    0.63134,    0.63185,   0.63236,   0.63287,  0.63338,   0.63389, 
                                   0.63438,   0.63415,   0.63392,   0.63369,    0.63346,    0.63323,   0.63300,   0.63277,  0.63254,   0.63231,  
                                   0.63205])

#These eff numbers for DEta=2p0 from varun (Read README)
eff_qstar_DPhi_2p5   = array('d',  [ 0.488302,  0.49915,   0.51,      0.52,     0.53,      0.54,      0.55,      0.56,      0.57,      0.5793, 
                                     0.58861,   0.59115,   0.5937,    0.59634,   0.599,     0.60195,   0.6049,    0.60745,   0.61,      0.5993,
                                     0.615333,  0.61636,   0.6174,    0.61845,   0.6195,    0.62065,   0.6218,    0.62295,   0.6241,    0.62511,
                                     0.626118,  0.62584,   0.625565,  0.625537,  0.62551,   0.62538,   0.62525,   0.625065,  0.62488,   0.62474,
                                     0.624607])


xstimesEff_DPhi_1p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DPhi_1p5)])
xstimesEff_DPhi_2p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DPhi_2p0)])
xstimesEff_DPhi_2p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DPhi_2p5)])

#print xstimesEff 

graph_exp_DPhi_1p5 = TGraph(len(masses_tev), masses_tev, xs_exp_limits_DPhi_1p5)
graph_exp_DPhi_1p5.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_DPhi_1p5.GetYaxis().SetTitle("#sigma #times BR #times A #times #epsilon [pb]")
#graph_exp_DPhi_1p5.GetYaxis().SetRangeUser(1e-4,1)
graph_exp_DPhi_1p5.GetYaxis().SetRangeUser(5e-4,5e-3)
graph_exp_DPhi_1p5.GetXaxis().SetNdivisions(510)
#graph_exp_DPhi_1p5.GetXaxis().SetLimits(0.6,6.0) 
graph_exp_DPhi_1p5.GetXaxis().SetLimits(4.1,4.3) 

graph_exp_DPhi_1p5.GetYaxis().CenterTitle()
graph_exp_DPhi_1p5.GetYaxis().SetLabelSize(0.05)
graph_exp_DPhi_1p5.GetYaxis().SetTitleOffset(1.1)
graph_exp_DPhi_1p5.GetXaxis().CenterTitle()
graph_exp_DPhi_1p5.GetXaxis().SetLabelSize(0.05)
graph_exp_DPhi_1p5.GetXaxis().SetTitleOffset(1.1)

graph_exp_DPhi_1p5.SetMarkerStyle(24)
graph_exp_DPhi_1p5.SetMarkerColor(3)
graph_exp_DPhi_1p5.SetMarkerSize(0.7)
graph_exp_DPhi_1p5.SetLineWidth(2)
graph_exp_DPhi_1p5.SetLineColor(3)

graph_exp_DPhi_2p0 = TGraph(len(masses_tev), masses_tev, xs_exp_limits_DPhi_2p0)
graph_exp_DPhi_2p0.SetMarkerStyle(27)
graph_exp_DPhi_2p0.SetMarkerColor(4)
graph_exp_DPhi_2p0.SetMarkerSize(0.7)
graph_exp_DPhi_2p0.SetLineWidth(2)
graph_exp_DPhi_2p0.SetLineColor(4)

graph_exp_DPhi_2p5 = TGraph(len(masses_tev), masses_tev, xs_exp_limits_DPhi_2p5)
graph_exp_DPhi_2p5.SetMarkerStyle(26)
graph_exp_DPhi_2p5.SetMarkerColor(2)
graph_exp_DPhi_2p5.SetMarkerSize(0.7)
graph_exp_DPhi_2p5.SetLineWidth(2)
graph_exp_DPhi_2p5.SetLineColor(2)


graph_qstar_DPhi_1p5 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DPhi_1p5)
graph_qstar_DPhi_1p5.SetLineWidth(2)
graph_qstar_DPhi_1p5.SetLineColor(3)
graph_qstar_DPhi_1p5.SetLineStyle(2)

graph_qstar_DPhi_2p0 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DPhi_2p0)
graph_qstar_DPhi_2p0.SetLineWidth(2)
graph_qstar_DPhi_2p0.SetLineColor(4)
graph_qstar_DPhi_2p0.SetLineStyle(2)

graph_qstar_DPhi_2p5 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DPhi_2p5)
graph_qstar_DPhi_2p5.SetLineWidth(2)
graph_qstar_DPhi_2p5.SetLineColor(2)
graph_qstar_DPhi_2p5.SetLineStyle(2)


c = TCanvas("c", "",800,800)
c.cd()

graph_exp_DPhi_1p5.Draw("ALP")
graph_exp_DPhi_2p0.Draw("LP")
graph_exp_DPhi_2p5.Draw("LP")

graph_qstar_DPhi_1p5.Draw("L")
graph_qstar_DPhi_2p0.Draw("L")
graph_qstar_DPhi_2p5.Draw("L")


legend = TLegend(.55,.69,.85,.92)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.025)
legend.SetHeader('95% CL Upper Expected Limits')
legend.AddEntry(graph_exp_DPhi_1p5, "Expected (#Delta #phi > 1.5)","lp")
legend.AddEntry(graph_exp_DPhi_2p0, "Expected (#Delta #phi > 2.0)","lp")
legend.AddEntry(graph_exp_DPhi_2p5, "Expected (#Delta #phi > 2.5)","lp")

legend.AddEntry(graph_qstar_DPhi_1p5, "q* (#Delta #phi > 1.5)", "l")
legend.AddEntry(graph_qstar_DPhi_2p0, "q* (#Delta #phi > 2.0)", "l")
legend.AddEntry(graph_qstar_DPhi_2p5, "q* (#Delta #phi > 2.5)", "l")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.035)
l1.DrawLatex(0.20,0.43, "CMS Preliminary")
#l1.DrawLatex(0.215,0.37, "CMS (unpublished)")
l1.DrawLatex(0.20,0.38, "#sqrt{s} = 13 TeV")
l1.DrawLatex(0.20,0.32, "#int Ldt = 2.5 fb^{-1}")

gPad.RedrawAxis();

#c.SetLogy()
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DPhiOptimization_ExpVsQstar.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DPhiOptimization_ExpVsQstar.png')

c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DPhiOptimization_ExpVsQstar_zoomed.pdf')
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DPhiOptimization_ExpVsQstar_zoomed.png')

#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.eps')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.png')


