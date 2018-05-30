#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.045, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.045, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
t_m = 0.06  ##top margin
b_m = 0.11   ##botton margin
l_m = 0.115  ##left margin
r_m = 0.04  ##right margin
gStyle.SetPadTopMargin(t_m)
gStyle.SetPadBottomMargin(b_m)
gStyle.SetPadLeftMargin(l_m)
gStyle.SetPadRightMargin(r_m)
gROOT.ForceStyle()

masses_gev = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0, 5600.0, 5800.0, 6000.0, 6200.0])

masses_tev = array('d', [ 0.001 * float(b) for b in masses_gev ])

xs_obs_limits_qstar_asymt = array('d', [0.0229369, 0.0111276, 0.0065780, 0.0041319, 0.0030301, 0.0034745, 0.0023898, 0.0019194, 0.0006016, 0.0008710, 0.0010699, 0.0008776, 0.0003817, 0.0002301, 0.0001880, 0.0001314, 0.0001468, 0.0002210, 0.0002349, 0.0001977, 0.0001340, 0.0001061, 9.441e-05, 8.925e-05, 8.469e-05, 7.939e-05, 7.935e-05])

xs_obs_limits_qstar_asymt_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_qstar_asymt ])

xs_obs_limits_qstar_hybrid = array('d', [0.024762056, 0.011707512, 0.006955963, 0.004236571, 0.003075956, 0.003327733, 0.002448240, 0.001942139, 0.000648622, 0.000881947, 0.001060457, 0.000879131, 0.000336838, 0.000289556, 0.000189539, 0.000163028, 0.000135616, 0.000227598, 0.000243655, 0.000209112, 0.000144294, 0.000113245, 0.000117031, 0.000107715, 0.000107786, 0.000107359, 0.000104589])

xs_obs_limits_qstar_hybrid_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_qstar_hybrid ])

masses_qstar =  array('d', [ 1.0,     1.1,      1.2,      1.3,      1.4,      
                             1.5,     1.6,      1.7,      1.8,      1.9,      
                             2.0,     2.1,      2.2,      2.3,      2.4,      
                             2.5,     2.6,      2.7,      2.8,      2.9,      
                             3.0,     3.1,      3.2,      3.3,      3.4,      
                             3.5,     3.6,      3.7,      3.8,      3.9,      
                             4.0,     4.1,      4.2,      4.3,      4.4,      
                             4.5,     4.6,      4.7,      4.8,      4.9,      
                             5.0,     5.1,      5.2,      5.3,      5.4, 
                             5.5,     5.6,      5.7,      5.8,      5.9,    6.0])

xs_qstar_f1p0 =  array('d', [1.635e+01, 1.068e+01, 7.036e+00, 4.827e+00, 3.363e+00, 2.373e+00, 1.709e+00, 1.243e+00, 9.292e-01, 6.885e-01, 5.244e-01, 3.951e-01, 3.010e-01, 2.322e-01, 1.805e-01, 1.402e-01, 1.099e-01, 8.650e-02, 6.787e-02, 5.372e-02, 4.273e-02, 3.391e-02, 2.720e-02, 2.186e-02, 1.744e-02, 1.417e-02, 1.126e-02, 9.062e-03, 7.276e-03, 5.911e-03, 4.814e-03, 3.870e-03, 3.156e-03, 2.554e-03, 2.057e-03, 1.656e-03, 1.354e-03, 1.089e-03, 8.813e-04, 7.214e-04, 5.836e-04, 4.734e-04, 3.807e-04, 3.108e-04, 2.517e-04, 2.051e-04, 1.650e-04, 1.339e-04, 1.072e-04, 8.685e-05, 7.085e-05])

eff_qstar_latest    = array('d', [ 0.397222, 0.402547, 0.407871, 0.413196, 0.41852, 0.423845, 0.429169, 0.434494, 0.439818, 0.445143, 0.450467, 0.45171, 0.452952, 0.454195, 0.455437, 0.45668, 0.457922, 0.459165, 0.460407, 0.46165, 0.462892, 0.462899, 0.462905, 0.462911, 0.462918, 0.462924, 0.462931, 0.462937, 0.462943, 0.46295, 0.462956, 0.4632, 0.463444, 0.463687, 0.463931, 0.464175, 0.464419, 0.464662, 0.464906, 0.46515, 0.465394, 0.464291, 0.463189, 0.462086, 0.460984, 0.459881, 0.458779, 0.45767, 0.456574, 0.455471, 0.454369])

xsEff_f1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_latest)])
xsEff_f1p0_fb = array('d', [1000.0 * float(b) for b in xsEff_f1p0])

graph_obs_asymt = TGraph(len(masses_tev),masses_tev,xs_obs_limits_qstar_asymt_fb)
#graph_obs_aysmt.SetFillColor(kYellow)
graph_obs_asymt.GetXaxis().SetTitle("Resonance Mass [TeV]")
graph_obs_asymt.GetYaxis().SetTitle("#sigma #times B #times A #times #epsilon [fb]")
graph_obs_asymt.GetYaxis().SetRangeUser(1e-2,1e+2)
graph_obs_asymt.GetXaxis().SetNdivisions(510)
graph_obs_asymt.GetXaxis().SetLimits(0.6,6.5)
graph_obs_asymt.SetMarkerStyle(20)
graph_obs_asymt.SetMarkerColor(2)
graph_obs_asymt.SetLineWidth(2)
graph_obs_asymt.SetLineStyle(1)
graph_obs_asymt.SetLineColor(2)

graph_obs_asymt.GetYaxis().CenterTitle()
graph_obs_asymt.GetYaxis().SetLabelSize(0.04)
graph_obs_asymt.GetYaxis().SetTitleOffset(1.1)
graph_obs_asymt.GetXaxis().CenterTitle()
graph_obs_asymt.GetXaxis().SetLabelSize(0.04)
graph_obs_asymt.GetXaxis().SetTitleOffset(1.1)
graph_obs_asymt.GetXaxis().CenterTitle()

graph_obs_hybrid =TGraph(len(masses_tev),masses_tev,xs_obs_limits_qstar_hybrid_fb)
graph_obs_hybrid.SetMarkerStyle(24)
graph_obs_hybrid.SetMarkerColor(1)
graph_obs_hybrid.SetLineWidth(2)
graph_obs_hybrid.SetLineStyle(1)
graph_obs_hybrid.SetLineColor(1)

graph_qstar = TGraph(len(masses_qstar),masses_qstar,xsEff_f1p0_fb)
graph_qstar.SetLineWidth(2)                
graph_qstar.SetLineColor(12)
graph_qstar.SetLineStyle(8)

c = TCanvas("c", "",800,800)
c.cd()

graph_obs_asymt.Draw("ALP")
graph_obs_hybrid.Draw("LP")
graph_qstar.Draw("L")

##legend = TLegend(.55,.69,.85,.92)  ## for pas twiki
legend = TLegend(.54,.54,.76,.75)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.031)
legend.SetHeader('95% CL upper limits')
legend.AddEntry(graph_obs_asymt,"Obs. lmt (Asymptotic)","lp")
legend.AddEntry(graph_obs_hybrid,"Obs. lmt (Full CLs)","lp")
legend.AddEntry(graph_qstar,"Excited quark (f = 1.0)","l")
legend.Draw()

lumiTextSize = 0.6
lumiTextOffset = 0.2
lumi = TLatex()
lumi.SetNDC()
lumi.SetTextAngle(0)
lumi.SetTextColor(kBlack)
lumi.SetTextFont(42)
lumi.SetTextAlign(31)
lumi.SetTextSize(lumiTextSize*t_m)
lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "35.9 fb^{-1} (13 TeV)")

cmsTextFont = 61
cmsTextSize = 0.75
extraTextFont = 52
extraOverCmsTextSize = 0.76
relExtraDY = 1.2
##posX_ = l_m     + 0.045 * (1 - l_m - r_m)  ##Top left
##posX_ = l_m     + 0.5   * (1 - l_m - r_m)   ## Centre
posX_ = 1 - r_m + 0.045 * (1 - l_m - r_m)  ## Top right
posY_ = 1 - t_m - 0.035 * (1 - t_m - b_m)

cms =  TLatex()
cms.SetNDC()
cms.SetTextFont(cmsTextFont)
cms.SetTextSize(cmsTextSize * t_m)
cms.SetTextAlign(33)  ### 11-top left;  21-top centre;  31-top right
##cms.DrawLatex(posX_, posY_, "CMS")
cms.DrawLatex(0.9, 0.92, "CMS")

### if extra text (Unpublished or Preliminary)
cms.SetTextFont(extraTextFont);
cms.SetTextAlign(33);
cms.SetTextSize(extraOverCmsTextSize*cmsTextSize*t_m);
##cms.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t_m, extraText);
#cms.DrawLatex(0.9, 0.9 - relExtraDY*cmsTextSize*t_m, "Preliminary");

l1 = TLatex()
l1.SetNDC()
l1.SetTextAlign(13)
l1.SetTextFont(62)
l1.SetTextSize(0.04)
l1.SetTextFont(42)
l1.DrawLatex(0.78,0.85, "q*#rightarrow q#gamma") ##0.21,0.26


gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('LimitComp_qstar_asymptotic_vs_hybrid.pdf')
c.SaveAs('LimitComp_qstar_asymptotic_vs_hybrid.png')


