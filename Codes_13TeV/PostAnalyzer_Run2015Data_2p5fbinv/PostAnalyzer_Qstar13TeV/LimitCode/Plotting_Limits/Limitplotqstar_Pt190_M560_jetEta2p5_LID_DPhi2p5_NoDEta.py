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


BR = 1.0


masses_gev = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0])
masses_tev = array('d', [ 0.001 * float(b) for b in masses_gev ])

xs_obs_limits  = array('d', [0.0619262, 0.0747435, 0.0705053, 0.023007, 0.0173682, 0.0252859, 0.0150554, 0.0104444, 0.00593128, 0.00502294, 0.00456198, 0.00410041, 0.00361522, 0.00342761, 0.00302666, 0.00244801, 0.00193925, 0.00166212, 0.0015471, 0.00145717, 0.00146663, 0.00140875, 0.00138473])
xs_obs_limits_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits ])

xs_exp_limits = array('d', [0.0907457, 0.067453, 0.0450353, 0.0319344, 0.0238541, 0.0169504, 0.0129428, 0.00969425, 0.00761593, 0.00604852, 0.00466165, 0.00397718, 0.00333839, 0.00268328, 0.00230749, 0.0020925, 0.00192678, 0.00173963, 0.00159242, 0.00150717, 0.00149152, 0.00143949, 0.00139619])
xs_exp_limits_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits ])

masses_exp_gev = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0, 5400.0, 5200.0, 5000.0, 4800.0, 4600.0, 4400.0, 4200.0, 4000.0, 3800.0, 3600.0, 3400.0, 3200.0, 3000.0, 2800.0, 2600.0, 2400.0, 2200.0, 2000.0, 1800.0, 1600.0, 1400.0, 1200.0, 1000.0])
masses_exp_tev = array('d', [ 0.001 * float(b) for b in masses_exp_gev ])

xs_exp_limits_1sigma  = array('d', [0.0663396, 0.0440913, 0.03332, 0.0230641, 0.0165411, 0.0127336, 0.00915343, 0.00668142, 0.00526271, 0.00421082, 0.0033644, 0.00287679, 0.00242496, 0.00209394, 0.00181965, 0.00164835, 0.0015099, 0.00142974, 0.00137081, 0.00133732, 0.00134133, 0.00131532, 0.00129687, 0.00165772, 0.00175408, 0.00187442, 0.00192823, 0.00200858, 0.00216831, 0.00240765, 0.00279827, 0.00325373, 0.00396049, 0.00455003, 0.00561336, 0.00685461, 0.00851328, 0.0109248, 0.0136372, 0.0176993, 0.024301, 0.0333421, 0.0464012, 0.061782, 0.0884147, 0.123])
xs_exp_limits_1sigma_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits_1sigma ])

xs_exp_limits_2sigma = array('d', [0.0479359, 0.0339254, 0.0262679, 0.0169141, 0.012425, 0.00911116, 0.00689107, 0.00494983, 0.00405135, 0.00337902, 0.00278569, 0.00207617, 0.0017685, 0.00170219, 0.0015245, 0.00145112, 0.00137045, 0.00132927, 0.00130608, 0.00127666, 0.00128414, 0.00127601, 0.00126244, 0.00212793, 0.00228474, 0.0026295, 0.00273521, 0.0031693, 0.00336921, 0.00344483, 0.00379696, 0.00425057, 0.00496627, 0.00607073, 0.00750095, 0.00957491, 0.0123677, 0.0145106, 0.0203232, 0.0242716, 0.0340696, 0.0468813, 0.0618612, 0.0788479, 0.114697, 0.179771])
xs_exp_limits_2sigma_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits_2sigma ])


masses_qstar      = array('d', [ 1.0,       1.1,       1.2,       1.3,       1.4,       1.5,       1.6,       1.7,       1.8,       1.9,
                                 2.0,       2.1,       2.2,       2.3,       2.4,       2.5,       2.6,       2.7,       2.8,       2.9,
				 3.0,       3.1,       3.2,       3.3,       3.4,       3.5,       3.6,       3.7,       3.8,       3.9,
				 4.0,       4.1,       4.2,       4.3,       4.4,       4.5,       4.6,       4.7,       4.8,       4.9,
				 5.0])

#Eff for no deta cut
eff_qstar_LID    = array('d',  [ 0.53586, 0.54994, 0.56401,  0.57809,  0.59217,  0.60624,  0.62032,  0.63440,  0.64848,  0.66255,  
                                 0.67663, 0.68029, 0.68395,  0.68761,  0.69127,  0.69493,  0.69859,  0.70225,  0.70591,  0.70957, 
                                 0.71323, 0.71510, 0.71698,  0.71885,  0.72073,  0.72260,  0.72447,  0.72635,  0.72822,  0.73010,
                                 0.73197, 0.73250, 0.73302,  0.73355,  0.73407,  0.73460,  0.73512,  0.73565,  0.73617,  0.73670, 
                                 0.73722])

xs_qstar_f1p0     = array('d', [ 1.635e+01, 1.068e+01, 7.036e+00, 4.827e+00, 3.363e+00, 2.373e+00, 1.709e+00, 1.243e+00, 9.292e-01, 6.885e-01,
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


xsEff_f1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_LID)])
xsEff_f1p0_fb = array('d', [1000.0 * float(b) for b in xsEff_f1p0])

xsEff_f0p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p5, eff_qstar_LID)])
xsEff_f0p5_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p5])

xsEff_f0p1 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p1, eff_qstar_LID)])
xsEff_f0p1_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p1])


graph_exp_2sigma = TGraph(len(masses_exp_tev),masses_exp_tev,xs_exp_limits_2sigma_fb)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times B #times A #times #epsilon [fb]")
graph_exp_2sigma.GetYaxis().SetRangeUser(5e-1,300)
graph_exp_2sigma.GetXaxis().SetNdivisions(510)
graph_exp_2sigma.GetXaxis().SetLimits(0.6,6)

graph_exp_2sigma.GetYaxis().CenterTitle()
graph_exp_2sigma.GetYaxis().SetLabelSize(0.04)
graph_exp_2sigma.GetYaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetXaxis().CenterTitle()
graph_exp_2sigma.GetXaxis().SetLabelSize(0.04)
graph_exp_2sigma.GetXaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetXaxis().CenterTitle()

graph_exp_1sigma = TGraph(len(masses_exp_tev),masses_exp_tev,xs_exp_limits_1sigma_fb)
graph_exp_1sigma.SetFillColor(kGreen+1)

graph_exp = TGraph(len(masses_tev),masses_tev,xs_exp_limits_fb)
#graph_exp.SetMarkerStyle(24)
graph_exp.SetLineWidth(2)
graph_exp.SetLineStyle(2)
graph_exp.SetLineColor(4)

graph_obs = TGraph(len(masses_tev),masses_tev,xs_obs_limits_fb)
graph_obs.SetMarkerStyle(20)
graph_obs.SetLineWidth(2)
#graph_obs.SetLineStyle(1)
graph_obs.SetLineColor(1)

graph_qstar = TGraph(len(masses_qstar),masses_qstar,xsEff_f1p0_fb)
graph_qstar.SetLineWidth(2)                
graph_qstar.SetLineColor(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
graph_qstar.Draw("L")

##legend = TLegend(.55,.69,.85,.92)  ## for pas twiki
legend = TLegend(.58,.52,.89,.72)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetHeader('95% CL upper limits')
legend.AddEntry(graph_obs,"Observed limit","lp")
legend.AddEntry(graph_exp,"Expected limit","lp")
legend.AddEntry(graph_exp_1sigma,"Expected limit #pm 1#sigma","f")
legend.AddEntry(graph_exp_2sigma,"Expected limit #pm 2#sigma","f")
legend.Draw()

legend1 = TLegend(.16,.18,.47,.24)
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetFillStyle(0)
legend1.SetTextFont(42)
legend1.SetTextSize(0.035)
legend1.AddEntry(graph_qstar,"Excited quark (f = 1.0)","l")
legend1.Draw()

lumiTextSize = 0.6
lumiTextOffset = 0.2
lumi = TLatex()
lumi.SetNDC()
lumi.SetTextAngle(0)
lumi.SetTextColor(kBlack)
lumi.SetTextFont(42)
lumi.SetTextAlign(31)
lumi.SetTextSize(lumiTextSize*t_m)
lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "2.5 fb^{-1} (13 TeV)")

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
cms.DrawLatex(0.9, 0.9, "CMS")

### if extra text (Unpublished or Preliminary)
cms.SetTextFont(extraTextFont);
cms.SetTextAlign(33);
cms.SetTextSize(extraOverCmsTextSize*cmsTextSize*t_m);
##cms.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t_m, extraText);
cms.DrawLatex(0.9, 0.9 - relExtraDY*cmsTextSize*t_m, "Preliminary");

l1 = TLatex()
l1.SetNDC()
l1.SetTextAlign(13)
l1.SetTextFont(62)
l1.SetTextSize(0.04)
l1.SetTextFont(42)
l1.DrawLatex(0.785,0.805, "q*#rightarrow q#gamma") ##0.21,0.26


gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('Limit_f1p0_ObseExp_xsBRAccEff_Pt190_M560_jetEta2p5_LID_DPhi2p5_NoDEta.pdf')
#c.SaveAs('Limit_f1p0_ObseExp_xsBRAccEff_Pt190jetEta2p4dphi2p5dEta2p0M560LID.png')


