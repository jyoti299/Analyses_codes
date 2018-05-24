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



masses = array('d', [700.0, 1000.0, 1200.0, 1500.0, 1700.0, 2000.0, 2500.0, 3000.0, 4000.0])

xs_exp_limits = array('d', [0.0200381, 0.00936334, 0.00594011, 0.00356639, 0.00247994, 0.00154272, 0.000762209, 0.000546239, 0.000564629])

xs_exp_limits = array('d', [0.0200381, 0.00936334, 0.00594011, 0.00356639, 0.00247994, 0.00154272, 0.000762209, 0.000546239, 0.000564629])

xs_exp_limits = array('d', [0.0200381, 0.00936334, 0.00594011, 0.00356639, 0.00247994, 0.00154272, 0.000762209, 0.000546239, 0.000564629])


xs_exp_limits_1 = copy.copy()
xs_exp_limits_2 = copy.copy()
xs_exp_limits_3 = copy.copy()



masses_bstar  = array('d', [0.7,   1.0,   1.2,   1.5,    1.7,     2.0,   2.5])
xs_bstar_f1p0 = array('d', [2.096E-1, 2.291E-2, 6.52E-3, 1.214e-3, 4.321e-4,  1.013e-4,  1.072e-5 ])


eff_bstar_f1p0 = array('d', [0.33293, 0.32294, 0.29971, 0.267201,  0.250887,  0.228362,  0.198002])
eff_bstar_f1p0 = array('d', [0.33293, 0.32294, 0.29971, 0.267201,  0.250887,  0.228362,  0.198002])
eff_bstar_f1p0 = array('d', [0.33293, 0.32294, 0.29971, 0.267201,  0.250887,  0.228362,  0.198002])

eff_bstar_f1p0_1 = copy.copy()
eff_bstar_f1p0_2 = copy.copy()
eff_bstar_f1p0_3 = copy.copy()



xstimesEff = array('d', [float(b) * float(m) for b,m in zip(xs_bstar_f1p0, eff_bstar_f1p0)])
xstimesEff = array('d', [float(b) * float(m) for b,m in zip(xs_bstar_f1p0, eff_bstar_f1p0)])
xstimesEff = array('d', [float(b) * float(m) for b,m in zip(xs_bstar_f1p0, eff_bstar_f1p0)])

#print xstimesEff 

result = array('d',[0.7,1.0,1.2,1.5,1.7,2.0,2.5,3.0,4.0])

result_2sigma= array('d',[0.7,1.0,1.2,1.5,1.7,2.0,2.5,3.0,4.0,4.0,3.0,2.5,2.0,1.7,1.5,1.2,1.0,0.7])

graph_exp_2sigma = TGraph(len(masses_exp),result_2sigma,xs_exp_limits_2sigma)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetXaxis().SetTitle("b* Mass [TeV]")
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times A #times #epsilon [pb]")
#graph_exp_2sigma.GetYaxis().SetTitle("#sigma [pb]")
graph_exp_2sigma.GetYaxis().SetRangeUser(1e-05,0.3)
graph_exp_2sigma.GetXaxis().SetNdivisions(510)
graph_exp_2sigma.GetXaxis().SetLimits(0.3,4.3)

graph_exp_2sigma.GetYaxis().CenterTitle()
graph_exp_2sigma.GetYaxis().SetLabelSize(0.05)
graph_exp_2sigma.GetYaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetXaxis().CenterTitle()
graph_exp_2sigma.GetXaxis().SetLabelSize(0.05)
graph_exp_2sigma.GetXaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetXaxis().CenterTitle()

graph_exp_1sigma = TGraph(len(masses_exp),result_2sigma,xs_exp_limits_1sigma)
graph_exp_1sigma.SetFillColor(kGreen+1)

graph_exp = TGraph(len(masses),result,xs_exp_limits)
#graph_exp.SetMarkerStyle(24)
graph_exp.SetLineWidth(2)
graph_exp.SetLineStyle(2)
graph_exp.SetLineColor(4)

graph_obs = TGraph(len(masses),result,xs_obs_limits)
graph_obs.SetMarkerStyle(20)
graph_obs.SetLineWidth(2)
#graph_obs.SetLineStyle(1)
graph_obs.SetLineColor(1)

graph_bstar = TGraph(len(masses_bstar),masses_bstar,xstimesEff)
#graph_bstar = TGraph(len(masses_bstar),masses_bstar,xs_bstar_f1p0)
graph_bstar.SetLineWidth(2)                
graph_bstar.SetLineColor(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
graph_bstar.Draw("L")

legend = TLegend(.55,.69,.85,.92)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
legend.SetHeader('95% CL Upper Limits')
legend.AddEntry(graph_obs,"Observed Limit","lp")
legend.AddEntry(graph_exp,"Expected Limit","lp")
legend.AddEntry(graph_exp_1sigma,"Expected Limit #pm 1#sigma","f")
legend.AddEntry(graph_exp_2sigma,"Expected Limit #pm 2#sigma","f")
legend.AddEntry(graph_bstar,"Excited b-Quark ( f = 1.0 )","l")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.035)
l1.DrawLatex(0.215,0.44, "CMS Preliminary")
#l1.DrawLatex(0.215,0.37, "CMS (unpublished)")
l1.DrawLatex(0.20,0.38, "#sqrt{s} = 8 TeV")
l1.DrawLatex(0.20,0.32, "#int Ldt = 19.7 fb^{-1}")
l1.DrawLatex(0.20,0.25, "b* #rightarrow b#gamma")

gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xsAccEff_Limits_an.pdf')
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xsAccEff_Limits_an.eps')
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xsAccEff_Limits_an.png')

#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.eps')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.png')


