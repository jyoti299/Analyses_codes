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



masses = array('d', [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 7000.0])

xs_obs_limits = array('d', [0.0371239, 0.023107, 0.00376865, 0.00208556, 0.00158382, 0.00167191])

xs_exp_limits = array('d', [0.0746867, 0.0106247, 0.00399112, 0.00245635, 0.00177422, 0.00184932])

masses_exp = array('d', [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 7000.0, 7000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0])

xs_exp_limits_1sigma  = array('d', [0.0523496, 0.00611, 0.00282453, 0.00182564, 0.00151632, 0.00162367, 0.00236121, 0.00240354, 0.00329348, 0.00590026, 0.015902, 0.100877 ])

xs_exp_limits_2sigma = array('d', [0.00404762, 0.00391847, 0.00214452, 0.001579, 0.00144512, 0.00155052, 0.00327205, 0.00324095, 0.00418183, 0.00731276, 0.0227643, 0.14318 ])


masses_qstar  = array('d', [1.0,   2.0,   3.0,    4.0,     5.0,   7.0])
xs_qstar_f1p0 = array('d', [1.632E1,  5.213E-1,  4.272E-2,  4.8E-3,  5.835E-4, 8.66E-6])
eff_qstar_f1p0 = array('d', [0.4400, 0.5386, 0.5623, 0.5754, 0.5784, 0.5372])

xstimesEff = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_f1p0)])

#print xstimesEff 

result = array('d',[1.0,2.0,3.0,4.0,5.0,7.0])

result_2sigma= array('d',[1.0,2.0,3.0,4.0,5.0,7.0,7.0,5.0,4.0,3.0,2.0,1.0])

graph_exp_2sigma = TGraph(len(masses_exp),result_2sigma,xs_exp_limits_2sigma)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times A #times #epsilon [pb]")
#graph_exp_2sigma.GetYaxis().SetTitle("#sigma [pb]")
graph_exp_2sigma.GetYaxis().SetRangeUser(1e-04,0.7)
graph_exp_2sigma.GetXaxis().SetNdivisions(510)
graph_exp_2sigma.GetXaxis().SetLimits(0.3,7.5)

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

graph_qstar = TGraph(len(masses_qstar),masses_qstar,xstimesEff)
#graph_qstar = TGraph(len(masses_qstar),masses_qstar,xs_qstar_f1p0)
graph_qstar.SetLineWidth(2)                
graph_qstar.SetLineColor(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
graph_qstar.Draw("L")

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
legend.AddEntry(graph_qstar,"Excited Quark ( f = 1.0 )","l")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.035)
l1.DrawLatex(0.215,0.44, "CMS Preliminary")
#l1.DrawLatex(0.215,0.37, "CMS (unpublished)")
l1.DrawLatex(0.20,0.38, "#sqrt{s} = 13 TeV")
l1.DrawLatex(0.20,0.32, "#int Ldt = 2.2 fb^{-1}")
l1.DrawLatex(0.20,0.25, "q* #rightarrow q#gamma")

gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('ExcitedQuarksToGJ_f1p0_ObseExp_xsAccEff_Limits_an.pdf')
c.SaveAs('ExcitedQuarksToGJ_f1p0_ObseExp_xsAccEff_Limits_an.eps')
c.SaveAs('ExcitedQuarksToGJ_f1p0_ObseExp_xsAccEff_Limits_an.png')

#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.eps')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_ObseExp_xs_Limits_an.png')


