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

xs_obs_limits = array('d', [0.037955, 0.00414443, 0.00370767, 0.00209104, 0.00156958, 0.00166186])

xs_exp_limits = array('d', [0.0717511, 0.0810783, 0.00417976, 0.00225336, 0.00181125, 0.00185173])

masses_exp = array('d', [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 7000.0, 7000.0, 5000.0, 4000.0, 3000.0, 2000.0, 1000.0])

xs_exp_limits_1sigma  = array('d', [0.0470744, 0.0570128, 0.00294089, 0.00174258, 0.00152267, 0.00163709, 0.00245075, 0.00247928,  0.00302866, 0.00559832, 0.109487, 0.105772 ])

xs_exp_limits_2sigma = array('d', [0.00403188, 0.041052, 0.00241373, 0.00159855, 0.00146452, 0.00156238 , 0.00347851, 0.0036724, 0.004026,  0.0075283, 0.161284, 0.136349 ])


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


