#!/usr/bin/env python

import string, re
from ROOT import *
from array import array
import CMS_lumi


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadBottomMargin(0.14)
gROOT.ForceStyle()


masses = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0])

xs_stat = array('d', [0.445766, 0.531756, 0.660085, 0.462665, 0.259771, 0.176882, 0.11975, 0.0883606, 0.0911673, 0.091631, 0.0673097, 0.0576683, 0.0663765, 0.0715058, 0.0707653, 0.0679017, 0.0648869, 0.0605661, 0.0551873, 0.0536497, 0.0525322, 0.0449858, 0.0376862, 0.0386666, 0.0394637, 0.032315, 0.0210878, 0.0132809, 0.00920013, 0.00718502, 0.00635795, 0.00565515, 0.00524019, 0.00555239, 0.00664857, 0.00749758, 0.00812448, 0.00821765, 0.00778863, 0.0072654, 0.00689535, 0.00641743, 0.00594571, 0.00542539, 0.00487017, 0.00432282, 0.00385402, 0.00345907, 0.00325031, 0.0029336, 0.00266552, 0.00247031, 0.0023408, 0.00224978, 0.00218731, 0.0021452, 0.00214538, 0.0021549])
xs_stat_exp = array('d', [0.477582, 0.401964, 0.363716, 0.289744, 0.23842, 0.210805, 0.184404, 0.156756, 0.133975, 0.116447, 0.101665, 0.0899399, 0.0785615, 0.07051, 0.059203, 0.0554069, 0.0494558, 0.0459754, 0.0401039, 0.0355165, 0.031371, 0.0280482, 0.0250272, 0.0237096, 0.0207305, 0.0191532, 0.0166698, 0.0160243, 0.0137829, 0.0132228, 0.0124812, 0.0118259, 0.010205, 0.00943852, 0.0084869, 0.00801166, 0.00735571, 0.00677923, 0.0063039, 0.00570682, 0.00525888, 0.0049462, 0.00463206, 0.00440061, 0.00415374, 0.00370666, 0.00357487, 0.00338426, 0.00335521, 0.00328086, 0.00321736, 0.003065, 0.0029891, 0.00288091, 0.00276867, 0.00267989, 0.0026509, 0.0026118])

xs_sys_all = array('d', [0.895858, 0.965575, 0.900752, 0.697496, 0.423405, 0.26049, 0.181675, 0.132657, 0.125319, 0.119397, 0.0987246, 0.0877885, 0.0911367, 0.0928461, 0.0908348, 0.0858027, 0.0808088, 0.0741853, 0.0681326, 0.0642383, 0.060582, 0.0544358, 0.0476643, 0.0443734, 0.0423859, 0.037657, 0.0311503, 0.0211091, 0.0130474, 0.00900213, 0.00728843, 0.00637582, 0.00597609, 0.00623773, 0.00705963, 0.00787881, 0.0082644, 0.00835833, 0.00805248, 0.00763612, 0.00711483, 0.00660508, 0.00613461, 0.00565583, 0.00512362, 0.00461948, 0.00412781, 0.00375494, 0.00339019, 0.003035, 0.00278894, 0.00257018, 0.00243658, 0.00227797, 0.00219318, 0.00215883, 0.00214122, 0.0021581])
xs_sys_all_exp = array('d', [0.9946315, 0.8051145, 0.627147, 0.475644, 0.374002, 0.304605, 0.2522265, 0.208702, 0.183822, 0.161106, 0.1410845, 0.123595, 0.1101135, 0.0996241, 0.0878194, 0.0777768, 0.0679316, 0.05905805, 0.0536272, 0.0459997, 0.0409238, 0.0364848, 0.0330868, 0.0286179, 0.0252488, 0.02282465, 0.02144215, 0.0180953, 0.0170692, 0.0151538, 0.01402765, 0.01239115, 0.01140105, 0.01064755, 0.00952408, 0.00845977, 0.00801375, 0.007042805, 0.006565405, 0.00586786, 0.00550331, 0.00508579, 0.00470573, 0.00448573, 0.00417389, 0.003885045, 0.0036984, 0.00356998, 0.003465925, 0.00331277, 0.0031561, 0.003022025, 0.00295003, 0.00287352, 0.002794315, 0.002722985, 0.002670735, 0.00265583])


r_all_exp = array('d')
r_all = array('d')

for i in range(0,len(xs_stat)):
  r_all_exp.append( xs_sys_all_exp[i] / xs_stat_exp[i] )
  r_all.append( xs_sys_all[i] / xs_stat[i] )

g_all_exp = TGraph(len(masses),masses,r_all_exp)
g_all_exp.SetMarkerStyle(24)
g_all_exp.SetMarkerColor(kGreen+2)
g_all_exp.SetLineWidth(2)
g_all_exp.SetLineStyle(2)
g_all_exp.SetLineColor(kGreen+2)
g_all_exp.GetXaxis().SetTitle("qq resonance mass [GeV]")
g_all_exp.GetYaxis().SetTitle("Limit ratio")
g_all_exp.GetYaxis().SetTitleOffset(1.1)
g_all_exp.GetYaxis().SetRangeUser(0.5,3.5)
#g_all_exp.GetXaxis().SetNdivisions(1005)

g_all = TGraph(len(masses),masses,r_all)
g_all.SetMarkerStyle(20)
g_all.SetMarkerColor(kBlack)
g_all.SetLineWidth(2)
g_all.SetLineStyle(1)
g_all.SetLineColor(kBlack)

c = TCanvas("c", "",800,800)
c.cd()

g_all_exp.Draw("ALP")
g_all.Draw("LP")

legend = TLegend(.40,.65,.60,.75)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetMargin(0.20)
legend.AddEntry(g_all_exp, "All / Stat. only (expected)","lp")
legend.AddEntry(g_all, "All / Stat. only (observed)","lp")

legend.Draw()

#l1 = TLatex()
#l1.SetTextAlign(12)
#l1.SetTextFont(42)
#l1.SetNDC()
#l1.SetTextSize(0.04)
#l1.DrawLatex(0.18,0.43, "CMS Preliminary")
#l1.DrawLatex(0.18,0.35, "#intLdt = 5 fb^{-1}")
#l1.DrawLatex(0.19,0.30, "#sqrt{s} = 7 TeV")
#l1.DrawLatex(0.18,0.25, "|#eta| < 2.5, |#Delta#eta| < 1.3")
#l1.DrawLatex(0.18,0.20, "Wide Jets")

#draw the lumi text on the canvas
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "2.4 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

CMS_lumi.CMS_lumi(c, iPeriod, iPos)

c.SetGridx()
c.SetGridy()

c.SaveAs('xs_limit_ratio_DijetLimitCode_qq_Run2_13TeV_DATA_2445_invpb.eps')
