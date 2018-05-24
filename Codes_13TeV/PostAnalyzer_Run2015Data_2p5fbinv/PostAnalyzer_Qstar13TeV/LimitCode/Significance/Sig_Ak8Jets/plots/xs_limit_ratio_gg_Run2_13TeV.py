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

xs_stat = array('d', [0.836478, 0.774447, 1.23998, 0.959869, 0.526238, 0.310743, 0.219122, 0.14414, 0.117179, 0.131122, 0.127128, 0.0906086, 0.104278, 0.117512, 0.122472, 0.120322, 0.121221, 0.113741, 0.102283, 0.0975416, 0.0961606, 0.0925603, 0.0792465, 0.0676229, 0.0701066, 0.0671454, 0.0527954, 0.0336437, 0.0224766, 0.0150056, 0.0121905, 0.0103843, 0.00921951, 0.00866013, 0.00924763, 0.0110804, 0.0127599, 0.0139334, 0.0141779, 0.013463, 0.0126392, 0.0119217, 0.0112133, 0.010475, 0.00967835, 0.00879789, 0.00802754, 0.00750105, 0.00693862, 0.00675818, 0.00620987, 0.00612767, 0.00588462, 0.00575895, 0.00572472, 0.00575174, 0.00590138, 0.00613816])
xs_stat_exp = array('d', [0.807865, 0.749822, 0.636492, 0.501165, 0.442055, 0.373683, 0.308796, 0.274763, 0.2461, 0.21873, 0.197548, 0.159459, 0.144599, 0.12325, 0.112418, 0.101445, 0.092047, 0.0772591, 0.0728687, 0.0636055, 0.0546508, 0.0506898, 0.0467854, 0.0444067, 0.0383595, 0.0337996, 0.0299116, 0.0287613, 0.0256873, 0.0245892, 0.0222354, 0.0203452, 0.0190944, 0.0182157, 0.015392, 0.0142076, 0.0136367, 0.0126589, 0.012195, 0.0112466, 0.0102815, 0.00954721, 0.00920796, 0.00848258, 0.00804014, 0.00762756, 0.0071216, 0.00727543, 0.00694033, 0.00698688, 0.00702847, 0.0068209, 0.00682605, 0.00684432, 0.00685966, 0.00686981, 0.00704506, 0.00730533])

xs_sys_all = array('d', [2.0615, 1.94582, 2.06492, 1.62597, 1.05113, 0.619752, 0.409027, 0.278679, 0.223427, 0.231844, 0.213956, 0.17448, 0.182167, 0.187449, 0.181862, 0.17469, 0.163485, 0.151662, 0.139721, 0.129114, 0.120302, 0.11357, 0.100043, 0.0874448, 0.081939, 0.0751703, 0.0658594, 0.0515029, 0.0351077, 0.0223095, 0.0158919, 0.0128366, 0.0113868, 0.0106171, 0.0108656, 0.0121578, 0.0134528, 0.0142962, 0.0143024, 0.0139918, 0.0132415, 0.0124992, 0.0118243, 0.0110701, 0.0101539, 0.00942484, 0.00875206, 0.00818491, 0.00757354, 0.00701268, 0.0064779, 0.00638241, 0.00611084, 0.00592149, 0.0057839, 0.00583632, 0.00589959, 0.00615538])
xs_sys_all_exp = array('d', [2.35524, 1.718395, 1.428665, 1.09346, 0.8535145, 0.6807505, 0.5645625, 0.4462245, 0.3892175, 0.338744, 0.296201, 0.261048, 0.227171, 0.203175, 0.18259, 0.1602235, 0.14302, 0.125324, 0.1061265, 0.09317755, 0.0843383, 0.07174945, 0.064988, 0.05691025, 0.0509401, 0.0453084, 0.04059635, 0.0373627, 0.0317369, 0.0300091, 0.0270188, 0.02477945, 0.0223704, 0.0199402, 0.0187306, 0.01701985, 0.0151752, 0.01351145, 0.01255925, 0.01157715, 0.0106024, 0.009756425, 0.00934829, 0.00875005, 0.00851081, 0.008080615, 0.00765728, 0.007629135, 0.007262815, 0.007117835, 0.007197575, 0.00695401, 0.006931495, 0.00696994, 0.00698784, 0.007031325, 0.007249375, 0.00743711])


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
g_all_exp.GetXaxis().SetTitle("gg resonance mass [GeV]")
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

c.SaveAs('xs_limit_ratio_DijetLimitCode_gg_Run2_13TeV_DATA_2445_invpb.eps')
