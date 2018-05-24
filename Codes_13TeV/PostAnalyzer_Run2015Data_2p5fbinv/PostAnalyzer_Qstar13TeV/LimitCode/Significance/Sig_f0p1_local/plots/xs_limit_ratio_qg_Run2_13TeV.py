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

xs_stat = array('d', [0.602272, 0.64491, 0.865673, 0.616196, 0.343269, 0.212624, 0.140649, 0.103605, 0.0918694, 0.101727, 0.0892428, 0.0678928, 0.0825205, 0.0915746, 0.0943468, 0.0911699, 0.0866712, 0.0832746, 0.0744257, 0.0706772, 0.0713686, 0.065973, 0.0549593, 0.0506374, 0.0501734, 0.0443592, 0.0335244, 0.0213411, 0.0135937, 0.0100408, 0.00828029, 0.00734691, 0.00693662, 0.00679592, 0.00720566, 0.00848763, 0.00966361, 0.0100682, 0.0101104, 0.00961391, 0.00901056, 0.00843457, 0.00787787, 0.007268, 0.0068383, 0.0062159, 0.00564396, 0.00512824, 0.00467353, 0.00421974, 0.00397535, 0.00369705, 0.00362929, 0.00343399, 0.00331401, 0.00323391, 0.00324306, 0.00327224])
xs_stat_exp = array('d', [0.604644, 0.474845, 0.431861, 0.371758, 0.309556, 0.251925, 0.210053, 0.196129, 0.167196, 0.158303, 0.131379, 0.108406, 0.100809, 0.0946049, 0.0795457, 0.0718638, 0.0629254, 0.0589563, 0.0533615, 0.047158, 0.0424394, 0.0381913, 0.0333004, 0.0309208, 0.0274924, 0.024903, 0.0237022, 0.0209057, 0.0194533, 0.0179081, 0.0165851, 0.0141789, 0.0139069, 0.0127362, 0.0119621, 0.0104715, 0.00996637, 0.00912483, 0.00839968, 0.0077472, 0.00685672, 0.00680356, 0.00636493, 0.0060644, 0.00583616, 0.00540399, 0.00510707, 0.00492442, 0.0046922, 0.00473624, 0.00445597, 0.00434573, 0.00420864, 0.00409594, 0.00399577, 0.00389151, 0.00387871, 0.00392686])

xs_sys_all = array('d', [1.3959, 1.44829, 1.31726, 1.00642, 0.637327, 0.378596, 0.248195, 0.179449, 0.156801, 0.157078, 0.140802, 0.123014, 0.126673, 0.129163, 0.12916, 0.123556, 0.115661, 0.10677, 0.097489, 0.090278, 0.0860718, 0.0789864, 0.0697045, 0.0619197, 0.0574263, 0.0508895, 0.0431843, 0.0315208, 0.0207677, 0.0137897, 0.0105188, 0.00888977, 0.00784571, 0.00769529, 0.00820165, 0.00918646, 0.0099582, 0.0102429, 0.0102928, 0.00988389, 0.00943999, 0.00886349, 0.0083116, 0.00767947, 0.00708125, 0.00644123, 0.00599399, 0.00549601, 0.0049915, 0.00448342, 0.00424979, 0.00392172, 0.00367954, 0.00349917, 0.00334032, 0.00325643, 0.00322229, 0.00324652])
xs_sys_all_exp = array('d', [1.5228, 1.141925, 0.9300815, 0.689262, 0.546956, 0.444851, 0.362029, 0.3022745, 0.2590265, 0.2264055, 0.202871, 0.173006, 0.1562905, 0.140351, 0.126195, 0.1108165, 0.09769975, 0.0839188, 0.07267645, 0.06334305, 0.0558337, 0.05026125, 0.04454515, 0.0396875, 0.03530295, 0.0312174, 0.0287047, 0.02557055, 0.0228086, 0.020473, 0.01940965, 0.01746235, 0.01593015, 0.01452195, 0.01314315, 0.0118794, 0.0107423, 0.009848545, 0.00888611, 0.00831464, 0.00746645, 0.006944, 0.00649144, 0.00627163, 0.006093345, 0.005412895, 0.005109905, 0.00507699, 0.00489189, 0.0046593, 0.00453367, 0.004363205, 0.0042402, 0.004188755, 0.004090695, 0.00398848, 0.00398725, 0.00402766])


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
g_all_exp.GetXaxis().SetTitle("qg resonance mass [GeV]")
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

c.SaveAs('xs_limit_ratio_DijetLimitCode_qg_Run2_13TeV_DATA_2445_invpb.eps')
