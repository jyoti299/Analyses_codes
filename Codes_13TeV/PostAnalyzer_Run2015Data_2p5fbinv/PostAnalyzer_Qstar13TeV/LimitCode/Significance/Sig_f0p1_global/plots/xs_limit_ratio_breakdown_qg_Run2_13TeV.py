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

xs_sys_all = array('d', [1.3959, 1.44829, 1.31726, 1.00642, 0.637327, 0.378596, 0.248195, 0.179449, 0.156801, 0.157078, 0.140802, 0.123014, 0.126673, 0.129163, 0.12916, 0.123556, 0.115661, 0.10677, 0.097489, 0.090278, 0.0860718, 0.0789864, 0.0697045, 0.0619197, 0.0574263, 0.0508895, 0.0431843, 0.0315208, 0.0207677, 0.0137897, 0.0105188, 0.00888977, 0.00784571, 0.00769529, 0.00820165, 0.00918646, 0.0099582, 0.0102429, 0.0102928, 0.00988389, 0.00943999, 0.00886349, 0.0083116, 0.00767947, 0.00708125, 0.00644123, 0.00599399, 0.00549601, 0.0049915, 0.00448342, 0.00424979, 0.00392172, 0.00367954, 0.00349917, 0.00334032, 0.00325643, 0.00322229, 0.00324652])
xs_sys_lumi = array('d', [0.594873, 0.648152, 0.872382, 0.615117, 0.342542, 0.212035, 0.138339, 0.0995874, 0.0883757, 0.096518, 0.0864851, 0.0670714, 0.0808903, 0.0887106, 0.0930663, 0.0906983, 0.08704, 0.0834016, 0.0746245, 0.0697667, 0.071423, 0.066049, 0.0551166, 0.0507978, 0.0504566, 0.0446149, 0.0333252, 0.0209162, 0.0134808, 0.0097246, 0.00804557, 0.00717518, 0.00665063, 0.00650089, 0.00713424, 0.00845211, 0.00953664, 0.00999928, 0.0100676, 0.00957548, 0.009004, 0.00844683, 0.0079343, 0.00722663, 0.00663336, 0.00605794, 0.00551934, 0.00501037, 0.00456732, 0.00415564, 0.00391043, 0.00364222, 0.0034475, 0.00328219, 0.0031706, 0.00311496, 0.00313076, 0.00314383])
xs_sys_jes = array('d', [0.604779, 0.660484, 0.863743, 0.6346, 0.356208, 0.218863, 0.140999, 0.100439, 0.0899677, 0.0965544, 0.08684, 0.0696732, 0.0816664, 0.0898407, 0.0934652, 0.0917876, 0.0875332, 0.0841347, 0.0767135, 0.071816, 0.0716696, 0.0668486, 0.0575484, 0.0527056, 0.0509839, 0.0458684, 0.0371859, 0.0253632, 0.0159672, 0.0109505, 0.00861923, 0.00753376, 0.00691684, 0.0067643, 0.007351, 0.00850387, 0.00945971, 0.00982862, 0.00989599, 0.00964673, 0.00914131, 0.00856937, 0.00799325, 0.00733573, 0.00689399, 0.00626784, 0.00579177, 0.00525806, 0.00488021, 0.00439661, 0.00411143, 0.00382207, 0.00362227, 0.00341613, 0.00328506, 0.00319744, 0.00315642, 0.00318087])
xs_sys_jer = array('d', [0.593902, 0.658109, 0.869302, 0.611817, 0.340045, 0.213608, 0.138463, 0.0992404, 0.0881938, 0.0976631, 0.0858611, 0.0693968, 0.0810925, 0.0892428, 0.0942916, 0.0911835, 0.0873916, 0.0834038, 0.0758937, 0.0710476, 0.0719659, 0.0670417, 0.0563951, 0.0509355, 0.0502198, 0.0446204, 0.0327828, 0.0208577, 0.0134656, 0.00990647, 0.00807757, 0.00730242, 0.00671285, 0.00655493, 0.00724514, 0.0084301, 0.00961693, 0.0100016, 0.0100887, 0.00957921, 0.00901861, 0.00837822, 0.00790313, 0.00728371, 0.0066288, 0.00599279, 0.0054422, 0.00501238, 0.00459157, 0.00412951, 0.00389859, 0.00367268, 0.00351432, 0.00329745, 0.00323237, 0.00313178, 0.00314195, 0.00317233])
xs_sys_allexceptbkg = array('d', [0.606518, 0.67277, 0.881931, 0.636914, 0.360274, 0.222554, 0.142504, 0.101661, 0.0910753, 0.0980752, 0.0872266, 0.0716921, 0.0827352, 0.090132, 0.0944572, 0.092679, 0.0890855, 0.084805, 0.0780477, 0.0732637, 0.0725176, 0.0680824, 0.0589475, 0.0537052, 0.0517402, 0.0462893, 0.0372981, 0.0255326, 0.016154, 0.0109968, 0.00878908, 0.00760465, 0.00694678, 0.00684742, 0.0074726, 0.00859647, 0.00944546, 0.0097633, 0.00995764, 0.0096666, 0.00919886, 0.00868526, 0.0080872, 0.00747738, 0.00684485, 0.00635919, 0.00575654, 0.00530048, 0.0049206, 0.00445544, 0.00420546, 0.00390629, 0.00367791, 0.00344523, 0.00333409, 0.00320046, 0.00323719, 0.00324114])
xs_sys_bkg = array('d', [1.32749, 1.26908, 1.26506, 0.887255, 0.532602, 0.342323, 0.234023, 0.170733, 0.148354, 0.154467, 0.137864, 0.110305, 0.121266, 0.123918, 0.12591, 0.118193, 0.108704, 0.102976, 0.0905079, 0.0835408, 0.0815178, 0.0743441, 0.0625858, 0.0568327, 0.054615, 0.047838, 0.0362492, 0.0240357, 0.0157822, 0.0115733, 0.00942389, 0.00834118, 0.00743891, 0.0073188, 0.00775025, 0.0089576, 0.0100234, 0.0104213, 0.0103883, 0.00980813, 0.00920145, 0.00859254, 0.00801442, 0.00731748, 0.00672396, 0.00610263, 0.00556593, 0.00503144, 0.00461667, 0.00418602, 0.00393147, 0.00369379, 0.00348629, 0.00334396, 0.00321781, 0.00312172, 0.00311622, 0.00314165])


r_all = array('d')
r_lumi = array('d')
r_jes = array('d')
r_jer = array('d')
r_allexceptbkg = array('d')
r_bkg = array('d')


for i in range(0,len(xs_stat)):
  r_all.append( xs_sys_all[i] / xs_stat[i] )
  r_lumi.append( xs_sys_lumi[i] / xs_stat[i] )
  r_jes.append( xs_sys_jes[i] / xs_stat[i] )
  r_jer.append( xs_sys_jer[i] / xs_stat[i] )
  r_allexceptbkg.append( xs_sys_allexceptbkg[i] / xs_stat[i] )
  r_bkg.append( xs_sys_bkg[i] / xs_stat[i] )

g_all = TGraph(len(masses),masses,r_all)
g_all.SetMarkerStyle(20)
g_all.SetMarkerColor(kBlack)
g_all.SetLineWidth(2)
g_all.SetLineStyle(1)
g_all.SetLineColor(kBlack)
g_all.GetXaxis().SetTitle("qg resonance mass [GeV]")
g_all.GetYaxis().SetTitle("Limit ratio")
g_all.GetYaxis().SetTitleOffset(1.1)
g_all.GetYaxis().SetRangeUser(0.5,3.)
#g_all.GetXaxis().SetNdivisions(1005)

g_lumi = TGraph(len(masses),masses,r_lumi)
g_lumi.SetMarkerStyle(24)
g_lumi.SetMarkerColor(kGreen+2)
g_lumi.SetLineWidth(2)
g_lumi.SetLineStyle(2)
g_lumi.SetLineColor(kGreen+2)

g_jes = TGraph(len(masses),masses,r_jes)
g_jes.SetMarkerStyle(25)
g_jes.SetMarkerColor(kBlue)
g_jes.SetLineWidth(2)
g_jes.SetLineStyle(3)
g_jes.SetLineColor(kBlue)

g_jer = TGraph(len(masses),masses,r_jer)
g_jer.SetMarkerStyle(26)
g_jer.SetMarkerColor(45)
g_jer.SetLineWidth(2)
g_jer.SetLineStyle(4)
g_jer.SetLineColor(45)

g_allexceptbkg = TGraph(len(masses),masses,r_allexceptbkg)
g_allexceptbkg.SetMarkerStyle(27)
g_allexceptbkg.SetMarkerColor(kMagenta)
g_allexceptbkg.SetLineWidth(2)
g_allexceptbkg.SetLineStyle(5)
g_allexceptbkg.SetLineColor(kMagenta)

g_bkg = TGraph(len(masses),masses,r_bkg)
g_bkg.SetMarkerStyle(30)
g_bkg.SetMarkerColor(kRed)
g_bkg.SetLineWidth(2)
g_bkg.SetLineStyle(6)
g_bkg.SetLineColor(kRed)

c = TCanvas("c", "",800,800)
c.cd()

g_all.Draw("ALP")
g_lumi.Draw("LP")
g_jes.Draw("LP")
g_jer.Draw("LP")
g_allexceptbkg.Draw("LP")
g_bkg.Draw("LP")

legend = TLegend(.40,.56,.60,.80)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetMargin(0.20)
legend.AddEntry(g_all, "All / Stat. only","lp")
legend.AddEntry(g_lumi, "Lumi / Stat. only","lp")
legend.AddEntry(g_jes, "JES / Stat. only","lp")
legend.AddEntry(g_jer, "JER / Stat. only","lp")
legend.AddEntry(g_allexceptbkg, "All except bkg / Stat. only","lp")
legend.AddEntry(g_bkg, "Bkg / Stat. only","lp")

legend.Draw()

#draw the lumi text on the canvas
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "2.4 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

CMS_lumi.CMS_lumi(c, iPeriod, iPos)

c.SetGridx()
c.SetGridy()

c.SaveAs('xs_limit_ratio_breakdown_DijetLimitCode_qg_Run2_13TeV_DATA_2445_invpb.eps')
