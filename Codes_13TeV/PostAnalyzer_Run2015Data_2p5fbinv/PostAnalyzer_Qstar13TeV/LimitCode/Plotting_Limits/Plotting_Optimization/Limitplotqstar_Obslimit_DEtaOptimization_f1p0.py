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

xs_obs_limits_DEta_1p0 = array('d', [0.0342358, 0.0356627, 0.0211046, 0.0100935, 0.0210988, 0.0158, 0.00689509, 0.00411782, 0.00260126, 0.00280136, 0.00260232, 0.00189327, 0.00155613, 0.0014512, 0.00136332, 0.0013196, 0.00131451, 0.00127745, 0.00124765, 0.00125209, 0.00126532, 0.00125505, 0.00125493])

xs_obs_limits_DEta_1p2 = array('d', [0.0373834, 0.0376027, 0.0224994, 0.0113554, 0.020893, 0.0186862, 0.00785202, 0.00439992, 0.00275128, 0.00275036, 0.00253032, 0.00185464, 0.00151145, 0.00140891, 0.00135684, 0.0013376, 0.00130222, 0.00128581, 0.00126573, 0.00124648, 0.00125829, 0.00126057, 0.00123793])

xs_obs_limits_DEta_1p5 = array('d', [0.0429857, 0.0394882, 0.0243286, 0.016551, 0.0223863, 0.0205461, 0.00895466, 0.00561192, 0.00325719, 0.00292358, 0.00271489, 0.00228476, 0.00180051, 0.00153467, 0.00141721, 0.00136463, 0.00133804, 0.0013118, 0.0012878, 0.00126209, 0.0012948, 0.00125725, 0.00126219])

#these numbers for DEta=2p0 from varun (Read README)
xs_obs_limits_DEta_2p0 = array('d', [0.0400255, 0.0335698, 0.0366197, 0.0309246, 0.0323951, 0.0275996, 0.0134593, 0.00852663, 0.00386713, 0.00356322, 0.00333334, 0.002743, 0.00223023, 0.00214148, 0.00196558, 0.00171071, 0.00153965, 0.00146057, 0.00140387, 0.00137167, 0.00137552, 0.00136663, 0.0014239])

xs_obs_limits_DEta_2p2 = array('d', [0.0710324, 0.0303067, 0.0410112, 0.0255984, 0.0248121, 0.0247373, 0.0115842, 0.00758625, 0.00367181, 0.00327551, 0.00313894, 0.00267846, 0.00219208, 0.00209202, 0.00196737, 0.0017157, 0.00153134, 0.00144477, 0.00139459, 0.00135946, 0.00134918, 0.00134156, 0.00132523])


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


eff_qstar_DEta_1p0     = array('d', [ 0.31241,   0.31752,   0.32325,   0.32910,   0.33498,   0.34173,   0.34756,   0.34314,   0.35919,   0.36474,  
                                      0.36959,   0.37056,   0.37195,   0.37323,   0.37498,   0.37745,   0.37905,   0.38036,   0.38212,   0.38345,
                                      0.38423,   0.38445,   0.38464,   0.38487,   0.38502,   0.38534,   0.38557,   0.38579,   0.38599,   0.38601,
                                      0.38628,   0.38640,   0.38654,   0.38662,   0.38674,   0.38687,   0.38698,   0.38711,   0.38724,   0.38741,
                                      0.38752])

eff_qstar_DEta_1p2     = array('d', [0.36104,   0.36774,    0.37432,   0.37998,   0.38687,   0.39123,   0.39886,   0.40512,   0.41145,  0.41867,
                                     0.42690,   0.42923,    0.43129,   0.43345,   0.43534,   0.43723,   0.43934,   0.44112,   0.44324,  0.44487,
                                     0.44548,   0.44556,    0.44564,   0.44572,   0.44580,   0.44587,   0.44597,   0.44605,   0.44613,  0.44619,
                                     0.44627,   0.44633,    0.44639,   0.44645,   0.44651,   0.44657,   0.44663,   0.44670,   0.44676,  0.44682,
                                     0.44688])

eff_qstar_DEta_1p5     = array('d', [0.42083,   0.42867,   0.43678,   0.44456,    0.45234,   0.46057,   0.46879,   0.47689,   0.48502,  0.49334,
                                     0.50099,   0.50312,   0.50498,   0.50704,    0.50923,   0.51145,   0.51367,   0.51545,   0.51789,  0.51894, 
                                     0.52068,   0.52102,   0.52142,   0.52174,    0.52209,   0.52251,   0.52299,   0.52335,   0.52398,  0.52448,
                                     0.52491,   0.52504,   0.52516,   0.52532,    0.52545,   0.52559,   0.52571,   0.52586,   0.52599,  0.52613,
                                     0.52628])
      

#These eff numbers for DEta=2p0 from varun (Read README)
eff_qstar_DEta_2p0     = array('d',  [ 0.488302,  0.49915,   0.51,      0.52,     0.53,      0.54,      0.55,      0.56,      0.57,      0.5793, 
                                       0.58861,   0.59115,   0.5937,    0.59634,   0.599,     0.60195,   0.6049,    0.60745,   0.61,      0.5993,
                                       0.615333,  0.61636,   0.6174,    0.61845,   0.6195,    0.62065,   0.6218,    0.62295,   0.6241,    0.62511,
                                       0.626118,  0.62584,   0.625565,  0.625537,  0.62551,   0.62538,   0.62525,   0.625065,  0.62488,   0.62474,
                                       0.624607])

eff_qstar_DEta_2p2     = array('d', [0.50496,   0.51503,   0.52592,   0.53654,     0.54739,   0.55816,   0.56904,   0.57929,  0.59145,   0.60234,
                                     0.61446,   0.61752,   0.62035,   0.62319,     0.62588,   0.62823,   0.63027,   0.63303,  0.63578,   0.63845,   
                                     0.64136,   0.64256,   0.64381,   0.64502,     0.64623,   0.64753,   0.64882,   0.65001,  0.65140,   0.65266,    
                                     0.65389,   0.65385,   0.65381,   0.65377,     0.65372,   0.65367,   0.65363,   0.65359,  0.65355,   0.65351, 
                                     0.65347])


xstimesEff_DEta_1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DEta_1p0)])
xstimesEff_DEta_1p2 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DEta_1p2)])
xstimesEff_DEta_1p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DEta_1p5)])
xstimesEff_DEta_2p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DEta_2p0)])
xstimesEff_DEta_2p2 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_DEta_2p2)])

#print xstimesEff 

graph_obs_DEta_1p0 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_DEta_1p0)
graph_obs_DEta_1p0.GetXaxis().SetTitle("q* Mass [TeV]")
graph_obs_DEta_1p0.GetYaxis().SetTitle("#sigma #times BR #times A #times #epsilon [pb]")
#graph_obs_DEta_1p0.GetYaxis().SetRangeUser(1e-4,1)
graph_obs_DEta_1p0.GetYaxis().SetRangeUser(5e-4,5e-3)
graph_obs_DEta_1p0.GetXaxis().SetNdivisions(510)
#graph_obs_DEta_1p0.GetXaxis().SetLimits(0.6,6.0) 
graph_obs_DEta_1p0.GetXaxis().SetLimits(3.8,4.9) 

graph_obs_DEta_1p0.GetYaxis().CenterTitle()
graph_obs_DEta_1p0.GetYaxis().SetLabelSize(0.05)
graph_obs_DEta_1p0.GetYaxis().SetTitleOffset(1.1)
graph_obs_DEta_1p0.GetXaxis().CenterTitle()
graph_obs_DEta_1p0.GetXaxis().SetLabelSize(0.05)
graph_obs_DEta_1p0.GetXaxis().SetTitleOffset(1.1)

graph_obs_DEta_1p0.SetMarkerStyle(24)
graph_obs_DEta_1p0.SetMarkerColor(12)
graph_obs_DEta_1p0.SetMarkerSize(0.7)
graph_obs_DEta_1p0.SetLineWidth(2)
graph_obs_DEta_1p0.SetLineColor(12)

graph_obs_DEta_1p2 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_DEta_1p2)
graph_obs_DEta_1p2.SetMarkerStyle(25)
graph_obs_DEta_1p2.SetMarkerColor(3)
graph_obs_DEta_1p2.SetMarkerSize(0.7)
graph_obs_DEta_1p2.SetLineWidth(2)
graph_obs_DEta_1p2.SetLineColor(3)

graph_obs_DEta_1p5 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_DEta_1p5)
graph_obs_DEta_1p5.SetMarkerStyle(26)
graph_obs_DEta_1p5.SetMarkerColor(6)
graph_obs_DEta_1p5.SetMarkerSize(0.7)
graph_obs_DEta_1p5.SetLineWidth(2)
graph_obs_DEta_1p5.SetLineColor(6)

graph_obs_DEta_2p0 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_DEta_2p0)
graph_obs_DEta_2p0.SetMarkerStyle(27)
graph_obs_DEta_2p0.SetMarkerColor(2)
graph_obs_DEta_2p0.SetMarkerSize(0.7)
graph_obs_DEta_2p0.SetLineWidth(2)
graph_obs_DEta_2p0.SetLineColor(2)

graph_obs_DEta_2p2 = TGraph(len(masses_tev), masses_tev, xs_obs_limits_DEta_2p2)
graph_obs_DEta_2p2.SetMarkerStyle(28)
graph_obs_DEta_2p2.SetMarkerColor(4)
graph_obs_DEta_2p2.SetMarkerSize(0.7)
graph_obs_DEta_2p2.SetLineWidth(2)
graph_obs_DEta_2p2.SetLineColor(4)


graph_qstar_DEta_1p0 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DEta_1p0)
graph_qstar_DEta_1p0.SetLineWidth(2)
graph_qstar_DEta_1p0.SetLineColor(12)
graph_qstar_DEta_1p0.SetLineStyle(2)

graph_qstar_DEta_1p2 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DEta_1p2)
graph_qstar_DEta_1p2.SetLineWidth(2)
graph_qstar_DEta_1p2.SetLineColor(3)
graph_qstar_DEta_1p2.SetLineStyle(2)

graph_qstar_DEta_1p5 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DEta_1p5)
graph_qstar_DEta_1p5.SetLineWidth(2)
graph_qstar_DEta_1p5.SetLineColor(6)
graph_qstar_DEta_1p5.SetLineStyle(2)

graph_qstar_DEta_2p0 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DEta_2p0)
graph_qstar_DEta_2p0.SetLineWidth(2)
graph_qstar_DEta_2p0.SetLineColor(2)
graph_qstar_DEta_2p0.SetLineStyle(2)

graph_qstar_DEta_2p2 = TGraph(len(masses_qstar), masses_qstar, xstimesEff_DEta_2p2)
graph_qstar_DEta_2p2.SetLineWidth(2)
graph_qstar_DEta_2p2.SetLineColor(4)
graph_qstar_DEta_2p2.SetLineStyle(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_obs_DEta_1p0.Draw("ALP")
graph_obs_DEta_1p2.Draw("LP")
graph_obs_DEta_1p5.Draw("LP")
graph_obs_DEta_2p0.Draw("LP")
graph_obs_DEta_2p2.Draw("LP")

graph_qstar_DEta_1p0.Draw("L")
graph_qstar_DEta_1p2.Draw("L")
graph_qstar_DEta_1p5.Draw("L")
graph_qstar_DEta_2p0.Draw("L")
graph_qstar_DEta_2p2.Draw("L")


legend = TLegend(.55,.69,.85,.92)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.025)
legend.SetHeader('95% CL Upper Expected Limits')
legend.AddEntry(graph_obs_DEta_1p0, "Observed (#Delta #eta < 1.0)","lp" )
legend.AddEntry(graph_obs_DEta_1p2, "Observed (#Delta #eta < 1.2)","lp")
legend.AddEntry(graph_obs_DEta_1p5, "Observed (#Delta #eta < 1.5)","lp")
legend.AddEntry(graph_obs_DEta_2p0, "Observed (#Delta #eta < 2.0)","lp")
legend.AddEntry(graph_obs_DEta_2p2, "Observed (#Delta #eta < 2.2)","lp")

legend.AddEntry(graph_qstar_DEta_1p0, "q* (#Delta #eta < 1.0)", "l")
legend.AddEntry(graph_qstar_DEta_1p2, "q* (#Delta #eta < 1.2)", "l")
legend.AddEntry(graph_qstar_DEta_1p5, "q* (#Delta #eta < 1.5)", "l")
legend.AddEntry(graph_qstar_DEta_2p0, "q* (#Delta #eta < 2.0)", "l")
legend.AddEntry(graph_qstar_DEta_2p2, "q* (#Delta #eta < 2.2)", "l")
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

c.SetLogy()
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ObsVsQstar.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ObsVsQstar.png')

c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ObsVsQstar_zoomed.pdf')
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ObsVsQstar_zoomed.png')



