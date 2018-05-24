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


masses_gev = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0])
masses_tev = array('d', [ 0.001 * float(b) for b in masses_gev ])

xs_obs_limits  = array('d', [0.0400255, 0.0335698, 0.0366197, 0.0309246, 0.0323951, 0.0275996, 0.0134593, 0.00852663, 0.00386713, 0.00356322, 0.00333334, 0.002743, 0.00223023, 0.00214148, 0.00196558, 0.00171071, 0.00153965, 0.00146057, 0.00140387, 0.00137167, 0.00137552, 0.00136663, 0.0014239])
xs_obs_limits_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits ])

xs_exp_limits = array('d', [0.0753898, 0.0481741, 0.0331191, 0.0228193, 0.0168115, 0.0110654, 0.00836895, 0.00696791, 0.00545557, 0.00469033, 0.00389036, 0.00333766, 0.00296295, 0.00251753, 0.00229563, 0.00213392, 0.00202398, 0.00186637, 0.00176224, 0.00162506, 0.00160681, 0.00155983, 0.00151451])
xs_exp_limits_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits ])

masses_exp_gev = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0, 5400.0, 5200.0, 5000.0, 4800.0, 4600.0, 4400.0, 4200.0, 4000.0, 3800.0, 3600.0, 3400.0, 3200.0, 3000.0, 2800.0, 2600.0, 2400.0, 2200.0, 2000.0, 1800.0, 1600.0, 1400.0, 1200.0, 1000.0])
masses_exp_tev = array('d', [ 0.001 * float(b) for b in masses_exp_gev ])

xs_exp_limits_1sigma  = array('d', [0.0556791, 0.034296, 0.0226546, 0.0168, 0.0122775, 0.00834437, 0.00585152, 0.00488808, 0.00400916, 0.00347334, 0.00292542, 0.00246554, 0.00215404, 0.00190562, 0.00173677, 0.00157048, 0.00148674, 0.00145768, 0.00140231, 0.0013732, 0.00136946, 0.0013614, 0.00135007, 0.00200538, 0.00207958, 0.00212182, 0.00209974, 0.00228635, 0.00247011, 0.0026456, 0.00295388, 0.00303269, 0.00341805, 0.00413127, 0.00470363, 0.00560917, 0.00713875, 0.00843934, 0.00947565, 0.0124935, 0.0154206, 0.0238926, 0.033745, 0.0456901, 0.0637495, 0.104935])
xs_exp_limits_1sigma_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits_1sigma ])

xs_exp_limits_2sigma = array('d', [0.0402204, 0.0229533, 0.0190181, 0.0126229, 0.00898654, 0.0055623, 0.0045272, 0.0037742, 0.00333205, 0.00281111, 0.00227257, 0.00183896, 0.0017639, 0.00153834, 0.00143864, 0.00139636, 0.00133767, 0.00133109, 0.0012932, 0.0012943, 0.00129601, 0.00128195, 0.00129225, 0.00248205, 0.00247471, 0.00268562, 0.00269724, 0.00302142, 0.00330211, 0.00347472, 0.0038537, 0.00392481, 0.00454296, 0.00533445, 0.00613284, 0.00811228, 0.00881012, 0.0104857, 0.013546, 0.0160148, 0.0201772, 0.032697, 0.0421027, 0.0670475, 0.0852177, 0.13611])
xs_exp_limits_2sigma_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits_2sigma ])


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

xs_qstar_f0p5     = array('d', [ 4.137e+00, 2.642e+00, 1.768e+00, 1.217e+00, 8.445e-01, 6.012e-01, 4.345e-01, 3.179e-01, 2.342e-01, 1.765e-01,
                                 1.328e-01, 1.005e-01, 7.712e-02, 5.922e-02, 4.583e-02, 3.601e-02, 2.799e-02, 2.206e-02, 1.746e-02, 1.378e-02,
				 1.096e-02, 8.642e-03, 7.002e-03, 5.531e-03, 4.407e-03, 3.554e-03, 2.860e-03, 2.302e-03, 1.851e-03, 1.488e-03,
				 1.211e-03, 9.753e-04, 7.847e-04, 6.374e-04, 5.156e-04, 4.187e-04, 3.360e-04, 2.728e-04, 2.189e-04, 1.770e-04,
				 1.437e-04])

xs_qstar_f0p1     = array('d', [ 1.655e-01, 1.057e-01, 7.134e-02, 4.932e-02, 3.421e-02, 2.440e-02, 1.750e-02, 1.284e-02, 9.433e-03, 7.075e-03,
                                 5.298e-03, 4.025e-03, 3.107e-03, 2.374e-03, 1.861e-03, 1.431e-03, 1.130e-03, 8.902e-04, 7.051e-04, 5.527e-04,
				 4.363e-04, 3.511e-04, 2.784e-04, 2.232e-04, 1.791e-04, 1.435e-04, 1.148e-04, 9.267e-05, 7.459e-05, 6.014e-05,
				 4.852e-05, 3.902e-05, 3.157e-05, 2.536e-05, 2.058e-05, 1.677e-05, 1.344e-05, 1.087e-05, 8.690e-06, 7.102e-06, 
				 5.739e-06])

###------------------------------------------------------------------------------------------------
xs_qstar_f0p2     = array('d', [ 6.639e-01, 4.269e-01, 2.888e-01, 1.947e-01, 1.356e-01, 9.700e-02, 6.970e-02, 5.095e-02, 3.808e-02, 2.833e-02,
                                 2.153e-02, 1.628e-02, 1.234e-02, 9.542e-03, 7.398e-03, 5.779e-03, 4.491e-03, 3.534e-03, 2.783e-03, 2.214e-03,
				 1.764e-03, 1.388e-03, 1.120e-03, 8.970e-04, 7.071e-04, 5.739e-04, 4.573e-04, 3.708e-04, 2.976e-04, 2.406e-04,
				 1.945e-04, 1.563e-04, 1.268e-04, 1.021e-04, 8.215e-05, 6.647e-05, 5.349e-05, 4.343e-05, 3.515e-05, 2.863e-05,
				 2.259e-05])

xs_qstar_f0p3     = array('d', [ 1.485e+00, 9.758e-01, 6.363e-01, 4.373e-01, 3.095e-01, 2.172e-01, 1.565e-01, 1.142e-01, 8.533e-02, 6.345e-02,
                                 4.800e-02, 3.649e-02, 2.772e-02, 2.126e-02, 1.665e-02, 1.292e-02, 1.010e-02, 7.987e-03, 6.289e-03, 4.987e-03,
				 3.964e-03, 3.147e-03, 2.495e-03, 2.002e-03, 1.599e-03, 1.292e-03, 1.038e-03, 8.221e-04, 6.671e-04, 5.420e-04,
				 4.346e-04, 3.484e-04, 2.854e-04, 2.303e-04, 1.861e-04, 1.481e-04, 1.199e-04, 9.773e-05, 7.918e-05, 6.399e-05,
				 5.114e-05])

xs_qstar_f0p4     = array('d', [ 2.620e+00, 1.693e+00, 1.145e+00, 7.827e-01, 5.403e-01, 3.810e-01, 2.779e-01, 2.037e-01, 1.507e-01, 1.122e-01,
                                 8.484e-02, 6.493e-02, 4.930e-02, 3.778e-02, 2.942e-02, 2.286e-02, 1.791e-02, 1.425e-02, 1.116e-02, 8.801e-03,
				 6.993e-03, 5.556e-03, 4.401e-03, 3.517e-03, 2.837e-03, 2.296e-03, 1.841e-03, 1.475e-03, 1.196e-03, 9.596e-04,
				 7.733e-04, 6.208e-04, 5.032e-04, 4.069e-04, 3.287e-04, 2.677e-04, 2.129e-04, 1.726e-04, 1.397e-04, 1.139e-04,
				 9.133e-05])

xs_qstar_f0p6     = array('d', [ 5.885e+00, 3.825e+00, 2.559e+00, 1.738e+00, 1.225e+00, 8.697e-01, 6.238e-01, 4.568e-01, 3.357e-01, 2.521e-01,
                                 1.898e-01, 1.437e-01, 1.095e-01, 8.543e-02, 6.586e-02, 5.122e-02, 3.982e-02, 3.130e-02, 2.499e-02, 1.956e-02,
				 1.550e-02, 1.261e-02, 9.986e-03, 7.940e-03, 6.315e-03, 5.108e-03, 4.122e-03, 3.326e-03, 2.677e-03, 2.158e-03,
				 1.730e-03, 1.409e-03, 1.120e-03, 9.219e-04, 7.393e-04, 5.939e-04, 4.837e-04, 3.926e-04, 3.166e-04, 2.585e-04,
				 2.086e-04])

xs_qstar_f0p7     = array('d', [ 7.946e+00, 5.163e+00, 3.483e+00, 2.394e+00, 1.656e+00, 1.174e+00, 8.550e-01, 6.142e-01, 4.576e-01, 3.419e-01,
                                 2.583e-01, 1.961e-01, 1.482e-01, 1.167e-01, 8.892e-02, 6.964e-02, 5.433e-02, 4.290e-02, 3.364e-02, 2.680e-02,
				 2.125e-02, 1.693e-02, 1.358e-02, 1.071e-02, 8.698e-03, 6.930e-03, 5.574e-03, 4.434e-03, 3.635e-03, 2.918e-03,
				 2.343e-03, 1.893e-03, 1.542e-03, 1.249e-03, 1.009e-03, 8.111e-04, 6.582e-04, 5.356e-04, 4.311e-04, 3.499e-04,
				 2.817e-04])

xs_qstar_f0p8     = array('d', [ 1.051e+01, 6.789e+00, 4.466e+00, 3.078e+00, 2.154e+00, 1.538e+00, 1.106e+00, 8.069e-01, 5.905e-01, 4.475e-01,
                                 3.363e-01, 2.528e-01, 1.954e-01, 1.485e-01, 1.154e-01, 9.041e-02, 7.069e-02, 5.579e-02, 4.415e-02, 3.467e-02,
				 2.756e-02, 2.197e-02, 1.735e-02, 1.406e-02, 1.120e-02, 8.976e-03, 7.284e-03, 5.850e-03, 4.731e-03, 3.816e-03,
				 3.083e-03, 2.478e-03, 2.018e-03, 1.634e-03, 1.330e-03, 1.068e-03, 8.581e-04, 6.984e-04, 5.651e-04, 4.593e-04,
				 3.693e-04])

xs_qstar_f0p9     = array('d', [ 1.324e+01, 8.468e+00, 5.714e+00, 3.892e+00, 2.730e+00, 1.929e+00, 1.397e+00, 1.018e+00, 7.498e-01, 5.616e-01,
                                 4.202e-01, 3.205e-01, 2.450e-01, 1.894e-01, 1.462e-01, 1.146e-01, 8.950e-02, 6.994e-02, 5.557e-02, 4.399e-02,
				 3.461e-02, 2.783e-02, 2.219e-02, 1.795e-02, 1.416e-02, 1.137e-02, 9.077e-03, 7.411e-03, 5.976e-03, 4.825e-03,
				 3.845e-03, 3.131e-03, 2.545e-03, 2.051e-03, 1.667e-03, 1.354e-03, 1.096e-03, 8.852e-04, 7.141e-04, 5.826e-04,
				 4.664e-04])


#Taken same for all couplings as not much difference
eff_qstar_LID    = array('d',  [ 0.488302,  0.49915,   0.51,      0.52,     0.53,      0.54,      0.55,      0.56,      0.57,      0.5793, 
                                 0.58861,   0.59115,   0.5937,    0.59634,   0.599,     0.60195,   0.6049,    0.60745,   0.61,      0.5993,
                                 0.615333,  0.61636,   0.6174,    0.61845,   0.6195,    0.62065,   0.6218,    0.62295,   0.6241,    0.62511,
                                 0.626118,  0.62584,   0.625565,  0.625537,  0.62551,   0.62538,   0.62525,   0.625065,  0.62488,   0.62474,
                                 0.624607])


xsEff_f1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_LID)])
xsEff_f1p0_fb = array('d', [1000.0 * float(b) for b in xsEff_f1p0])

xsEff_f0p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p5, eff_qstar_LID)])
xsEff_f0p5_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p5])

xsEff_f0p1 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p1, eff_qstar_LID)])
xsEff_f0p1_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p1])

xsEff_f0p2 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p2, eff_qstar_LID)])
xsEff_f0p2_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p2])

xsEff_f0p3 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p3, eff_qstar_LID)])
xsEff_f0p3_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p3])

xsEff_f0p4 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p4, eff_qstar_LID)])
xsEff_f0p4_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p4])

xsEff_f0p6 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p6, eff_qstar_LID)])
xsEff_f0p6_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p6])

xsEff_f0p7 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p7, eff_qstar_LID)])
xsEff_f0p7_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p7])

xsEff_f0p8 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p8, eff_qstar_LID)])
xsEff_f0p8_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p8])

xsEff_f0p9 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p9, eff_qstar_LID)])
xsEff_f0p9_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p9])


graph_exp_2sigma = TGraph(len(masses_exp_tev),masses_exp_tev,xs_exp_limits_2sigma_fb)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times B #times A #times #epsilon [fb]")
graph_exp_2sigma.GetYaxis().SetRangeUser(1e-1,300)
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

graph_qstar_f1p0 = TGraph(len(masses_qstar),masses_qstar,xsEff_f1p0_fb)
graph_qstar_f1p0.SetLineWidth(2)                
graph_qstar_f1p0.SetLineColor(2)

graph_qstar_f0p5 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p5_fb)
graph_qstar_f0p5.SetLineWidth(2)                
graph_qstar_f0p5.SetLineColor(12)

graph_qstar_f0p1 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p1_fb)
graph_qstar_f0p1.SetLineWidth(2)                
graph_qstar_f0p1.SetLineColor(4)

graph_qstar_f0p2 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p2_fb)
graph_qstar_f0p2.SetLineWidth(2)                
graph_qstar_f0p2.SetLineColor(42)

graph_qstar_f0p3 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p3_fb)
graph_qstar_f0p3.SetLineWidth(2)                
graph_qstar_f0p3.SetLineColor(6)

graph_qstar_f0p4 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p4_fb)
graph_qstar_f0p4.SetLineWidth(2)                
graph_qstar_f0p4.SetLineColor(7)

graph_qstar_f0p6 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p6_fb)
graph_qstar_f0p6.SetLineWidth(2)                
graph_qstar_f0p6.SetLineColor(36)

graph_qstar_f0p7 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p7_fb)
graph_qstar_f0p7.SetLineWidth(2)                
graph_qstar_f0p7.SetLineColor(9)

graph_qstar_f0p8 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p8_fb)
graph_qstar_f0p8.SetLineWidth(2)                
graph_qstar_f0p8.SetLineColor(32)

graph_qstar_f0p9 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p9_fb)
graph_qstar_f0p9.SetLineWidth(2)                
graph_qstar_f0p9.SetLineColor(46)


c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
graph_qstar_f1p0.Draw("L")
graph_qstar_f0p5.Draw("L")
graph_qstar_f0p1.Draw("L")
graph_qstar_f0p2.Draw("L")
graph_qstar_f0p3.Draw("L")
graph_qstar_f0p4.Draw("L")
graph_qstar_f0p6.Draw("L")
graph_qstar_f0p7.Draw("L")
graph_qstar_f0p8.Draw("L")
graph_qstar_f0p9.Draw("L")


legend = TLegend(.65,.62,.90,.80)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.023)
legend.SetHeader('95% CL Upper Limits')
legend.AddEntry(graph_obs,"Observed limit","lp")
legend.AddEntry(graph_exp,"Expected limit","lp")
legend.AddEntry(graph_exp_1sigma,"Expected limit #pm 1#sigma","f")
legend.AddEntry(graph_exp_2sigma,"Expected limit #pm 2#sigma","f")
legend.Draw()

legend1 = TLegend(.17,.17,.42,.50)
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetFillStyle(0)
legend1.SetTextFont(42)
legend1.SetTextSize(0.025)
legend1.SetHeader('Diff q* Couplings')
legend1.AddEntry(graph_qstar_f0p1,"f = 0.1","l")
legend1.AddEntry(graph_qstar_f0p2,"f = 0.2","l")
legend1.AddEntry(graph_qstar_f0p3,"f = 0.3","l")
legend1.AddEntry(graph_qstar_f0p4,"f = 0.4","l")
legend1.AddEntry(graph_qstar_f0p5,"f = 0.5","l")
legend1.AddEntry(graph_qstar_f0p6,"f = 0.6","l")
legend1.AddEntry(graph_qstar_f0p7,"f = 0.7","l")
legend1.AddEntry(graph_qstar_f0p8,"f = 0.8","l")
legend1.AddEntry(graph_qstar_f0p9,"f = 0.9","l")
legend1.AddEntry(graph_qstar_f1p0,"f = 1.0","l")
legend1.Draw()


l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.025)
l1.DrawLatex(0.70,0.90, "CMS Preliminary")
#l1.DrawLatex(0.215,0.37, "CMS (unpublished)")
l1.DrawLatex(0.63,0.86, "#sqrt{s} = 13 TeV")
l1.DrawLatex(0.76,0.86, "#int Ldt = 2.5 fb^{-1}")

gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_withAllXS.pdf')
c.SaveAs('ExcitedbQuarksToGbJ_f1p0_withAllXS.png')

#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ExpVsQstar_zoomed.pdf')
#c.SaveAs('ExcitedbQuarksToGbJ_f1p0_DEtaOptimization_ExpVsQstar_zoomed.png')



