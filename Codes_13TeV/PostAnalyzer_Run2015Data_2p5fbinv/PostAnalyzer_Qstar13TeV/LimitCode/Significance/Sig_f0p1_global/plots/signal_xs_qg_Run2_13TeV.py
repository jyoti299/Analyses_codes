#!/usr/bin/env python

import sys, os, subprocess, string, re
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


masses = array('d')
sig = array('d')
sig_ex = array('d')
sig_ey = array('d')

mass_min = 1500
mass_max = 7200

########################################################
## Uncomment this part if running the limit code

### for running the limit code
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(mass)

  #cmd = "./stats " + str(int(mass)) + " qg"
  #print "Running: " + cmd
  #proc = subprocess.Popen( cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
  #output = proc.communicate()[0]
  #if proc.returncode != 0:
    #print output
    #sys.exit(1)
  ##print output

  #outputlines = output.split("\n")

  #foundLine = False

  #for line in outputlines:
    #if "Fitted signal xs:" in line:
      ##print line
      #foundLine = True
      #sig.append(float(line.split()[3]))
      #sig_ey.append(float(line.split()[5]))
      #break

  #if not foundLine:
    #print "Signal fit failed for m="+str(mass)+" GeV"
    #sig.append(0.)
    #sig_ey.append(0.)

  #sig_ex.append(0.)

#------------------------------------------------------

### for reading the limit code log files
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(mass)

  #log_file = open("stats_" + str(int(mass)) + "_qg.log",'r')
  #outputlines = log_file.readlines()
  #log_file.close()

  #foundLine = False

  #for line in outputlines:
    #if "Fitted signal xs:" in line:
      ##print line
      #foundLine = True
      #sig.append(float(line.split()[3]))
      #sig_ey.append(float(line.split()[5]))
      #break

  #if not foundLine:
    #print "Signal fit failed for m="+str(mass)+" GeV"
    #sig.append(0.)
    #sig_ey.append(0.)

  #sig_ex.append(0.)

##------------------------------------------------------

#print "masses =", masses
#print "sig =", sig
#print "sig_ex =", sig_ex
#print "sig_ey =", sig_ey

##
########################################################

########################################################
## Comment out this part if running the limit code

masses = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0])
sig = array('d', [-0.021354, 0.181856, 0.507565, 0.301059, 0.0481507, -0.0745822, -0.173741, -0.219941, -0.182699, -0.0981657, -0.0838652, -0.0989725, -0.0333594, -0.000836467, 0.0192645, 0.0267208, 0.0304998, 0.0335913, 0.0298443, 0.0290469, 0.0354452, 0.0337065, 0.0255528, 0.0241572, 0.0264757, 0.0228566, 0.0128836, 3.55709e-05, -0.0101949, -0.0160754, -0.0174712, -0.0164995, -0.0154132, -0.0126426, -0.00746283, -0.0026017, 0.000291682, 0.0013748, 0.00231994, 0.00222654, 0.0020308, 0.00189199, 0.00173025, 0.00117754, 0.000933006, 0.00062805, 0.000268176, -7.43922e-05, -0.000390788, -0.000737904, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
sig_ex = array('d', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
sig_ey = array('d', [0.501878, 0.390429, 0.32588, 0.271627, 0.226779, 0.191722, 0.16158, 0.136454, 0.115358, 0.0977854, 0.0837166, 0.0716966, 0.0618743, 0.0535481, 0.0466339, 0.0407558, 0.0436593, 0.0316492, 0.0340568, 0.024885, 0.0221832, 0.0198079, 0.0211345, 0.0159719, 0.0142135, 0.0143907, 0.0115227, 0.011568, 0.00909428, 0.00889743, 0.00700818, 0.0067786, 0.00559803, 0.0051141, 0.00468355, 0.00461555, 0.00409716, 0.00385913, 0.0036742, 0.00345634, 0.00321343, 0.00304897, 0.00271898, 0.00247725, 0.00228723, 0.0021062, 0.0019325, 0.00176778, 0.0015361, 0.00136623, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

##
########################################################

sig_pos = array('d')
sig_exl = array('d')
sig_exh = array('d')
sig_eyl = array('d')
sig_eyh = array('d')

# create final arrays
for i in range(0,len(masses)):
  sig_pos.append(sig[i] if sig[i]>0. else 0.)
  sig_exl.append(sig_ex[i])
  sig_exh.append(sig_ex[i])
  sig_eyl.append(sig_ey[i] if sig[i]>0. else 0.)
  sig_eyh.append(sig_ey[i])

graph_sig = TGraphAsymmErrors(len(masses),masses,sig_pos,sig_exl,sig_exh,sig_eyl,sig_eyh)
graph_sig.GetXaxis().SetTitle("qg resonance mass [GeV]")
graph_sig.GetYaxis().SetTitle("Signal cross section [pb]")
graph_sig.GetYaxis().SetTitleOffset(1.2)
graph_sig.GetYaxis().SetRangeUser(1e-4,1e1)
graph_sig.SetMarkerStyle(20)
graph_sig.SetMarkerColor(1)
graph_sig.SetLineWidth(2)
graph_sig.SetLineStyle(1)
graph_sig.SetLineColor(1)

c = TCanvas("c", "",800,800)
c.cd()

graph_sig.Draw("AP")

#draw the lumi text on the canvas
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "2.4 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

CMS_lumi.CMS_lumi(c, iPeriod, iPos)

gPad.RedrawAxis()

c.SetLogy()
c.SetGridx()
c.SetGridy()
c.SaveAs('signal_xs_DijetLimitCode_qg_Run2_13TeV_DATA_2445_invpb.eps')
