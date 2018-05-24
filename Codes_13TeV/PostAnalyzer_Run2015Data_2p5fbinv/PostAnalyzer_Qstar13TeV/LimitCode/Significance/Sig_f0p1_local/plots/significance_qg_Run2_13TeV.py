#!/usr/bin/env python

from ROOT import gROOT, gStyle, gPad, TF1, TH1F, TCanvas, TArrow, TGraph, kWhite, kRed, kBlue
from array import array
import CMS_lumi


gROOT.SetBatch(1)
gStyle.SetOptStat(111111)
gStyle.SetOptFit(1111)
#gStyle.SetOptTitle(0)
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
significances = array('d')

mass_min = 1500
mass_max = 7200

##------------------------------------------------------
## for reading the limit code log files
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(float(mass))

  #sig = TH1F("sig","Significance: qg m=" + str(int(mass)) + " GeV;Sig = sgn(S)#sqrt{-2ln(L_{B}/L_{S+B})};Entries/bin", 40, -4., 4.)

  #log_file = open("stats_" + str(int(mass)) + "_qg.log",'r')
  #outputlines = log_file.readlines()
  #log_file.close()

  #sig_obs = 0.

  #for line in outputlines:
    #if "Significance" in line:
      #if "(data)" in line:
        ##print line.split()[-1]
        #sig_obs = float(line.split()[-1])
      #else:
        ##print line.split()[-1]
        #sig.Fill( float(line.split()[-1]) )

  #if sig_obs<0.:
    #significances.append(0.)
  #else:
    #significances.append(sig_obs)


  ##c = TCanvas("c", "",800,800)
  ##c.cd()

  ##sig.Draw("hist")

  ##func = TF1("func", "gaus", -4., 4)
  ##sig.Fit("func","R")
  ##func.Draw("same")

  ##arrow = TArrow(sig_obs,0.3*sig.GetMaximum(),sig_obs,0.,0.05,"|>");
  ##arrow.SetLineWidth(2)
  ##arrow.Draw()

  ##gPad.Update()

  ##sig.SaveAs('significance_' + str(int(mass)) + '_qg.root')
  ##c.SaveAs('significance_' + str(int(mass)) + '_qg.eps')

#print "masses =", masses
#print "significances =", significances

##------------------------------------------------------

masses = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0])
significances = array('d', [0.0, 0.361683, 1.17934, 0.937784, 0.20214, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.351713, 0.54517, 0.698591, 0.877327, 0.877771, 1.01375, 1.27896, 1.35268, 1.12581, 1.1746, 1.46086, 1.41497, 0.872541, 0.00435573, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0856762, 0.433186, 0.542791, 0.55899, 0.556067, 0.573986, 0.583046, 0.538931, 0.459697, 0.332036, 0.152038, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

##------------------------------------------------------

graph_sig = TGraph(len(masses),masses,significances)
graph_sig.GetXaxis().SetTitle("qg resonance mass [GeV]")
graph_sig.GetYaxis().SetTitle("Significance (local)")
graph_sig.GetYaxis().SetTitleOffset(1.2)
graph_sig.GetYaxis().SetRangeUser(0.,3.0)
graph_sig.SetLineWidth(2)
graph_sig.SetLineColor(kRed)
graph_sig.SetMarkerStyle(21)
graph_sig.SetMarkerSize(1)
graph_sig.SetMarkerColor(kBlue)
#graph_sig.GetXaxis().SetNdivisions(1005)


gStyle.SetOptTitle(0)
c = TCanvas("c", "",800,800)
c.cd()

graph_sig.Draw("ALP")

# draw the lumi text on the canvas
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "2.4 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

CMS_lumi.CMS_lumi(c, iPeriod, iPos)

gPad.RedrawAxis()

c.SetGridx()
c.SetGridy()
c.SaveAs('significance_DijetLimitCode_qg_Run2_13TeV_DATA_2445_invpb.eps')
