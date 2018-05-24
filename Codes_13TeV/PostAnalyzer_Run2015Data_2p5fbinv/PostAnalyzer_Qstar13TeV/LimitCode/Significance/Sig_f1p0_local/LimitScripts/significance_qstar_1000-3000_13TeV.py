#!/usr/bin/env python

from ROOT import gROOT, gStyle, gPad, TF1, TH1F, TCanvas, TArrow, TGraph, TGaxis, kWhite, kRed, kBlue, kMagenta, kGreen
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
gStyle.SetPadRightMargin(0.07)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadBottomMargin(0.14)
gROOT.ForceStyle()


masses = array('d')
significances = array('d')

mass_min = 1000
mass_max = 3000

##------------------------------------------------------
## for reading the limit code log files
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(float(mass))

  #sig = TH1F("sig","Significance: qq m=" + str(int(mass)) + " GeV;Sig = sgn(S)#sqrt{-2ln(L_{B}/L_{S+B})};Entries/bin", 40, -4., 4.)

  #log_file = open("stats_" + str(int(mass)) + "_qq.log",'r')
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

  ##sig.SaveAs('significance_' + str(int(mass)) + '_qq.root')
  ##c.SaveAs('significance_' + str(int(mass)) + '_qq.eps')

#print "masses =", masses
#print "significances =", significances

##------------------------------------------------------

masses = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0])
#Got from the code run
#significances = array('d', [-0.0910847, -0.620492, -0.182553, -1.27005, -1.13241, -0.50456, -0.681315, -1.31427, 0.0395607, 2.35209, 3.23294, 0.662911, -0.437585, 0.945798, 1.08131, 0.242099, -1.13291, 0.0, -0.819911, -0.638817, 0.0]) 
#Making all the negative values to zero
#74X
#significances = array('d', [0, 0, 0, 0, 0, 0.923534, 0, 0, 0.827628, 2.4343, 2.10423, 0.447207, 0.668827, 1.13151, 0.0508347, 0, 0.0, 0, 0, 0, 0])
#76X
significances = array('d', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0395607, 2.35209, 3.23294, 0.662911, 0.0, 0.945798, 1.08131, 0.242099, 0.0, 0.0, 0.0, 0.0, 0.0])

##------------------------------------------------------

graph_sig = TGraph(len(masses),masses,significances)
graph_sig.GetXaxis().SetTitle("q* resonance mass [GeV]")
graph_sig.GetYaxis().SetTitle("Significance (local)")
graph_sig.GetYaxis().SetTitleOffset(1.2)
graph_sig.GetYaxis().SetRangeUser(0.,4.0)
graph_sig.GetYaxis().SetLabelFont(42)
graph_sig.GetYaxis().SetLabelSize(0.04)
graph_sig.GetXaxis().SetLabelSize(0.04)
graph_sig.GetYaxis().SetTitleFont(42)
graph_sig.GetXaxis().SetTitleFont(42)
graph_sig.GetYaxis().SetTitleSize(0.05)
graph_sig.GetXaxis().SetTitleSize(0.05)
graph_sig.SetLineWidth(2)
graph_sig.SetLineColor(kRed)
graph_sig.SetMarkerStyle(29)
graph_sig.SetMarkerSize(1.6)
graph_sig.SetMarkerColor(kBlue)
#graph_sig.GetXaxis().SetNdivisions(1005)

RightYaxis = TGaxis(3200.0, 0.0, 3200.0, 4.0, 0, 4, 510, "+L")
RightYaxis.SetLabelFont(42)
RightYaxis.SetLabelSize(0.04)
RightYaxis.SetTitle("")
RightYaxis.SetTitleFont(42)
RightYaxis.SetTitleSize(0.05)
RightYaxis.SetTitleOffset(1.1)
RightYaxis.CenterTitle()


gStyle.SetOptTitle(0)
c = TCanvas("c", "",800,800)
c.cd()

graph_sig.Draw("ALP")
RightYaxis.Draw()
# draw the lumi text on the canvas
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "2.67 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

CMS_lumi.CMS_lumi(c, iPeriod, iPos)

gPad.RedrawAxis()

c.SetGridx()
c.SetGridy()
c.SaveAs('significance_qstar_1000-3000_13TeV_DATA_2502_invpb_f1p0_76X.pdf')
c.SaveAs('significance_qstar_1000-3000_13TeV_DATA_2502_invpb_f1p0_76X.png')
