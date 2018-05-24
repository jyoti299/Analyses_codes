#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.045, "XYZ") #0.06
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.045, "XYZ") #0.06
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
t_m = 0.06  ##top margin 0.55
b_m = 0.11   ##botton margin 0.13
l_m = 0.12  ##left margin 0.13
r_m = 0.12  ##right margin 0.06
gStyle.SetPadTopMargin(t_m) ##0.06
gStyle.SetPadBottomMargin(b_m) ##0.15
gStyle.SetPadLeftMargin(l_m)  ###0.15
gStyle.SetPadRightMargin(r_m)  ###0.05
gROOT.ForceStyle()

### 13 TeV
Coupling_Obs = array('d', [ ])

ExclMass_Obs = array('d', [ ])

grCoupling_Obs = TGraph(len(Coupling_Obs), ExclMass_Obs, Coupling_Obs)
grCoupling_Obs.SetFillColor(kYellow-4)
grCoupling_Obs.SetLineColor(kGreen-6)
grCoupling_Obs.SetLineWidth(2)

grCoupling_Obs.GetYaxis().SetTitle("Couplings [f = f_{s} = f']")
grCoupling_Obs.GetYaxis().CenterTitle()
grCoupling_Obs.GetYaxis().SetTitleSize(0.05)
grCoupling_Obs.GetYaxis().SetTitleOffset(1.1) ##1.1
grCoupling_Obs.GetYaxis().SetLabelSize(0.045) ##0.05
grCoupling_Obs.GetYaxis().SetLabelOffset(0.01) 
grCoupling_Obs.GetYaxis().SetRangeUser(0.0,1.0)


RightYaxis = TGaxis( 6., 0.0, 6., 1.0,0,1,510,"+L")
RightYaxis.SetLabelFont(42)
RightYaxis.SetTitle("M_{q*} / #Lambda")
RightYaxis.SetTitleFont(42)
RightYaxis.SetTitleSize(0.05)
RightYaxis.SetTitleOffset(1.1)
RightYaxis.CenterTitle()

grCoupling_Obs.GetXaxis().SetTitle("q* Mass [TeV]")
grCoupling_Obs.GetXaxis().CenterTitle()
grCoupling_Obs.GetXaxis().SetTitleSize(0.05)
grCoupling_Obs.GetXaxis().SetTitleOffset(1.1)
grCoupling_Obs.GetXaxis().SetLabelSize(0.04) ##0.05
grCoupling_Obs.GetXaxis().SetLabelOffset(0.008)
grCoupling_Obs.GetXaxis().SetLimits(1.,6.)
grCoupling_Obs.GetXaxis().SetNdivisions(510)


grCoupling_Exp = TGraph(len(Coupling_Exp), ExclMass_Exp, Coupling_Exp)
grCoupling_Exp.SetLineColor(kBlue+3)
grCoupling_Exp.SetLineStyle(9)
grCoupling_Exp.SetLineWidth(2)

c = TCanvas("c", "",800,800)
c.cd()

grCoupling_Obs.Draw("AFL")
RightYaxis.Draw()
    
legend = TLegend(0.36,0.23,0.7,0.31)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)  ##35
legend.SetHeader('Observed 95% CL Exclusion')
legend.AddEntry(grCoupling_Obs, "#sqrt{s} = 13 TeV", "F")
legend.Draw()

lumiTextSize = 0.6
lumiTextOffset = 0.2
lumi = TLatex()
lumi.SetNDC()
lumi.SetTextAngle(0)
lumi.SetTextColor(kBlack)
lumi.SetTextFont(42)
lumi.SetTextAlign(31)
lumi.SetTextSize(lumiTextSize*t_m)
lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "2.67 fb^{-1} (13 TeV)")

cmsTextFont = 61
cmsTextSize = 0.75
##posX_ = l_m     + 0.045 * (1 - l_m - r_m)  ##Top left
##posX_ = l_m     + 0.5   * (1 - l_m - r_m)   ## Centre
posX_ = 1 - r_m + 0.045 * (1 - l_m - r_m)  ## Top right
posY_ = 1 - t_m - 0.035 * (1 - t_m - b_m)

cms =  TLatex()
cms.SetNDC()
cms.SetTextFont(cmsTextFont)
cms.SetTextSize(cmsTextSize * t_m)
cms.SetTextAlign(33)  ### 11-top left;  21-top centre;  31-top right
##cms.DrawLatex(posX_, posY_, "CMS")
###cms.DrawLatex(0.845, 0.88, "CMS")
cms.DrawLatex(0.82, 0.88, "CMS")
## if additional material need to be added: Uncomment below:
cms.SetTextFont(52);
cms.SetTextAlign(33);             
cms.SetTextSize(0.76*cmsTextSize*t_m);
##cms.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t_m, extraText);
cms.DrawLatex(0.825, 0.88 - 1.2*cmsTextSize*t_m, "Preliminary");


l1 = TLatex()
l1.SetNDC()
l1.SetTextAlign(13)
l1.SetTextSize(0.04) ##0.035
l1.SetTextFont(42)
##l1.DrawLatex(0.735,0.79, "q*#rightarrow q#gamma")
l1.DrawLatex(0.7,0.76, "q*#rightarrow q#gamma")

gPad.RedrawAxis();
c.SaveAs('Limit_.pdf')
