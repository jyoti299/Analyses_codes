#!/usr/bin/env python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi

CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "1 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

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
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()

#######################################################
def getQuantiles(limits, median, onesigma, twosigma):

    nit = len(limits)
    if nit==0:
        return

    # sort the list with limits
    limits = sorted(limits)
    limits_array = array('d',limits)

    # median for the expected limit
    median[0] = TMath.Median(nit, limits_array)

    # quantiles for the expected limit bands
    prob = array('d',[]) # array with quantile boundaries
    prob.append(0.021)
    prob.append(0.159)
    prob.append(0.841)
    prob.append(0.979)

    # array for the results

    quantiles = array('d',[0., 0., 0., 0.])
    TMath.Quantiles(nit, 4, limits_array, quantiles, prob) # evaluate quantiles

    onesigma[0] = quantiles[1]
    onesigma[1] = quantiles[2]
    twosigma[0] = quantiles[0]
    twosigma[1] = quantiles[3]

    return
#######################################################

masses = array('d')
sig = array('d')
masses_band = array('d')
sig_1sigma = array('d')
sig_1sigma_up = array('d')
sig_2sigma = array('d')
sig_2sigma_up = array('d')

mass_min = 1200
mass_max = 6000
ntoys = 200

########################################################
## Uncomment this part if running the limit code

### for running the limit code
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(mass)
  #masses_band.append(mass)

  #significances = []

  #for toy in range(1,ntoys+1):
    #cmd = "./stats " + str(int(mass)) + " gg 0 0 " + str(toy)
    #print "Running: " + cmd
    #proc = subprocess.Popen( cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
    #output = proc.communicate()[0]
    #if proc.returncode != 0:
      #print output
      #sys.exit(1)
    ##print output

    #outputlines = output.split("\n")

    #for line in outputlines:
      #if "Significance(data)" in line:
        ##print line
        #significances.append(float(line.split()[-1]))

  #median = [0.]
  #onesigma = [0., 0.]
  #twosigma = [0., 0.]

  #getQuantiles(significances, median, onesigma, twosigma)

  #sig.append(float(median[0]))
  #sig_1sigma.append(float(onesigma[0]))
  #sig_1sigma_up.append(float(onesigma[1]))
  #sig_2sigma.append(float(twosigma[0]))
  #sig_2sigma_up.append(float(twosigma[1]))
##------------------------------------------------------

#for i in range(0,len(masses)):
  #masses_band.append( masses[len(masses)-i-1] )
  #sig_1sigma.append( sig_1sigma_up[len(masses)-i-1] )
  #sig_2sigma.append( sig_2sigma_up[len(masses)-i-1] )


#print "masses:"
#print masses
#print "sig:"
#print sig
#print ""
#print "masses_band:"
#print masses_band
#print "sig_1sigma:"
#print sig_1sigma
#print "sig_2sigma:"
#print sig_2sigma

##
########################################################

########################################################
## Comment out this part if running the limit code

masses = array('d', [1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0])
sig = array('d', [0.18072749999999999, -0.46114, -0.29333050000000005, -0.0113902, 0.304138, 0.407108, 0.5327124999999999, 0.5266095, 0.6315459999999999, 0.482504, 0.3272425, 0.1363095, 0.0950312, -0.17319449999999997, -0.3506555, -0.48411150000000003, -0.513308, -0.7109315, -0.865345, -1.053785, -1.0634700000000001, -1.017605, -0.8916715, -0.9372905, -0.8120205, -0.647459, -0.534994, -0.424716, -0.2656155, -0.08349545, 0.2381985, 0.494488, 0.9146719999999999, 1.27687, 1.7425950000000001, 2.296575, 2.742965, 3.1292850000000003, 3.318485, 3.172495, 2.8363750000000003, 2.4979050000000003, 2.06827, 1.529305, 1.1226150000000001, 0.7245325, 0.35265599999999997, 0.0414658, -0.247311])

masses_band = array('d', [1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6000.0, 5900.0, 5800.0, 5700.0, 5600.0, 5500.0, 5400.0, 5300.0, 5200.0, 5100.0, 5000.0, 4900.0, 4800.0, 4700.0, 4600.0, 4500.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0])
sig_1sigma = array('d', [-0.76097942, -1.36975989, -1.2239281800000001, -1.01351101, -0.6282624139999999, -0.601404027, -0.698770528, -0.572673864, -0.6929832829999999, -0.716412018, -0.724918003, -0.838146677, -0.932035285, -1.18362188, -1.16188223, -1.3812459899999996, -1.6270255799999997, -1.7506956099999997, -1.9366155699999998, -1.9306860599999998, -1.9720846999999995, -2.02606516, -2.08757697, -2.1698528399999994, -2.04628917, -1.8209224699999997, -1.57372724, -1.4212375499999998, -1.24780596, -1.00898308, -0.6901418309999998, -0.4284106949999999, -0.022007613519999636, 0.359664458, 0.7802130720000005, 1.3143696500000002, 1.72611131, 1.99032204, 2.09690856, 2.05720001, 1.76288815, 1.4494674600000002, 0.9798569220000002, 0.5275938810000002, 0.01709345220000031, -0.46774341499999983, -0.78689666, -1.0342562699999998, -1.2941329199999998, 0.781649388, 1.0796862199999984, 1.347869319999998, 1.6765876299999998, 1.9844193099999998, 2.5278112200000002, 3.06814911, 3.615168259999999, 4.0366792799999995, 4.241194649999999, 4.26093078, 4.02628481, 3.7263014399999994, 3.3482520900000003, 2.7889026699999997, 2.3426443499999996, 1.8308230999999997, 1.3335363, 1.07606862, 0.766365777, 0.6181976269999998, 0.47798655199999995, 0.3477848769999994, 0.12284800999999931, 0.129262266, 0.17583514999999933, 0.07570354059999938, -0.05476961520000026, -0.0788715501, 0.0394586168999989, 0.12559076099999975, 0.33103544299999993, 0.30476252699999984, 0.5607433179999999, 0.5834689939999996, 0.7703687399999992, 0.8270082569999998, 1.14222531, 1.3949328199999989, 1.4330281099999995, 1.4443118499999996, 1.4948288599999997, 1.3535776299999998, 1.4995090999999992, 1.1808759099999997, 0.8811262609999996, 0.5139568749999992, 0.6000548939999999, 1.1942229599999994])
sig_2sigma = array('d', [-1.84565345, -2.20397083, -2.32323244, -2.2123769199999996, -1.89687996, -1.76684659, -1.8644351399999999, -1.68239322, -1.57142886, -1.3508323, -1.71673706, -1.96892182, -2.35504632, -2.1990581999999996, -2.20739983, -2.48637262, -2.95059642, -3.3106218, -3.0866486, -3.01376419, -2.9551583, -2.9517799399999998, -3.23732754, -3.2908139500000004, -3.23849761, -2.9528299099999997, -2.8528085699999997, -2.27190045, -2.06916584, -2.00679343, -1.82357654, -1.48918104, -0.9153237909999999, -0.6990133789999999, -0.166460434, 0.381097164, 0.951038327, 1.23628846, 1.36052391, 1.26600394, 0.9557143379999999, 0.3442203850000001, -0.140761124, -0.486536176, -0.9569343770000001, -1.19426209, -1.7078748, -1.88936471, -2.566568, 1.63703776, 2.06938042, 2.45418514, 2.7483304299999998, 3.1710082799999997, 3.6401182999999993, 4.03757705, 4.4088167899999995, 5.011359479999999, 5.3332393, 5.37586551, 4.9870987499999995, 4.54228996, 4.03004521, 3.58185027, 3.3681134999999998, 2.84526594, 2.33755307, 2.1395874699999995, 1.6236116299999999, 1.3919140299999997, 1.2188086399999998, 1.3872082499999996, 1.4532013199999994, 1.2380241519999995, 0.8360614289999999, 0.9037561629999996, 0.62725538, 0.9815045769999997, 1.0917246249999997, 0.9821873869999997, 1.1438587999999998, 1.3351659299999998, 1.54940241, 1.90825551, 2.04756795, 1.8659728099999997, 1.8455242399999998, 2.3768013299999997, 2.33516308, 2.2492541899999994, 2.4220601999999993, 2.46774863, 2.35240541, 2.3451681499999992, 1.8019056399999998, 1.3680189999999999, 1.2702152, 2.21760668])

##
########################################################

graph_sig_2sigma = TGraph(len(masses_band),masses_band,sig_2sigma)
graph_sig_2sigma.SetFillColor(kYellow)
graph_sig_2sigma.GetXaxis().SetTitle("gg resonance mass [GeV]")
graph_sig_2sigma.GetYaxis().SetTitle("Significance")
graph_sig_2sigma.GetYaxis().SetRangeUser(-3.,5.)
#graph_sig_2sigma.GetXaxis().SetNdivisions(1005)

graph_sig_1sigma = TGraph(len(masses_band),masses_band,sig_1sigma)
graph_sig_1sigma.SetFillColor(kGreen+1)

graph_sig = TGraph(len(masses),masses,sig)
#graph_sig.SetMarkerStyle(24)
graph_sig.SetLineWidth(2)
graph_sig.SetLineStyle(2)
graph_sig.SetLineColor(4)

c = TCanvas("c", "",800,800)
c.cd()

graph_sig_2sigma.Draw("AF")
graph_sig_1sigma.Draw("F")
graph_sig.Draw("L")

legend = TLegend(.20,.75,.50,.80)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
legend.AddEntry(graph_sig,"Median significance","lp")
legend.Draw()

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c, iPeriod, iPos)

gPad.RedrawAxis()

#c.SetLogy()
c.SaveAs('significance_blindExercise_gg_Run2_13TeV.pdf')
