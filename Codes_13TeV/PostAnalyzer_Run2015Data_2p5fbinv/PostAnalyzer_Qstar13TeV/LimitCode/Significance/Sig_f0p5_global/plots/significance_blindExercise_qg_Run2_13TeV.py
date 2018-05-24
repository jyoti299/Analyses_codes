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
    #cmd = "./stats " + str(int(mass)) + " qg 0 0 " + str(toy)
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
sig = array('d', [0.03894795, -0.4974445, -0.3262405, 0.03041775, 0.3182955, 0.3918105, 0.4440465, 0.5293300000000001, 0.5717734999999999, 0.380628, 0.2026525, 0.163444, 0.03272695, -0.223557, -0.330528, -0.42696100000000003, -0.520624, -0.723382, -0.815175, -1.015215, -1.02704, -0.994483, -0.9192834999999999, -0.9013055, -0.756575, -0.664498, -0.5416955, -0.40694600000000003, -0.3170145, -0.028118499999999998, 0.26863, 0.6533340000000001, 0.9715605, 1.39898, 1.888825, 2.43004, 2.892695, 3.212575, 3.2760800000000003, 3.0345, 2.70158, 2.40736, 1.940305, 1.421975, 1.00839, 0.5847595, 0.19636, -0.0669791, -0.345885])

masses_band = array('d', [1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6000.0, 5900.0, 5800.0, 5700.0, 5600.0, 5500.0, 5400.0, 5300.0, 5200.0, 5100.0, 5000.0, 4900.0, 4800.0, 4700.0, 4600.0, 4500.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0])
sig_1sigma = array('d', [-0.9747270349999999, -1.45272447, -1.2861774599999998, -1.03513039, -0.6671459419999997, -0.6250242619999999, -0.6295515059999998, -0.587794465, -0.722877239, -0.656293995, -0.717535177, -0.7695662289999999, -0.974853715, -1.15592178, -1.10116749, -1.45737557, -1.5529647899999999, -1.7396211999999998, -1.9144879699999997, -1.9201495199999998, -1.9706828499999998, -2.01219467, -2.0125455699999995, -2.12647958, -1.97407829, -1.7355094899999999, -1.52218779, -1.4023154, -1.19708694, -0.916201214, -0.6256168629999999, -0.3003123229999999, 0.07983452940000005, 0.5195992950000002, 0.9340079960000003, 1.46934922, 1.8338081600000002, 2.02619395, 2.11224437, 1.97528408, 1.6635302900000002, 1.3338104300000002, 0.8328362830000005, 0.39234495900000044, -0.11651775599999997, -0.5557266879999999, -0.8749717779999999, -1.1192683999999997, -1.3494438699999995, 0.5996066449999998, 0.9180043319999993, 1.2011684599999997, 1.5461718699999991, 1.91371219, 2.3096979399999995, 2.9319648999999997, 3.50808346, 3.9056573500000002, 4.24087985, 4.244649429999999, 4.11623374, 3.8547463499999997, 3.4213590199999997, 2.9706807299999998, 2.5294801899999997, 1.98574816, 1.5079159699999987, 1.0858147999999994, 0.823789619, 0.6215958029999991, 0.5065011629999998, 0.3667503379999995, 0.2549052159999998, 0.1924620539999997, 0.2390784689999999, 0.14706755499999963, -0.000673562600000596, -0.010612612640000133, 0.09244252189999977, 0.05933940319999887, 0.2627950849999996, 0.3682857229999999, 0.5829611099999997, 0.5845461279999998, 0.7806337289999996, 0.826499212, 1.0768465799999993, 1.3455141499999996, 1.2949264899999995, 1.4867566399999992, 1.4532506499999998, 1.29192391, 1.4091285099999997, 1.0741193199999977, 0.8794341969999999, 0.587432063, 0.489517164, 1.06419798])
sig_2sigma = array('d', [-1.99442542, -2.27262251, -2.7013034599999997, -2.2677496699999997, -1.88097366, -1.68478421, -2.01563415, -1.6790445100000002, -1.48312625, -1.5531232499999998, -1.85584379, -1.96400613, -2.3932799, -2.1283087700000003, -2.3091856799999997, -2.6829483599999997, -3.07525446, -3.17758734, -3.1452932899999997, -2.8214307400000003, -2.84690783, -2.92094447, -3.1021299800000004, -3.23927776, -3.2375126499999998, -2.93816214, -2.64220588, -2.14409218, -2.17546664, -2.03498752, -1.7462971999999999, -1.42327738, -0.8741884999999999, -0.5894312239999999, 0.016863520820000023, 0.600511914, 1.0541292800000002, 1.3516112900000001, 1.37782957, 1.18344477, 0.831544177, 0.237665807, -0.295606148, -0.617307793, -0.992824312, -1.3938893499999998, -1.8563015999999999, -1.95370727, -2.3503296799999998, 1.4525638299999994, 1.9037144199999998, 2.27995869, 2.6927685099999996, 3.0878431599999994, 3.5564725399999997, 3.91355282, 4.26127977, 4.81333407, 5.28240078, 5.3943886800000005, 5.15086164, 4.7482356, 4.2068072, 3.7109125199999995, 3.40999129, 2.97247297, 2.4502246199999997, 2.18798103, 1.8465763199999994, 1.4141895199999999, 1.1820524399999999, 1.2839520999999998, 1.4518471099999999, 1.332888697999999, 0.9479988069999997, 0.888269508, 0.7714936549999997, 0.8150547109999998, 1.1094358209999995, 1.13415935, 1.1596882899999998, 1.30283223, 1.57061831, 1.8905948600000002, 2.10792361, 1.9532274399999998, 1.9329075500000001, 2.06198408, 2.44062165, 2.31572666, 2.28677083, 2.3172663499999997, 2.36177528, 2.3498348499999997, 1.7129987599999998, 1.35407588, 1.28299098, 2.31814421])

##
########################################################

graph_sig_2sigma = TGraph(len(masses_band),masses_band,sig_2sigma)
graph_sig_2sigma.SetFillColor(kYellow)
graph_sig_2sigma.GetXaxis().SetTitle("qg resonance mass [GeV]")
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
c.SaveAs('significance_blindExercise_qg_Run2_13TeV.pdf')
