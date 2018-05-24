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
    #cmd = "./stats " + str(int(mass)) + " qq 0 0 " + str(toy)
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
sig = array('d', [-0.1651715, -0.4792245, -0.17837750000000002, 0.06661995, 0.28867, 0.3844475, 0.3461785, 0.454203, 0.515432, 0.30227499999999996, 0.1950425, 0.20323049999999998, 0.07379785, -0.0652005, -0.247884, -0.3321225, -0.4129525, -0.6089370000000001, -0.7618419999999999, -0.8723905000000001, -0.9153125, -0.921064, -0.8867970000000001, -0.9390339999999999, -0.7908435, -0.7357484999999999, -0.6133315, -0.5397235, -0.412436, -0.194977, 0.1718675, 0.505324, 0.8386695, 1.3375400000000002, 2.02533, 2.515335, 2.9987899999999996, 3.24848, 3.1581900000000003, 2.819135, 2.482455, 1.97229, 1.4534799999999999, 1.06561, 0.6333085, 0.2950315, 0.01331883, -0.10046859999999999, -0.248907])

masses_band = array('d', [1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6000.0, 5900.0, 5800.0, 5700.0, 5600.0, 5500.0, 5400.0, 5300.0, 5200.0, 5100.0, 5000.0, 4900.0, 4800.0, 4700.0, 4600.0, 4500.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0])
sig_1sigma = array('d', [-1.23895322, -1.36275123, -1.2943618399999999, -1.0832120599999997, -0.7708680449999998, -0.7124653160000001, -0.6988050719999999, -0.6855135229999998, -0.763364135, -0.6648585819999999, -0.7218913659999995, -0.8292521660000001, -0.9150543889999999, -1.01722746, -1.07049343, -1.35019104, -1.44516351, -1.6824491199999998, -1.67634324, -1.7783254199999998, -1.85524762, -1.9281565299999999, -2.06861698, -2.13980815, -1.9295129599999998, -1.7567801799999998, -1.5885208299999998, -1.4947843199999997, -1.35246695, -1.10413636, -0.7942110219999999, -0.43915807899999976, -0.044558140099999935, 0.30913716800000035, 0.9137247130000001, 1.55828122, 1.8355388000000001, 2.00305554, 2.02783931, 1.72528844, 1.39597561, 0.854264412, 0.36215211700000005, -0.2465726139999996, -0.6197072189999999, -0.8231409970000001, -0.9450105019999999, -1.13290648, -1.21135847, 0.3844071619999998, 0.7114934469999998, 1.05083062, 1.2882724699999986, 1.5911122299999998, 1.9141064700000001, 2.40603618, 3.00320582, 3.488922259999999, 4.035114049999999, 4.272492389999999, 4.2726589399999995, 3.8771765499999997, 3.4892128799999993, 3.0167144199999996, 2.532785089999999, 1.9546975, 1.4265213099999996, 1.0135963199999996, 0.7736131819999997, 0.5462755679999995, 0.3639968269999993, 0.287216451999999, 0.1758325659999992, 0.004584781899999162, 0.15124982099999978, 0.23107593099999946, 0.06477727129999944, 0.04976670569999997, 0.13906881499999943, 0.2068792959999997, 0.455983413, 0.5491599599999999, 0.6652721999999991, 0.7747180209999989, 0.7943101889999999, 0.9347047909999993, 1.0313074799999997, 1.2906043699999992, 1.3570487699999998, 1.4038006899999989, 1.3954761299999998, 1.2433902899999993, 1.3682077000000001, 1.1770228999999999, 0.9773947169999997, 0.5754074599999998, 0.5178327959999995, 1.007572230999999])
sig_2sigma = array('d', [-1.9797867400000002, -2.09364429, -2.2077735699999996, -1.9860571699999998, -1.9356898500000002, -2.03997674, -1.7291949999999998, -1.73432485, -1.6945088599999998, -1.75618863, -1.9334036499999998, -2.0745534500000002, -2.01790446, -1.92341463, -2.24425164, -2.66190312, -3.12101978, -2.99744251, -2.576704, -2.68759962, -2.6822189, -3.0107950199999998, -3.2780712800000003, -3.24919995, -3.17671167, -3.03866929, -2.39371202, -2.38525529, -2.40368451, -2.15481436, -2.00024139, -1.43366131, -1.1235906629999999, -0.527513687, 0.10426035740000002, 0.6838297059999999, 1.01088977, 1.31283198, 1.20971238, 0.8269548290000001, 0.327868882, -0.289212988, -0.509341388, -1.0194830529999999, -1.4620201899999998, -1.6593975699999999, -1.9162418300000001, -2.1751451599999996, -2.08802865, 1.1936889599999998, 1.6282454999999996, 1.98038443, 2.2787639, 2.5212123799999997, 3.0547293799999995, 3.45237717, 3.8751015399999997, 4.39497966, 5.02587727, 5.428108, 5.5592031099999994, 5.17862293, 4.52390444, 3.7712036099999997, 3.24100829, 3.0616124499999993, 2.4493413, 2.1321680699999996, 1.8059909099999993, 1.3896491999999998, 1.10740258, 1.01722209, 1.2527587799999995, 1.383929289999999, 1.1199507999999998, 0.8955645939999999, 0.8664684919999998, 0.718644305, 1.2697459399999997, 1.5384010899999998, 1.3599760399999998, 1.33049311, 1.7117116899999998, 1.8688401799999999, 2.2227596599999995, 2.0889950099999997, 1.9282291099999995, 2.10277974, 2.37547602, 2.3456858299999994, 2.13734678, 2.1250930599999998, 2.3301146, 2.1542206200000003, 1.93960679, 1.36700592, 1.38875182, 1.9104304599999997])

##
########################################################

graph_sig_2sigma = TGraph(len(masses_band),masses_band,sig_2sigma)
graph_sig_2sigma.SetFillColor(kYellow)
graph_sig_2sigma.GetXaxis().SetTitle("qq resonance mass [GeV]")
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
c.SaveAs('significance_blindExercise_qq_Run2_13TeV.pdf')
