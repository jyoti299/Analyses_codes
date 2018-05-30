#!/usr/bin/env python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
import numpy as np
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
xs_obs_limits = array('d')
xs_exp_limits = array('d')
masses_exp = array('d')
xs_exp_limits_1sigma = array('d')
xs_exp_limits_1sigma_up = array('d')
xs_exp_limits_2sigma = array('d')
xs_exp_limits_2sigma_up = array('d')


syst = True
#syst = False

useTeV = True
#useTeV = False

mass_min = 1500
mass_max = 7200

########################################################
## Uncomment this part if running the limit code


### for running the limit code
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(float(mass))
  #masses_exp.append(float(mass))

  #cmd = "./stats " + str(int(mass)) + " gg | tee stats_" + str(int(mass)) + "_gg.log"
  #print "Running: " + cmd
  #proc = subprocess.Popen( cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
  #output = proc.communicate()[0]
  #if proc.returncode != 0:
    #print output
    #sys.exit(1)
  ##print output

  #outputlines = output.split("\n")

  #for line in outputlines:
    #if re.search("observed bound =", line):
      #xs_obs_limits.append(float(line.split()[6]))
    #if re.search("median:", line):
      #xs_exp_limits.append(float(line.split()[1]))
    #if re.search("1 sigma band:", line):
      #xs_exp_limits_1sigma.append(float(line.split()[4]))
      #xs_exp_limits_1sigma_up.append(float(line.split()[6]))
    #if re.search("2 sigma band:", line):
      #xs_exp_limits_2sigma.append(float(line.split()[4]))
      #xs_exp_limits_2sigma_up.append(float(line.split()[6]))

##------------------------------------------------------

### for reading the limit code log files
#for mass in range(mass_min,mass_max+100,100):

  #masses.append(float(mass))
  #masses_exp.append(float(mass))

  #log_file = open("stats_" + str(int(mass)) + "_gg.log",'r')
  #outputlines = log_file.readlines()
  #log_file.close()

  #for line in outputlines:
    #if re.search("observed bound =", line):
      #xs_obs_limits.append(float(line.split()[6]))
    #if re.search("median:", line):
      #xs_exp_limits.append(float(line.split()[1]))
    #if re.search("1 sigma band:", line):
      #xs_exp_limits_1sigma.append(float(line.split()[4]))
      #xs_exp_limits_1sigma_up.append(float(line.split()[6]))
    #if re.search("2 sigma band:", line):
      #xs_exp_limits_2sigma.append(float(line.split()[4]))
      #xs_exp_limits_2sigma_up.append(float(line.split()[6]))

##------------------------------------------------------

#for i in range(0,len(masses)):
  #masses_exp.append( masses[len(masses)-i-1] )
  #xs_exp_limits_1sigma.append( xs_exp_limits_1sigma_up[len(masses)-i-1] )
  #xs_exp_limits_2sigma.append( xs_exp_limits_2sigma_up[len(masses)-i-1] )


#print "masses =", masses
#print "xs_obs_limits =", xs_obs_limits
#print "xs_exp_limits =", xs_exp_limits
#print ""
#print "masses_exp =", masses_exp
#print "xs_exp_limits_1sigma =", xs_exp_limits_1sigma
#print "xs_exp_limits_2sigma =", xs_exp_limits_2sigma

##
########################################################

########################################################
## Comment out this part if running the limit code

masses = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0])
xs_obs_limits = array('d', [0.836478, 0.774447, 1.23998, 0.959869, 0.526238, 0.310743, 0.219122, 0.14414, 0.117179, 0.131122, 0.127128, 0.0906086, 0.104278, 0.117512, 0.122472, 0.120322, 0.121221, 0.113741, 0.102283, 0.0975416, 0.0961606, 0.0925603, 0.0792465, 0.0676229, 0.0701066, 0.0671454, 0.0527954, 0.0336437, 0.0224766, 0.0150056, 0.0121905, 0.0103843, 0.00921951, 0.00866013, 0.00924763, 0.0110804, 0.0127599, 0.0139334, 0.0141779, 0.013463, 0.0126392, 0.0119217, 0.0112133, 0.010475, 0.00967835, 0.00879789, 0.00802754, 0.00750105, 0.00693862, 0.00675818, 0.00620987, 0.00612767, 0.00588462, 0.00575895, 0.00572472, 0.00575174, 0.00590138, 0.00613816])
xs_exp_limits = array('d', [0.807865, 0.749822, 0.636492, 0.501165, 0.442055, 0.373683, 0.308796, 0.274763, 0.2461, 0.21873, 0.197548, 0.159459, 0.144599, 0.12325, 0.112418, 0.101445, 0.092047, 0.0772591, 0.0728687, 0.0636055, 0.0546508, 0.0506898, 0.0467854, 0.0444067, 0.0383595, 0.0337996, 0.0299116, 0.0287613, 0.0256873, 0.0245892, 0.0222354, 0.0203452, 0.0190944, 0.0182157, 0.015392, 0.0142076, 0.0136367, 0.0126589, 0.012195, 0.0112466, 0.0102815, 0.00954721, 0.00920796, 0.00848258, 0.00804014, 0.00762756, 0.0071216, 0.00727543, 0.00694033, 0.00698688, 0.00702847, 0.0068209, 0.00682605, 0.00684432, 0.00685966, 0.00686981, 0.00704506, 0.00730533])

masses_exp = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0, 7200.0, 7100.0, 7000.0, 6900.0, 6800.0, 6700.0, 6600.0, 6500.0, 6400.0, 6300.0, 6200.0, 6100.0, 6000.0, 5900.0, 5800.0, 5700.0, 5600.0, 5500.0, 5400.0, 5300.0, 5200.0, 5100.0, 5000.0, 4900.0, 4800.0, 4700.0, 4600.0, 4500.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0])
xs_exp_limits_1sigma = array('d', [0.335161, 0.324007, 0.321838, 0.286281, 0.244971, 0.197079, 0.180205, 0.161407, 0.133198, 0.123439, 0.111531, 0.0981806, 0.0907468, 0.0797318, 0.0713594, 0.0619089, 0.0561893, 0.0500096, 0.0466732, 0.0440751, 0.0365057, 0.0329189, 0.0314334, 0.0272091, 0.0248765, 0.0215591, 0.0203825, 0.0193157, 0.0167753, 0.0161453, 0.0144415, 0.0138956, 0.0131834, 0.0113471, 0.0106243, 0.00948066, 0.00890732, 0.00886268, 0.00798395, 0.00759016, 0.0071914, 0.00690173, 0.00669298, 0.00640444, 0.00613202, 0.00569241, 0.00544685, 0.00550536, 0.00547637, 0.00539133, 0.00558063, 0.00529489, 0.00523485, 0.00533668, 0.00538193, 0.00545733, 0.00560797, 0.00592521, 0.00967265, 0.00941806, 0.00932741, 0.0092946, 0.00939976, 0.0095617, 0.00969715, 0.00958647, 0.0100849, 0.0103257, 0.0104162, 0.0104984, 0.0113076, 0.0116728, 0.0119872, 0.0138671, 0.0147656, 0.0156666, 0.0170438, 0.0177269, 0.0186168, 0.019711, 0.0222257, 0.0242314, 0.0273703, 0.0284973, 0.031896, 0.0347668, 0.0380603, 0.0399196, 0.0454472, 0.049391, 0.0538697, 0.0570296, 0.0699744, 0.0767439, 0.078204, 0.0886945, 0.103308, 0.120837, 0.127646, 0.145325, 0.161374, 0.180015, 0.204109, 0.245011, 0.282155, 0.331893, 0.368612, 0.424991, 0.493316, 0.590081, 0.666704, 0.819146, 0.949194, 1.29105, 1.60552, 2.11643])
xs_exp_limits_2sigma = array('d', [0.20407, 0.209338, 0.195197, 0.174303, 0.15069, 0.123191, 0.111347, 0.0899456, 0.0885408, 0.0818745, 0.0709853, 0.0642, 0.0604732, 0.0561892, 0.0465196, 0.0422077, 0.0399515, 0.0344268, 0.0330462, 0.0293144, 0.0275114, 0.0242617, 0.0230381, 0.0192778, 0.0174884, 0.0161142, 0.0145942, 0.0138156, 0.0130697, 0.0112503, 0.0107593, 0.0100603, 0.0087881, 0.00825536, 0.00767079, 0.00718252, 0.00678177, 0.00672023, 0.00649652, 0.0064116, 0.00612197, 0.00564294, 0.00557085, 0.00531687, 0.00493033, 0.00458629, 0.00448363, 0.00431258, 0.00431956, 0.00428019, 0.00423466, 0.00422984, 0.00426799, 0.00433075, 0.00442372, 0.00453681, 0.00470056, 0.00494792, 0.0150724, 0.0145872, 0.013956, 0.013989, 0.0140444, 0.0141739, 0.0145428, 0.0136323, 0.0144526, 0.0152436, 0.0154896, 0.0148061, 0.016124, 0.0165779, 0.0164022, 0.0198081, 0.0207638, 0.0223936, 0.0231302, 0.0239143, 0.0268797, 0.0287656, 0.0322244, 0.0365646, 0.039878, 0.0421563, 0.0459401, 0.0502642, 0.0560174, 0.0573951, 0.0650935, 0.0720316, 0.0829751, 0.083066, 0.0944203, 0.109722, 0.109963, 0.126134, 0.142736, 0.16952, 0.184624, 0.215184, 0.23732, 0.263481, 0.335052, 0.358513, 0.432019, 0.501457, 0.541806, 0.63875, 0.750352, 0.921764, 1.0043, 1.25544, 1.56756, 1.91173, 2.38783, 3.80777])

if syst:
  masses = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0])
  xs_obs_limits = array('d', [2.0615, 1.94582, 2.06492, 1.62597, 1.05113, 0.619752, 0.409027, 0.278679, 0.223427, 0.231844, 0.213956, 0.17448, 0.182167, 0.187449, 0.181862, 0.17469, 0.163485, 0.151662, 0.139721, 0.129114, 0.120302, 0.11357, 0.100043, 0.0874448, 0.081939, 0.0751703, 0.0658594, 0.0515029, 0.0351077, 0.0223095, 0.0158919, 0.0128366, 0.0113868, 0.0106171, 0.0108656, 0.0121578, 0.0134528, 0.0142962, 0.0143024, 0.0139918, 0.0132415, 0.0124992, 0.0118243, 0.0110701, 0.0101539, 0.00942484, 0.00875206, 0.00818491, 0.00757354, 0.00701268, 0.0064779, 0.00638241, 0.00611084, 0.00592149, 0.0057839, 0.00583632, 0.00589959, 0.00615538])
  xs_exp_limits = array('d', [2.35524, 1.718395, 1.428665, 1.09346, 0.8535145, 0.6807505, 0.5645625, 0.4462245, 0.3892175, 0.338744, 0.296201, 0.261048, 0.227171, 0.203175, 0.18259, 0.1602235, 0.14302, 0.125324, 0.1061265, 0.09317755, 0.0843383, 0.07174945, 0.064988, 0.05691025, 0.0509401, 0.0453084, 0.04059635, 0.0373627, 0.0317369, 0.0300091, 0.0270188, 0.02477945, 0.0223704, 0.0199402, 0.0187306, 0.01701985, 0.0151752, 0.01351145, 0.01255925, 0.01157715, 0.0106024, 0.009756425, 0.00934829, 0.00875005, 0.00851081, 0.008080615, 0.00765728, 0.007629135, 0.007262815, 0.007117835, 0.007197575, 0.00695401, 0.006931495, 0.00696994, 0.00698784, 0.007031325, 0.007249375, 0.00743711])

  masses_exp = array('d', [1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0, 6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0, 7000.0, 7100.0, 7200.0, 7200.0, 7100.0, 7000.0, 6900.0, 6800.0, 6700.0, 6600.0, 6500.0, 6400.0, 6300.0, 6200.0, 6100.0, 6000.0, 5900.0, 5800.0, 5700.0, 5600.0, 5500.0, 5400.0, 5300.0, 5200.0, 5100.0, 5000.0, 4900.0, 4800.0, 4700.0, 4600.0, 4500.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0])
  xs_exp_limits_1sigma = array('d', [1.00337116, 0.866955067, 0.752425447, 0.690813204, 0.557126493, 0.452343571, 0.369974379, 0.291965862, 0.254603308, 0.22262009, 0.198962406, 0.180808615, 0.159799945, 0.146661554, 0.128943636, 0.112611294, 0.103093084, 0.0907000155, 0.0771243702, 0.0661439466, 0.0573709982, 0.0508367378, 0.0456360786, 0.0399667206, 0.0339849389, 0.030288527, 0.027896682, 0.0254099014, 0.0226699241, 0.0213114804, 0.0192591892, 0.0175941976, 0.0157615114, 0.013838817, 0.0130338082, 0.0117812166, 0.0108143444, 0.00934186314, 0.00896829673, 0.00828872597, 0.00771801739, 0.00718462537, 0.0067390125, 0.00646815324, 0.00618935429, 0.00586042785, 0.00586463704, 0.00584641849, 0.0057361607, 0.00567584782, 0.00558652705, 0.00553129486, 0.00546272191, 0.00547211088, 0.00548410478, 0.00552754675, 0.00564953595, 0.00584748945, 0.0100845357, 0.00993503261, 0.00981003639, 0.00975574841, 0.00972731409, 0.00989015291, 0.010033865, 0.0103405915, 0.0105689159, 0.0105038322, 0.0112085277, 0.0110980886, 0.0114751799, 0.012443272, 0.0131314737, 0.014051216, 0.0149478967, 0.0156688959, 0.0170780353, 0.0178464862, 0.0191880225, 0.0220464487, 0.0244916478, 0.0271920821, 0.0301781311, 0.0332523471, 0.0363692962, 0.0402007849, 0.0439346766, 0.0482859262, 0.0538109068, 0.0589490221, 0.0640103592, 0.0695765626, 0.0817522438, 0.0935356018, 0.104712428, 0.126824538, 0.135229674, 0.149406551, 0.175783888, 0.198555658, 0.222101061, 0.252941695, 0.304163464, 0.340469463, 0.379226478, 0.425511036, 0.501697685, 0.614025374, 0.685345584, 0.82294336, 0.990558043, 1.26525887, 1.69771927, 2.19449563, 3.04527554, 4.17068182])
  xs_exp_limits_2sigma = array('d', [0.54888898, 0.52204808, 0.474718403, 0.397639556, 0.351394671, 0.302690054, 0.252403694, 0.213393025, 0.189859488, 0.155648801, 0.141244894, 0.131521788, 0.115663227, 0.10740684, 0.0945093062, 0.0815758087, 0.0683574154, 0.0598277833, 0.0521576438, 0.0459937582, 0.0446251456, 0.0395625195, 0.0340091703, 0.0293762369, 0.0231166017, 0.0201350766, 0.0206927148, 0.018663006, 0.0166369141, 0.0154194806, 0.014386231, 0.0122479487, 0.0109092083, 0.0101478307, 0.0098818491, 0.00879845911, 0.00771553858, 0.00757435346, 0.00718372213, 0.00666078628, 0.00646523702, 0.00586224672, 0.0053227103, 0.005199051, 0.00494559699, 0.00471175705, 0.00503444006, 0.00485785365, 0.00479673824, 0.00460783626, 0.00438148375, 0.00444699432, 0.00444588156, 0.00454888576, 0.00460047311, 0.00468674348, 0.00487382335, 0.005121591, 0.0156015531, 0.0149743229, 0.0144452808, 0.0139797751, 0.0140042234, 0.0145114841, 0.0147583004, 0.0147875288, 0.0143771972, 0.014809522, 0.0158853178, 0.0153294932, 0.0162746371, 0.0167309479, 0.0175705089, 0.0192033326, 0.0203578441, 0.021372262, 0.023406303, 0.0251327853, 0.0268561262, 0.0307997795, 0.0346109915, 0.039699849, 0.0429967449, 0.0472399862, 0.0520788044, 0.0574647094, 0.0619649635, 0.0666650444, 0.076035464, 0.0775039182, 0.0875381938, 0.0988259657, 0.114154504, 0.132315214, 0.136144743, 0.170988029, 0.180299121, 0.199779352, 0.221410487, 0.259954698, 0.312874745, 0.362959971, 0.410513564, 0.479408071, 0.519822835, 0.578186354, 0.688012017, 0.841772799, 0.978688174, 1.18332872, 1.29074209, 1.61914179, 2.27758681, 3.17004091, 4.26499675, 5.54217356])

##
########################################################

massesS8 = array('d', [1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,3000.0,3100.0,3200.0,3300.0,3400.0,3500.0,3600.0,3700.0,3800.0,3900.0,4000.0,4100.0,4200.0,4300.0,4400.0,4500.0,4600.0,4700.0,4800.0,4900.0,5000.0,5100.0,5200.0,5300.0,5400.0,5500.0,5600.0,5700.0,5800.0,5900.0,6000.0])
xsS8 = array('d', [5.15E+02,2.93E+02,1.73E+02,1.11E+02,6.68E+01,4.29E+01,2.86E+01,1.90E+01,1.30E+01,8.71E+00,6.07E+00,4.32E+00,2.99E+00,2.14E+00,1.53E+00,1.09E+00,8.10E-01,5.83E-01,4.38E-01,3.25E-01,2.43E-01,1.78E-01,1.37E-01,1.03E-01,7.66E-02,5.76E-02,4.46E-02,3.42E-02,2.60E-02,1.94E-02,1.50E-02,1.20E-02,9.12E-03,6.99E-03,5.47E-03,4.19E-03,3.21E-03,2.53E-03,1.90E-03,1.50E-03,1.18E-03,9.13E-04,7.07E-04,5.60E-04,4.35E-04,3.36E-04,2.59E-04,2.09E-04,1.59E-04,1.21E-04,9.38E-05])

xs_max = 2e+01
idx = 0

for i, xs in enumerate(xsS8):
  if xs < xs_max:
    idx = i
    break

if useTeV:
  masses = array('d',(np.array(masses.tolist())/1000.).tolist())
  masses_exp = array('d',(np.array(masses_exp.tolist())/1000.).tolist())
  massesS8 = array('d',(np.array(massesS8.tolist())/1000.).tolist())

graph_xsS8 = TGraph(len(massesS8[idx:-1]),massesS8[idx:-1],xsS8[idx:-1])
graph_xsS8.SetLineWidth(3)
graph_xsS8.SetLineStyle(6)
graph_xsS8.SetLineColor(6)

graph_exp_2sigma = TGraph(len(masses_exp),masses_exp,xs_exp_limits_2sigma)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetXaxis().SetTitle("gg resonance mass [GeV]")
if useTeV:
  graph_exp_2sigma.GetXaxis().SetTitle("gg resonance mass [TeV]")
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #it{B} #it{A} [pb]")
graph_exp_2sigma.GetYaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetYaxis().SetRangeUser(1e-03,1e+02)
#graph_exp_2sigma.GetXaxis().SetNdivisions(1005)

graph_exp_1sigma = TGraph(len(masses_exp),masses_exp,xs_exp_limits_1sigma)
graph_exp_1sigma.SetFillColor(kGreen+1)

graph_exp = TGraph(len(masses),masses,xs_exp_limits)
#graph_exp.SetMarkerStyle(24)
graph_exp.SetLineWidth(3)
graph_exp.SetLineStyle(4)
graph_exp.SetLineColor(4)

graph_obs = TGraph(len(masses),masses,xs_obs_limits)
graph_obs.SetMarkerStyle(20)
graph_obs.SetLineWidth(3)
#graph_obs.SetLineStyle(1)
graph_obs.SetLineColor(1)


c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
graph_xsS8.Draw("L")

legend = TLegend(.59,.56,.94,.76)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetMargin(0.20)
legend.SetHeader('95% CL upper limits')
legend.AddEntry(graph_obs,"Observed","lp")
legend.AddEntry(graph_exp,"Expected","lp")
legend.AddEntry(graph_exp_1sigma,"#pm 1 std. deviation","F")
legend.AddEntry(graph_exp_2sigma,"#pm 2 std. deviation","F")
legend.Draw()

legendTh = TLegend(.59,.80,.94,.84)
legendTh.SetBorderSize(0)
legendTh.SetFillColor(0)
legendTh.SetFillStyle(0)
legendTh.SetTextFont(42)
legendTh.SetTextSize(0.035)
legendTh.SetMargin(0.20)
legendTh.AddEntry(graph_xsS8,"Color-octet scalar","l")
legendTh.Draw()

#l1 = TLatex()
#l1.SetTextAlign(12)
#l1.SetTextFont(42)
#l1.SetNDC()
#l1.SetTextSize(0.04)
#l1.SetTextSize(0.04)
#l1.DrawLatex(0.18,0.40, "CMS Preliminary")
#l1.DrawLatex(0.18,0.32, "#intLdt = 1 fb^{-1}")
#l1.DrawLatex(0.19,0.27, "#sqrt{s} = 13 TeV")

#draw the lumi text on the canvas
#CMS_lumi.extraText = "Preliminary"
CMS_lumi.extraText = ""
CMS_lumi.lumi_sqrtS = "2.4 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.15
iPeriod = 0

CMS_lumi.CMS_lumi(c, iPeriod, iPos)

gPad.RedrawAxis()

c.SetLogy()
c.SaveAs('xs_limit_DijetLimitCode_gg' + ('_NoSyst' if not syst else '') + '_Run2_13TeV_DATA_2445_invpb.eps')