#!/usr/bin/env python                                                                                                                         
import sys, os, subprocess, string, re
from ROOT import *
from array import array

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.045, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.045, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
t_m = 0.06  ##top margin                   
b_m = 0.11   ##botton margin                                        
l_m = 0.115  ##left margin                                                      
r_m = 0.04  ##right margin                                                         
gStyle.SetPadTopMargin(t_m)
gStyle.SetPadBottomMargin(b_m)
gStyle.SetPadLeftMargin(l_m)
gStyle.SetPadRightMargin(r_m)
gROOT.ForceStyle()


xs = array('d',  [10, 20, 30, 40, 50])
print xs
eff = array('d', [.1, .2, .3, .4, .5])
print eff
xs = array('d', [float(b) / float(m) for b,m in zip(xs, eff)])
print xs
