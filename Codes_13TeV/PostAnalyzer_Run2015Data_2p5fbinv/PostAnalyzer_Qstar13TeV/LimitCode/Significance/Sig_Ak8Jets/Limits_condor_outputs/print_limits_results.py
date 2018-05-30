#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array


BR = 1.0

masses = array('d')
Sig = array('d')

mass_start = 1000
mass_step = 100
steps = 20

for i in range(0,steps+1):

  ##mass = mass_start + float(i)*mass_step
  mass = mass_start + i*mass_step 
  masses.append(mass)

  f = open("stats_" + str(mass) + "_f1p0.log",'r')
  for output in iter(f):

      outputlines = output.split("\n")

      for line in outputlines:
        if re.search("Significance", line):
         Sig.append(float(line.split()[2]))

  f.close()


print "mass = "
print masses
print "Significance = "
print Sig