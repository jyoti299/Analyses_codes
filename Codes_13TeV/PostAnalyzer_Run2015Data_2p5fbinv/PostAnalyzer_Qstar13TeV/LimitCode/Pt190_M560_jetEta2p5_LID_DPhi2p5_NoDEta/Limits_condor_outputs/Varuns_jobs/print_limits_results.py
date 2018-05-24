#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array


BR = 1.0

masses = array('d')
xs_obs_limits = array('d')
xs_exp_limits = array('d')
masses_exp = array('d')
xs_exp_limits_1sigma = array('d')
xs_exp_limits_1sigma_up = array('d')
xs_exp_limits_2sigma = array('d')
xs_exp_limits_2sigma_up = array('d')

mass_start = 1000
mass_step = 200
steps = 22

for i in range(0,steps+1):

  ##mass = mass_start + float(i)*mass_step
  mass = mass_start + i*mass_step
  
  masses.append(mass)
  masses_exp.append(mass)

  f = open("stats_" + str(mass) + "_f1p0.log",'r')
  for output in iter(f):

      outputlines = output.split("\n")

      for line in outputlines:
        if re.search("observed bound =", line):
         xs_obs_limits.append(float(line.split()[6]))
#         print line.split()[6]
        if re.search("median:", line):
          xs_exp_limits.append(float(line.split()[1]))
        if re.search("1 sigma band:", line):
         xs_exp_limits_1sigma.append(float(line.split()[4]))
         xs_exp_limits_1sigma_up.append(float(line.split()[6]))
        if re.search("2 sigma band:", line):
         xs_exp_limits_2sigma.append(float(line.split()[4]))
         xs_exp_limits_2sigma_up.append(float(line.split()[6]))

  f.close()

for i in range(0,len(masses)):
  masses_exp.append( masses[len(masses)-i-1] )
  xs_exp_limits_1sigma.append( xs_exp_limits_1sigma_up[len(masses)-i-1] )
  xs_exp_limits_2sigma.append( xs_exp_limits_2sigma_up[len(masses)-i-1] )


print "masses = "
print masses
print ""
print "xs_obs_limits  = " 
print xs_obs_limits
print ""
print "xs_exp_limits = "
print xs_exp_limits
print ""
print "masses_exp = "
print masses_exp
print ""
print "xs_exp_limits_1sigma  = "
print xs_exp_limits_1sigma
print ""
print "xs_exp_limits_2sigma = "
print xs_exp_limits_2sigma
