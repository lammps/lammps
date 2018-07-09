#!/usr/bin/python2.7

import os,sys
from fractions import Fraction
from numpy import *

### PARAMETERS ###
sigma = 1.00
epsilon = 1.00

b_hh = 0.95 * sigma
b_ht = 0.95 * sigma
b_tt = 1.00 * sigma

r_init = 0.000001
r_max = sigma * 3.
r_space = 0.01
##################

### INPUTS ###
if len(sys.argv) == 2:
  w_cut = float(sys.argv[1])
else:
  w_cut = 1.6
# 1.6 seems to be 'good' for vesicles, bilayers 1.4
##############

def WCA_energy(b, r):
# Calculate WCA energy
  E_pot = 0
  val1 = math.pow((b / r), 12)
  val2 = -math.pow((b / r), 6)
  E_pot = 4 * epsilon * (val1 + val2 + 0.25)
  return E_pot

def WCA_forces(b, r):
# Calculate WCA forces
  Force = 0
  val1 = 24 * math.pow(b, 6) / math.pow(r, 7)
  val2 = -48 * math.pow(b, 12) / math.pow(r, 13)
  Force = -(val1 + val2)
  return Force

def Tail_energy(b, r, r_cut):
# Calculate extra Tail energy
  E_pot = 0
  if (r < r_cut):
    E_pot = -1 * epsilon
  else:
    val1 = math.cos((math.pi * (r - r_cut)) / (2 * w_cut))
    E_pot = -1 * epsilon * math.pow(val1, 2)
  return E_pot

def Tail_forces(b, r, r_cut):
  Force = 0
  if (r < r_cut):
    Force = 0;
  else:
    val1 = math.sin((math.pi * (r - r_cut)) / w_cut)
    Force = -math.pi * val1 / (2 * w_cut)
  return Force


##############
ofile = open('tabulated_potential.dat', 'w')
tot_potential_hh = zeros((int(r_max / r_space) + 1, 4))
tot_potential_ht = zeros((int(r_max / r_space) + 1, 4))
tot_potential_tt = zeros((int(r_max / r_space) + 1, 4))

# Setup up formatting & distances in all arrays
for i in range(int(r_max / r_space)+1):
  tot_potential_hh[:,0][i] = i+1
  tot_potential_ht[:,0][i] = i+1
  tot_potential_tt[:,0][i] = i+1
for i in range(1, int(r_max / r_space)+1):
  tot_potential_hh[:,1][i] = tot_potential_hh[:,1][i-1] + r_space
  tot_potential_ht[:,1][i] = tot_potential_ht[:,1][i-1] + r_space
  tot_potential_tt[:,1][i] = tot_potential_tt[:,1][i-1] + r_space
tot_potential_hh[:,1][0] = r_init
tot_potential_ht[:,1][0] = r_init
tot_potential_tt[:,1][0] = r_init



ofile.write("# Tabulated potential for Cooke 3-bead lipid model, Wc = %f\n\n" % w_cut)
num = len(tot_potential_hh[:,0])

### Calcaulte first potential, H-H
ofile.write("HEAD_HEAD\n")
r_cut = 2**Fraction('1/6') * b_hh
rmax = int(r_cut / r_space)
ofile.write("N %d R %f %f\n\n" % (num, r_init, r_max))
ofile.write("1 %f %f %f\n" % (tot_potential_hh[:,1][0], tot_potential_hh[:,2][0], tot_potential_hh[:,3][0]))

for i in range(1, rmax+1):
  tot_potential_hh[:,2][i] = WCA_energy(b_hh, tot_potential_hh[:,1][i])
  tot_potential_hh[:,3][i] = WCA_forces(b_hh, tot_potential_hh[:,1][i])

for i in range(1, int(r_max / r_space)+1):
  ofile.write("%d %f %f %f\n" % (i+1, tot_potential_hh[:,1][i], tot_potential_hh[:,2][i], tot_potential_hh[:,3][i]))
ofile.write("\n")



### Calcaulte second potential, H-T
ofile.write("HEAD_TAIL\n")
r_cut = 2**Fraction('1/6') * b_ht
rmax = int(r_cut / r_space)
ofile.write("N %d R %f %f\n\n" % (num, r_init, r_max))
ofile.write("1 %f %f %f\n" % (tot_potential_ht[:,1][0], tot_potential_ht[:,2][0], tot_potential_ht[:,3][0]))

for i in range(1, rmax+1):
  tot_potential_ht[:,2][i] = WCA_energy(b_ht, tot_potential_ht[:,1][i])
  tot_potential_ht[:,3][i] = WCA_forces(b_ht, tot_potential_ht[:,1][i])

for i in range(1, int(r_max / r_space)+1):
  ofile.write("%d %f %f %f\n" % (i+1, tot_potential_ht[:,1][i], tot_potential_ht[:,2][i], tot_potential_ht[:,3][i]))
ofile.write("\n")



### Calcaulte third potential, T-T
# Also include extra tail-tail attraction term
ofile.write("TAIL_TAIL\n")
r_cut = 2**Fraction('1/6') * b_tt
rmax = int(r_cut / r_space)
ofile.write("N %d R %f %f\n\n" % (num, r_init, r_max))
ofile.write("1 %f %f %f\n" % (tot_potential_tt[:,1][0], tot_potential_tt[:,2][0], tot_potential_tt[:,3][0]))

for i in range(1, rmax+1):
  tot_potential_tt[:,2][i] = WCA_energy(b_tt, tot_potential_tt[:,1][i])
  tot_potential_tt[:,3][i] = WCA_forces(b_tt, tot_potential_tt[:,1][i])

max2 = int( (r_cut + w_cut) / r_space)
for i in range(1, max2+1):
  tot_potential_tt[:,2][i] = tot_potential_tt[:,2][i] + Tail_energy(b_tt, tot_potential_tt[:,1][i], r_cut)
  tot_potential_tt[:,3][i] = tot_potential_tt[:,3][i] + Tail_forces(b_tt, tot_potential_tt[:,1][i], r_cut)

for i in range(1, int(r_max / r_space)+1):
  ofile.write("%d %f %f %f\n" % (i+1, tot_potential_tt[:,1][i], tot_potential_tt[:,2][i], tot_potential_tt[:,3][i]))
ofile.write("\n")


sys.exit()
