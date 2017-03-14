#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# mc.py
# Purpose: mimic operation of example/MC/in.mc via Python
# Syntax:  mc.py in.mc
#          in.mc = LAMMPS input script

from __future__ import print_function
import sys,random,math

# set these parameters
# make sure neigh skin (in in.mc) > 2*deltamove

nloop = 3000
deltaperturb = 0.2
deltamove = 0.1
kT = 0.05
random.seed(27848)

# parse command line

argv = sys.argv
if len(argv) != 2:
  print("Syntax: mc.py in.mc")
  sys.exit()

infile = sys.argv[1]

from lammps import lammps
lmp = lammps()

# run infile one line at a time
# just sets up MC problem

lines = open(infile,'r').readlines()
for line in lines: lmp.command(line)
lmp.command("variable e equal pe")

# run 0 to get energy of perfect lattice
# emin = minimum energy

lmp.command("run 0")

natoms = lmp.extract_global("natoms",0)
emin = lmp.extract_compute("thermo_pe",0,0) / natoms
lmp.command("variable emin equal $e")

# disorder the system
# estart = initial energy

x = lmp.extract_atom("x",3)

for i in range(natoms):
  x[i][0] += deltaperturb * (2*random.random()-1)
  x[i][1] += deltaperturb * (2*random.random()-1)

lmp.command("variable elast equal $e")
lmp.command("thermo_style custom step v_emin v_elast pe")
lmp.command("run 0")
x = lmp.extract_atom("x",3)
lmp.command("variable elast equal $e")
  
estart = lmp.extract_compute("thermo_pe",0,0) / natoms

# loop over Monte Carlo moves
# extract x after every run, in case reneighboring changed ptr in LAMMPS

elast = estart
naccept = 0

for i in range(nloop):
  iatom = random.randrange(0,natoms)
  x0 = x[iatom][0]
  y0 = x[iatom][1]

  x[iatom][0] += deltamove * (2*random.random()-1)
  x[iatom][1] += deltamove * (2*random.random()-1)

  lmp.command("run 1 pre no post no")
  x = lmp.extract_atom("x",3)
  e = lmp.extract_compute("thermo_pe",0,0) / natoms

  if e <= elast:
    elast = e
    lmp.command("variable elast equal $e")
    naccept += 1
  elif random.random() <= math.exp(natoms*(elast-e)/kT):
    elast = e
    lmp.command("variable elast equal $e")
    naccept += 1
  else:
    x[iatom][0] = x0
    x[iatom][1] = y0
    
# final energy and stats

lmp.command("variable nbuild equal nbuild")
nbuild = lmp.extract_variable("nbuild",None,0)

lmp.command("run 0")
estop = lmp.extract_compute("thermo_pe",0,0) / natoms

print("MC stats:")
print("  starting energy =",estart)
print("  final energy =",estop)
print("  minimum energy of perfect lattice =",emin)
print("  accepted MC moves =",naccept)
print("  neighbor list rebuilds =",nbuild)
