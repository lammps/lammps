#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# demo.py
# Purpose: illustrate use of many library interface commands
# Syntax:  demo.py
#          uses in.demo as LAMMPS input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 1:
  print "Syntax: demo.py"
  sys.exit()

me = 0
# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from lammps import lammps

lmp = lammps()

# test out various library functions after running in.demo

lmp.file("in.demo")

if me == 0: print "\nPython output:"

natoms = lmp.extract_global("natoms",0)
mass = lmp.extract_atom("mass",2)
x = lmp.extract_atom("x",3)
print "Natoms, mass, x[0][0] coord =",natoms,mass[1],x[0][0]

temp = lmp.extract_compute("thermo_temp",0,0)
print "Temperature from compute =",temp

eng = lmp.extract_variable("eng",None,0)
print "Energy from equal-style variable =",eng

vy = lmp.extract_variable("vy","all",1)
print "Velocity component from atom-style variable =",vy[1]

natoms = lmp.get_natoms()
print "Natoms from get_natoms =",natoms

xc = lmp.gather_atoms("x",1,3)
print "Global coords from gather_atoms =",xc[0],xc[1],xc[31]

xc[0] = xc[0] + 1.0
lmp.scatter_atoms("x",1,3,xc)

print "Changed x[0][0] via scatter_atoms =",x[0][0]

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
