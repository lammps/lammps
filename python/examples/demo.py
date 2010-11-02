#!/usr/local/bin/python -i
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
from lammps import LMPINT as INT
from lammps import LMPDOUBLE as DOUBLE
from lammps import LMPIPTR as IPTR
from lammps import LMPDPTR as DPTR
from lammps import LMPDPTRPTR as DPTRPTR

lmp = lammps()

# test out various library functions after running in.demo

lmp.file("in.demo")

if me == 0: print "\nPython output:"

natoms = lmp.extract_global("natoms",DOUBLE)
mass = lmp.extract_atom("mass",DPTR)
x = lmp.extract_atom("x",DPTRPTR)
print "Natoms, mass, x[0][0] coord =",natoms,mass[1],x[0][0]

temp = lmp.extract_compute("thermo_temp",0,0)
print "Temperature from compute =",temp

eng = lmp.extract_variable("eng",None,0)
print "Energy from equal-style variable =",eng

vy = lmp.extract_variable("vy","all",1)
print "Velocity component from atom-style variable =",vy[1]

natoms = lmp.get_natoms()
print "Natoms from get_natoms =",natoms

xc = lmp.get_coords()
print "Global coords from get_coords =",xc[0],xc[1],xc[31]

xc[0] = xc[0] + 1.0
lmp.put_coords(xc)

print "Changed x[0][0] via put_coords =",x[0][0]

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
