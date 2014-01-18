#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# simple.py
# Purpose: mimic operation of couple/simple/simple.cpp via Python
# Syntax:  simple.py in.lammps
#          in.lammps = LAMMPS input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 2:
  print "Syntax: simple.py in.lammps"
  sys.exit()

infile = sys.argv[1]

me = 0

# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from lammps import lammps
lmp = lammps()

# run infile one line at a time

lines = open(infile,'r').readlines()
for line in lines: lmp.command(line)

# run 10 more steps
# get coords from LAMMPS
# change coords of 1st atom
# put coords back into LAMMPS
# run a single step with changed coords

lmp.command("run 10")
x = lmp.gather_atoms("x",1,3)
epsilon = 0.1
x[0] += epsilon
lmp.scatter_atoms("x",1,3,x)
lmp.command("run 1");

f = lmp.extract_atom("f",3)
print "Force on 1 atom via extract_atom: ",f[0][0]

fx = lmp.extract_variable("fx","all",1)
print "Force on 1 atom via extract_variable:",fx[0]

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
