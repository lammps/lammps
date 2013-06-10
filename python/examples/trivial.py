#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# trivial.py
# Purpose: run a LAMMPS input script via Python
# Syntax:  trivial.py in.lammps
#          in.lammps = LAMMPS input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 2:
  print "Syntax: trivial.py in.lammps"
  sys.exit()

infile = sys.argv[1]

me = 0
# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from lammps import lammps
lmp = lammps()

# run infile all at once

lmp.file(infile)

# run infile one line at a time

#lines = open(infile,'r').readlines()
#for line in lines: lmp.command(line)

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
