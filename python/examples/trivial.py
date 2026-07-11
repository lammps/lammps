#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# trivial.py
# Purpose: run a LAMMPS input script via Python
# Syntax:  trivial.py in.lammps
#          in.lammps = LAMMPS input script

from __future__ import print_function
import sys

# parse command line

argv = sys.argv
if len(argv) != 2:
  print("Syntax: trivial.py in.lammps")
  sys.exit()

infile = sys.argv[1]

from lammps import lammps
lmp = lammps()

# run infile all at once

lmp.file(infile)

# or run infile one line at a time

#lines = open(infile,'r').readlines()
#for line in lines: lmp.command(line)
