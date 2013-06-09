#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# viz_atomeye.py
# Purpose: viz running LAMMPS simulation via AtomEye
# Syntax:  viz_atomeye.py in.lammps Nfreq Nsteps
#          in.lammps = LAMMPS input script
#          Nfreq = dump and viz shapshot every this many steps
#          Nsteps = run for this many steps

import sys,os

# set this to point to AtomEye version 3 executable
# first line if want AtomEye output to screen, 2nd line to file
#ATOMEYE3 = "/home/sjplimp/tools/atomeye3/A3.i686-20060530"
ATOMEYE3 = "/home/sjplimp/tools/atomeye3/A3.i686-20060530 > atomeye.out"

# parse command line

argv = sys.argv
if len(argv) != 4:
  print "Syntax: viz_atomeye.py in.lammps Nfreq Nsteps"
  sys.exit()

infile = sys.argv[1]
nfreq = int(sys.argv[2])
nsteps = int(sys.argv[3])

me = 0
# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from lammps import lammps
lmp = lammps()

# run infile all at once
# assumed to have no run command in it
# dump a file in extended CFG format for AtomEye

lmp.file(infile)
lmp.command("thermo %d" % nfreq)
lmp.command("dump python all cfg %d tmp.cfg.* id type xs ys zs" % nfreq)

# initial 0-step run to generate dump file and image

lmp.command("run 0 pre yes post no")
ntimestep = 0

# wrapper on GL window via Pizza.py gl tool
# just proc 0 handles reading of dump file and viz

if me == 0:
  a = os.popen(ATOMEYE3,'w')
  a.write("load_config tmp.cfg.0\n")
  a.flush()

# run nfreq steps at a time w/out pre/post, read dump snapshot, display it

while ntimestep < nsteps:
  lmp.command("run %d pre no post no" % nfreq)
  ntimestep += nfreq
  if me == 0:
    a.write("load_config tmp.cfg.%d\n" % ntimestep)
    a.flush()

lmp.command("run 0 pre no post yes")

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
