#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# viz_pymol.py
# Purpose: viz running LAMMPS simulation via PyMol
# Syntax:  viz_pymol.py in.lammps Nfreq Nsteps
#          in.lammps = LAMMPS input script
#          Nfreq = dump and viz shapshot every this many steps
#          Nsteps = run for this many steps

from __future__ import print_function
import sys
sys.path.append("./pizza")

# parse command line

argv = sys.argv
if len(argv) != 4:
  print("Syntax: viz_pymol.py in.lammps Nfreq Nsteps")
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
# dump a file in native LAMMPS dump format for Pizza.py dump tool

lmp.file(infile)
lmp.command("thermo %d" % nfreq)
lmp.command("dump python all atom %d tmp.dump" % nfreq)

# initial 0-step run to generate dump file and image

lmp.command("run 0 pre yes post no")
ntimestep = 0

# wrapper on PyMol
# just proc 0 handles reading of dump file and viz

if me == 0:
  import pymol
  pymol.finish_launching()

  from dump import dump
  from pdbfile import pdbfile
  from pymol import cmd as pm

  d = dump("tmp.dump",0)
  p = pdbfile(d)
  d.next()
  d.unscale()
  p.single(ntimestep)
  pm.load("tmp.pdb")
  pm.show("spheres","tmp")
  
# run nfreq steps at a time w/out pre/post, read dump snapshot, display it

while ntimestep < nsteps:
  lmp.command("run %d pre no post no" % nfreq)
  ntimestep += nfreq
  if me == 0:
    d.next()
    d.unscale()
    p.single(ntimestep)
    pm.load("tmp.pdb")
    pm.forward()
    
lmp.command("run 0 pre no post yes")

# uncomment if running in parallel via Pypar
#print("Proc %d out of %d procs has" % (me,nprocs), lmp)
#pypar.finalize()
