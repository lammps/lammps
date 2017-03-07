#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# viz_vmd.py
# Purpose: viz running LAMMPS simulation via VMD
# Syntax:  viz_vmd.py in.lammps Nfreq Nsteps
#          in.lammps = LAMMPS input script
#          Nfreq = dump and viz shapshot every this many steps
#          Nsteps = run for this many steps

from __future__ import print_function
import sys
sys.path.append("./pizza")

# parse command line

argv = sys.argv
if len(argv) != 4:
  print("Syntax: viz_vmd.py in.lammps Nfreq Nsteps")
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

lmp.file(infile)
lmp.command("thermo %d" % nfreq)
lmp.command("dump python all atom %d tmp.dump" % nfreq)

# initial 0-step run to generate dump file and image

lmp.command("run 0 pre yes post no")
ntimestep = 0

# wrapper on VMD window via Pizza.py vmd tool
# just proc 0 handles reading of dump file and viz

if me == 0:
  from vmd import vmd
  v = vmd()
  v('menu main off')
  v.rep('VDW')

  from dump import dump
  from pdbfile import pdbfile

  d = dump('tmp.dump',0)
  p = pdbfile(d)
  d.next()
  d.unscale()
  p.single(ntimestep)
  v.new('tmp.pdb','pdb')

# run nfreq steps at a time w/out pre/post, read dump snapshot, display it

while ntimestep < nsteps:
  lmp.command("run %d pre no post no" % nfreq)
  ntimestep += nfreq
  if me == 0:
    d.next()
    d.unscale()
    p.single(ntimestep)
    # add frame to current data set
    v.append('tmp.pdb','pdb')
    # delete all frame and add new.
    #v.update('tmp.dump')

lmp.command("run 0 pre no post yes")

if me == 0:
  v.flush()
  # uncomment the following, if you want to work with the viz some more.
  #v('menu main on')
  #print("type quit to terminate.")
  #v.enter()
  #v.stop()

# uncomment if running in parallel via Pypar
#print("Proc %d out of %d procs has" % (me,nprocs), lmp)
#pypar.finalize()
