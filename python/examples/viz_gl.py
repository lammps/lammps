#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# viz_gl.py
# Purpose: viz running LAMMPS simulation via GL tool in Pizza.py
# Syntax:  viz_gl.py in.lammps Nfreq Nsteps
#          in.lammps = LAMMPS input script
#          Nfreq = dump and viz shapshot every this many steps
#          Nsteps = run for this many steps

import sys
sys.path.append("./pizza")

# parse command line

argv = sys.argv
if len(argv) != 4:
  print "Syntax: viz_gl.py in.lammps Nfreq Nsteps"
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

# wrapper on GL window via Pizza.py gl tool
# just proc 0 handles reading of dump file and viz

if me == 0:
  import Tkinter
  tkroot = Tkinter.Tk()
  tkroot.withdraw()

  from dump import dump
  from gl import gl

  d = dump("tmp.dump",0)
  g = gl(d)
  d.next()
  d.unscale()
  g.zoom(1)
  g.shift(0,0)
  g.rotate(0,270)
  g.q(10)
  g.box(1)
  g.show(ntimestep)

# run nfreq steps at a time w/out pre/post, read dump snapshot, display it

while ntimestep < nsteps:
  lmp.command("run %d pre no post no" % nfreq)
  ntimestep += nfreq
  if me == 0:
    d.next()
    d.unscale()
    g.show(ntimestep)

lmp.command("run 0 pre no post yes")

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
