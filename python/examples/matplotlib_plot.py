#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# matplotlib_plot.py
# Purpose: plot Temp of running LAMMPS simulation via matplotlib
# Syntax:  plot.py in.lammps Nfreq Nsteps compute-ID
#          in.lammps = LAMMPS input script
#          Nfreq = plot data point every this many steps
#          Nsteps = run for this many steps
#          compute-ID = ID of compute that calculates temperature
#                       (or any other scalar quantity)

from __future__ import print_function
import sys
sys.path.append("./pizza")
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

# parse command line

argv = sys.argv
if len(argv) != 5:
  print("Syntax: plot.py in.lammps Nfreq Nsteps compute-ID")
  sys.exit()

infile = sys.argv[1]
nfreq = int(sys.argv[2])
nsteps = int(sys.argv[3])
compute = sys.argv[4]

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

# initial 0-step run to generate initial 1-point plot

lmp.command("run 0 pre yes post no")
value = lmp.extract_compute(compute,0,0)
ntimestep = 0
xaxis = [ntimestep]
yaxis = [value]

# create matplotlib plot
# just proc 0 handles plotting

if me == 0:
  fig = plt.figure()
  line, = plt.plot(xaxis, yaxis)
  plt.xlim([0, nsteps])
  plt.title(compute)
  plt.xlabel("Timestep")
  plt.ylabel("Temperature")
  plt.show(block=False)

# run nfreq steps at a time w/out pre/post, query compute, refresh plot
import time

while ntimestep < nsteps:
  lmp.command("run %d pre no post no" % nfreq)
  ntimestep += nfreq
  value = lmp.extract_compute(compute,0,0)
  xaxis.append(ntimestep)
  yaxis.append(value)
  if me == 0:
    line.set_xdata(xaxis)
    line.set_ydata(yaxis)
    ax = plt.gca()
    ax.relim()
    ax.autoscale_view(True, True, True)
    fig.canvas.draw()

lmp.command("run 0 pre no post yes")

# uncomment if running in parallel via Pypar
#print("Proc %d out of %d procs has" % (me,nprocs), lmp)
#pypar.finalize()

if sys.version_info[0] == 3:
    input("Press Enter to exit...")
else:
    raw_input("Press Enter to exit...")
