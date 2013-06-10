#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# vizplotgui_vmd.py
# Purpose: viz running LAMMPS simulation via VMD with plot and GUI
# Syntax:  vizplotgui_vmd.py in.lammps Nfreq compute-ID
#          in.lammps = LAMMPS input script
#          Nfreq = plot data point and viz shapshot every this many steps
#          compute-ID = ID of compute that calculates temperature
#                       (or any other scalar quantity)

# IMPORTANT: this script cannot yet be run in parallel via Pypar,
#            because I can't seem to do a MPI-style broadcast in Pypar

import sys,time
sys.path.append("./pizza")

# methods called by GUI

def run():
  global runflag
  runflag = 1
def stop():
  global runflag
  runflag = 0
def settemp(value):
  global temptarget
  temptarget = slider.get()
def quit():
  global breakflag
  breakflag = 1

# method called by timestep loop every Nfreq steps
# read dump snapshot and viz it, update plot with compute value

def update(ntimestep):
  d.next()
  d.unscale()
  p.single(ntimestep)
  v.append('tmp.pdb','pdb')
  value = lmp.extract_compute(compute,0,0)
  xaxis.append(ntimestep)
  yaxis.append(value)
  gn.plot(xaxis,yaxis)

# parse command line

argv = sys.argv
if len(argv) != 4:
  print "Syntax: vizplotgui_vmd.py in.lammps Nfreq compute-ID"
  sys.exit()

infile = sys.argv[1]
nfreq = int(sys.argv[2])
compute = sys.argv[3]

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

# initial 0-step run to generate initial 1-point plot, dump file, and image

lmp.command("run 0 pre yes post no")
value = lmp.extract_compute(compute,0,0)
ntimestep = 0
xaxis = [ntimestep]
yaxis = [value]

breakflag = 0
runflag = 0
temptarget = 1.0

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

# display GUI with run/stop buttons and slider for temperature

if me == 0:
  from Tkinter import *
  tkroot = Tk()
  tkroot.withdraw()
  root = Toplevel(tkroot)
  root.title("LAMMPS GUI")

  frame = Frame(root)
  Button(frame,text="Run",command=run).pack(side=LEFT)
  Button(frame,text="Stop",command=stop).pack(side=LEFT)
  slider = Scale(frame,from_=0.0,to=5.0,resolution=0.1,
                 orient=HORIZONTAL,label="Temperature")
  slider.bind('<ButtonRelease-1>',settemp)
  slider.set(temptarget)
  slider.pack(side=LEFT)
  Button(frame,text="Quit",command=quit).pack(side=RIGHT)
  frame.pack()
  tkroot.update()

# wrapper on GnuPlot via Pizza.py gnu tool

if me == 0:
  from gnu import gnu
  gn = gnu()
  gn.plot(xaxis,yaxis)
  gn.title(compute,"Timestep","Temperature")

# endless loop, checking status of GUI settings every Nfreq steps
# run with pre yes/no and post yes/no depending on go/stop status
# re-invoke fix langevin with new seed when temperature slider changes
# after re-invoke of fix langevin, run with pre yes

running = 0
temp = temptarget
seed = 12345

lmp.command("fix 2 all langevin %g %g 0.1 %d" % (temp,temp,seed))

while 1:
  if me == 0: tkroot.update()
  if temp != temptarget:
    temp = temptarget
    seed += me+1
    lmp.command("fix 2 all langevin %g %g 0.1 12345" % (temp,temp))
    running = 0
  if runflag and running:
    lmp.command("run %d pre no post no" % nfreq)
    ntimestep += nfreq
    if me == 0: update(ntimestep)
  elif runflag and not running:
    lmp.command("run %d pre yes post no" % nfreq)
    ntimestep += nfreq
    if me == 0: update(ntimestep)
  elif not runflag and running:
    lmp.command("run %d pre no post yes" % nfreq)
    ntimestep += nfreq
    if me == 0: update(ntimestep)
  if breakflag: break
  if runflag: running = 1
  else: running = 0
  time.sleep(0.01)

lmp.command("run 0 pre no post yes")

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#pypar.finalize()
