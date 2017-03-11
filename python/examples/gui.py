#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# gui.py
# Purpose: control a continuously running LAMMPS simulation via a Tkinter GUI
# Syntax:  gui.py in.lammps Nfreq
#          in.lammps = LAMMPS input script
#          Nfreq = query GUI every this many steps

# IMPORTANT: this script cannot yet be run in parallel via Pypar,
#            because I can't seem to do a MPI-style broadcast in Pypar

from __future__ import print_function
import sys,time

# methods called by GUI

def go():
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

# parse command line

argv = sys.argv
if len(argv) != 3:
  print("Syntax: gui.py in.lammps Nfreq")
  sys.exit()

infile = sys.argv[1]
nfreq = int(sys.argv[2])

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

# display GUI with go/stop/quit buttons and slider for temperature
# just proc 0 handles GUI

breakflag = 0
runflag = 0
temptarget = 1.0

if me == 0:
  try:
    from Tkinter import *
  except:
    from tkinter import *
  tkroot = Tk()
  tkroot.withdraw()
  root = Toplevel(tkroot)
  root.title("LAMMPS GUI")

  frame = Frame(root)
  Button(frame,text="Go",command=go).pack(side=LEFT)
  Button(frame,text="Stop",command=stop).pack(side=LEFT)
  slider = Scale(frame,from_=0.0,to=5.0,resolution=0.1,
                 orient=HORIZONTAL,label="Temperature")
  slider.bind('<ButtonRelease-1>',settemp)
  slider.set(temptarget)
  slider.pack(side=LEFT)
  Button(frame,text="Quit",command=quit).pack(side=RIGHT)
  frame.pack()
  tkroot.update()

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
    lmp.command("fix 2 all langevin %g %g 0.1 %d" % (temp,temp,seed))
    running = 0
  if runflag and running:
    lmp.command("run %d pre no post no" % nfreq)
  elif runflag and not running:
    lmp.command("run %d pre yes post no" % nfreq)
  elif not runflag and running:
    lmp.command("run %d pre no post yes" % nfreq)
  if breakflag: break
  if runflag: running = 1
  else: running = 0
  time.sleep(0.01)

# uncomment if running in parallel via Pypar
#print("Proc %d out of %d procs has" % (me,nprocs), lmp)
#pypar.finalize()
