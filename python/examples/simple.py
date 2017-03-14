#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# simple.py
# Purpose: mimic operation of examples/COUPLE/simple/simple.cpp via Python

# Serial syntax: simple.py in.lammps
#                in.lammps = LAMMPS input script

# Parallel syntax: mpirun -np 4 simple.py in.lammps
#                  in.lammps = LAMMPS input script
# also need to uncomment either Pypar or mpi4py sections below

from __future__ import print_function
import sys
import numpy as np
import ctypes

# parse command line

argv = sys.argv
if len(argv) != 2:
  print("Syntax: simple.py in.lammps")
  sys.exit()

infile = sys.argv[1]

me = 0

# uncomment this if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

# uncomment this if running in parallel via mpi4py
#from mpi4py import MPI
#me = MPI.COMM_WORLD.Get_rank()
#nprocs = MPI.COMM_WORLD.Get_size()

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
v = lmp.gather_atoms("v",1,3)
epsilon = 0.1
x[0] += epsilon
lmp.scatter_atoms("x",1,3,x)
lmp.command("run 1");

# extract force on single atom two different ways

f = lmp.extract_atom("f",3)
print("Force on 1 atom via extract_atom: ",f[0][0])

fx = lmp.extract_variable("fx","all",1)
print("Force on 1 atom via extract_variable:",fx[0])

# use commands_string() and commands_list() to invoke more commands

strtwo = "run 10\nrun 20"
lmp.commands_string(strtwo)

cmds = ["run 10","run 20"]
lmp.commands_list(cmds)

# delete all atoms
# create_atoms() to create new ones with old coords, vels
# initial thermo should be same as step 20

natoms = lmp.get_natoms()
type = natoms*[1]

lmp.command("delete_atoms group all");
lmp.create_atoms(natoms,None,type,x,v);
lmp.command("run 10");

# uncomment if running in parallel via Pypar
#print("Proc %d out of %d procs has" % (me,nprocs), lmp)
#pypar.finalize()

# uncomment if running in parallel via mpi4py
#print "Proc %d out of %d procs has" % (me,nprocs), lmp
#MPI.Finalize()
