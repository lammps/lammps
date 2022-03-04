#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# split.py
# Purpose: similar to simple.py, but first the world communicator
#   is split in two halves and LAMMPS is run only on one partition
# Syntax:  split.py in.lammps
#          in.lammps = LAMMPS input script

from __future__ import print_function
import sys

# parse command line

argv = sys.argv
if len(argv) != 2:
  print("Syntax: simple.py in.lammps")
  sys.exit()

infile = sys.argv[1]

me = 0

# this example *only* works with mpi4py version 2.0.0 or later

from mpi4py import MPI
comm = MPI.COMM_WORLD
me = comm.Get_rank()
nprocs = comm.Get_size()

# create two subcommunicators

if me < nprocs // 2:  color = 0
else: color = 1
  
split = comm.Split(color,key=0)

if color == 0:
  from lammps import lammps
  lmp = lammps(comm=split)

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
  epsilon = 0.1
  x[0] += epsilon
  lmp.scatter_atoms("x",1,3,x)
  lmp.command("run 1");

  f = lmp.extract_atom("f")
  print("Force on 1 atom via extract_atom: ",f[0][0])

  fx = lmp.extract_variable("fx","all",1)
  print("Force on 1 atom via extract_variable:",fx[0])
  print("Proc %d out of %d procs has" % (me,nprocs), lmp)
  print("Calculation on partition 0 complete")

else:
  # could run a 2nd calculation on second partition
  #   with different LAMMPS instance or another code
  # in this case, just sleep on second partition
  
  import time
  time.sleep(2)
  print("Calculation on partition 1 complete")

# shutdown mpi4py
  
comm.Barrier()
MPI.Finalize()
