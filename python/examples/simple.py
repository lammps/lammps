#!/usr/bin/env python -i
# preceding line should have path for Python on your machine

# simple.py
# Purpose: mimic operation of examples/COUPLE/simple/simple.cpp via Python

# Serial syntax: simple.py in.lammps
#                in.lammps = LAMMPS input script

# Parallel syntax: mpirun -np 4 simple.py in.lammps
#                  in.lammps = LAMMPS input script
# also need to uncomment mpi4py sections below

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

f = lmp.extract_atom("f")
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
atype = natoms*[1]

lmp.command("delete_atoms group all");
lmp.create_atoms(natoms,None,atype,x,v);
lmp.command("run 10");

############
# test of new gather/scatter and box extract/reset methods
# can try this in parallel and with/without atom_modify sort enabled

lmp.command("write_dump all custom tmp.simple id type x y z fx fy fz");

x = lmp.gather_atoms("x",1,3)
f = lmp.gather_atoms("f",1,3)

if me == 0: print("Gather XF:",x[3],x[9],f[3],f[9])

ids = lmp.gather_atoms_concat("id",0,1)
x = lmp.gather_atoms_concat("x",1,3)
f = lmp.gather_atoms_concat("f",1,3)

if me == 0: print("Gather concat XF:",ids[0],ids[1],x[0],x[3],f[0],f[3])

ids = (2*ctypes.c_int)()
ids[0] = 2
ids[1] = 4
x = lmp.gather_atoms_subset("x",1,3,2,ids)
f = lmp.gather_atoms_subset("f",1,3,2,ids)

if me == 0: print("Gather subset XF:",x[0],x[3],f[0],f[3])

x[0] = -1.0
x[1] = 0.0
x[2] = 0.0
x[3] = -2.0
x[4] = 0.0
x[5] = 0.0
ids[0] = 100
ids[1] = 200
lmp.scatter_atoms_subset("x",1,3,2,ids,x)

x = lmp.gather_atoms("x",1,3)
if me == 0: print("Gather post scatter subset:",
                  x[3],x[9],x[297],x[298],x[299],x[597],x[598],x[599])

boxlo,boxhi,xy,yz,xz,periodicity,box_change = lmp.extract_box()
if me == 0: print("Box info",boxlo,boxhi,xy,yz,xz,periodicity,box_change)

lmp.reset_box([0,0,0],[10,10,8],0,0,0)

boxlo,boxhi,xy,yz,xz,periodicity,box_change = lmp.extract_box()
if me == 0: print("Box info",boxlo,boxhi,xy,yz,xz,periodicity,box_change)

# uncomment if running in parallel via mpi4py
#print("Proc %d out of %d procs has" % (me,nprocs), lmp)
