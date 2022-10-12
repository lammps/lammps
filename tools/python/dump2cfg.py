#!/usr/bin/env python

# Script:  dump2cfg.py
# Purpose: convert a LAMMPS dump file to CFG format
# Syntax:  dump2cfg.py dumpfile Nid Ntype Nx Ny Nz cfgfile
#          dumpfile = LAMMPS dump file in native LAMMPS format
#          Nid,Ntype,Nx,Ny,Nz = columns #s for ID,type,x,y,z
#                               (usually 1,2,3,4,5)
#          cfgfile = new CFG file
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.insert(1,path)
from dump import dump
from cfg import cfg

if len(sys.argv) != 8:
  sys.exit("Syntax: dump2cfg.py dumpfile Nid Ntype Nx Ny Nz cfgfile")

dumpfile = sys.argv[1]
nid = int(sys.argv[2])
ntype = int(sys.argv[3])
nx = int(sys.argv[4])
ny = int(sys.argv[5])
nz = int(sys.argv[6])
cfgfile = sys.argv[7]

d = dump(dumpfile)
d.map(nid,"id",ntype,"type",nx,"x",ny,"y",nz,"z")
c = cfg(d)
c.one(cfgfile)
