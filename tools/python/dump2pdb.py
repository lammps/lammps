#!/usr/bin/env python

# Script:  dump2pdb.py
# Purpose: convert a LAMMPS dump file to PDB format
# Syntax:  dump2pdb.py dumpfile Nid Ntype Nx Ny Nz pdbfile template
#          dumpfile = LAMMPS dump file in native LAMMPS format
#          Nid,Ntype,Nx,Ny,Nz = columns #s for ID,type,x,y,z
#                               (usually 1,2,3,4,5)
#          pdbfile = new PDB file
#          template = PDB file to use as template for creating new PDB file
#            this arg is optional, if not used a generic PDB file is created
# Author:  Steve Plimpton (Sandia), sjplimp at gmail.com

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.insert(1,path)
from dump import dump
from pdbfile import pdbfile

if len(sys.argv) != 8 and len(sys.argv) != 9:
  sys.exit("Syntax: dump2pdb.py dumpfile Nid Ntype Nx Ny Nz pdbfile template")

dumpfile = sys.argv[1]
nid = int(sys.argv[2])
ntype = int(sys.argv[3])
nx = int(sys.argv[4])
ny = int(sys.argv[5])
nz = int(sys.argv[6])
pfile = sys.argv[7]
if len(sys.argv) == 9: template = sys.argv[8]
else: template = ""

d = dump(dumpfile)
d.map(nid,"id",ntype,"type",nx,"x",ny,"y",nz,"z")
if template: p = pdbfile(template,d)
else: p = pdbfile(d)
p.one(pfile)
