#!/usr/bin/env python

# Script:  dumpsort.py
# Purpose: sort the snapshots in a LAMMPS dump file by atom ID
# Syntax:  dumpsort.py oldfile N newfile
#          oldfile = old LAMMPS dump file in native LAMMPS format
#          N = column # for atom ID (usually 1)
#          newfile = new sorted LAMMPS dump file
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.insert(1,path)
from dump import dump

if len(sys.argv) != 4:
  sys.exit("Syntax: dumpsort.py oldfile N newfile")

oldfile = sys.argv[1]
ncolumn = int(sys.argv[2])
newfile = sys.argv[3]

d = dump(oldfile)
d.map(ncolumn,"id")
d.sort()
d.write(newfile)
