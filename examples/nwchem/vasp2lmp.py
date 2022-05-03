#!/bin/env python

# convert VASP POSCAR file(s) to LAMMPS data file(s)

# Syntax: python vasp2lmp.py vaspfile datafile
#   if both args are files, vaspfile is POSCAR input, datafile is output
#   if both args are dirs:
#     each POSCAR file in vasp dir is converted to LAMMPS data file in data dir
#     output filename = POSCAR file name prepended with "data."

# POSCAR file settings this script recognizes
#   comment on line 1 (ignored)
#   scale factor > 0.0 on line 2 for box edge vectors and Cartesian coords
#   3 box edge vectors on lines 3,4,5
#   element names on line 6
#   species counts on line 7
#   optional Selective Dynamics line 8 (ignored)
#   Cartesian or Direct line 8 (or 9 if Selective Dynamics)
#   any T/F flags at end of atom coords lines are ignored
#   any remaining lines are ignored

# ---------------------------
# error message

def error(str=None):
  if not str: print "Syntax: vasp2lmp.py POSCARfile datafile"
  else: print str
  sys.exit()
  
# ---------------------------
# convert a VASP POSCAR file to a LAMMPS data file

def convert(vaspfile,datafile):

  # read VASP POSCAR file
  
  lines = open(vaspfile,"r").readlines()

  comment = lines[0]
  scale = float(lines[1])
  if scale < 0.0: scale = 1.0
  
  avec = lines[2].split()
  avec = [scale*float(one) for one in avec]
  bvec = lines[3].split()
  bvec = [scale*float(one) for one in bvec]
  cvec = lines[4].split()
  cvec = [scale*float(one) for one in cvec]

  if avec[1] != 0.0 or avec[2] != 0.0 or bvec[2] != 0.0:
    error("Triclinic VASP box does not meet LAMMPS triclinic restrictions")
    
  elements = lines[5].split()
  counts = lines[6].split()
  if len(elements) != len(counts):
    error("POSCAR file elements and counts lines do not match")

  ntypes = len(elements)
  counts = [int(one) for one in counts]
  natoms = sum(counts)
  
  label = lines[7].strip().lower()
  firstatomline = 8

  if label[0] == 's':
    label = lines[8].strip().lower()
    firstatomline = 9
    
  if label[0] == 'c' or label[0] == 'k': fractional = 0
  elif label[0] == 'd': fractional = 1
  else: error("POSCAR file Cartesian or Direct line is missing")

  atomcoords = lines[firstatomline:firstatomline + natoms]

  xyz = []
  for line in atomcoords:
    coords = line.split()
    x = float(coords[0])
    y = float(coords[1])
    z = float(coords[2])
    if fractional == 1:
      xnew = x*avec[0] + y*bvec[0] + z*cvec[0]
      ynew = x*avec[1] + y*bvec[1] + z*cvec[1]
      znew = x*avec[2] + y*bvec[2] + z*cvec[2]
      x,y,z = xnew,ynew,znew
    else:
      x *= scale
      y *= scale
      z *= scale
    xyz.append((x,y,z))

  # write LAMMPS data file

  d = data()

  d.title = "LAMMPS data file from VASP file %s\n" % vaspfile
  d.headers["atoms"] = natoms
  d.headers["atom types"] = ntypes
  d.headers["xlo xhi"] = (0.0,avec[0])
  d.headers["ylo yhi"] = (0.0,bvec[1])
  d.headers["zlo zhi"] = (0.0,cvec[2])
  if (bvec[0] != 0.0) or (cvec[0] != 0.0) or (cvec[1] != 0.0):
    d.headers["xy xz yz"] = (bvec[0],cvec[0],cvec[1])

  id = 0
  itype = 0
  previous = 0
  
  atomlines = []
  for x,y,z in xyz:
    if id == counts[itype] + previous:
      previous += counts[itype]
      itype += 1
    line = "%d %d %s %s %s\n" % (id+1,itype+1,str(x),str(y),str(z))
    atomlines.append(line)
    id += 1

  d.sections["Atoms"] = atomlines

  d.write(datafile)

# ---------------------------
# main program
# ---------------------------

import sys,os,glob
from data import data

if (len(sys.argv) != 3): error()

vaspfile = sys.argv[1]
datafile = sys.argv[2]

vflag = os.path.isdir(vaspfile)
dflag = os.path.isdir(datafile)

if (vflag and not dflag) or (not vflag and dflag):
  error("Both comand-line args must be files or both must be dirs")
  sys.exit()
  
# if both are files, vaspfile is input, datafile is output

if not vflag and not dflag:
  print "Converting",vaspfile,"to",datafile
  convert(vaspfile,datafile)

# if both are dirs, then process all files in POSCAR dir

if vflag and dflag:
  vaspfiles = glob.glob("%s/*" % vaspfile)
  for infile in vaspfiles:
    outfile = ("%s/data." % datafile) + os.path.basename(infile)
    print "Converting",infile,"to",outfile
    convert(infile,outfile)
