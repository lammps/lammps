#!/usr/bin/env python

# make new dump file by combining snapshots from multiple NEB replica dumps
# Syntax: neb_combine.py -switch arg(s) -switch arg(s) ...
#         -o outfile = new dump file
#            each snapshot has NEB atoms from all replicas
#         -r dump1 dump2 ... = replica dump files of NEB atoms
#            can be in any order
#         -b dumpfile = background atoms (optional)
#            first snapshot in this file used as static non-NEB atoms

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.append(path)
from dump import dump

# parse args

outfile = ""
backfile = ""
rfiles = []

argv = sys.argv
iarg = 1
narg = len(argv)
while iarg < narg:
  if argv[iarg] == "-o":
    outfile = argv[iarg+1]
    iarg += 2
  elif argv[iarg] == "-b":
    backfile = argv[iarg+1]
    iarg += 2
  elif argv[iarg] == "-r":
    ilast = iarg + 1
    while ilast < narg and argv[ilast][0] != '-': ilast += 1
    rfiles = argv[iarg+1:ilast]
    iarg = ilast
  else: break

if iarg < narg or not outfile or not rfiles:
  print "Syntax: neb_combine.py -o outfile -b backfile -r dump1 dump2 ..."
  sys.exit()

if os.path.exists(outfile): os.remove(outfile)

# ntotal = total atoms in each snapshot
# reset IDs of atoms in each NEB dump file

ntotal = 0
d = []
for file in rfiles:
  one = dump(file)
  nnew = one.snaps[0].nselect
  idvec = range(ntotal+1,ntotal+nnew+1)
  one.setv("id",idvec)
  ntotal += nnew
  d.append(one)

# nback = additional atoms in each snapshot
# reset IDs of atoms in background file

if backfile:
  back = dump(backfile)
  t = back.time()
  back.tselect.one(t[0])
  nback = back.snaps[0].nselect
  idvec = range(ntotal+1,ntotal+nback+1)
  back.setv("id",idvec)
else: nback = 0
ntotal += nback

# write out each snapshot
# natoms = ntotal, by overwriting nselect
# add background atoms if requested

times = d[0].time()
for time in times:
  d[0].tselect.one(time)
  i = d[0].findtime(time)
  hold = d[0].snaps[i].nselect
  d[0].snaps[i].nselect = ntotal
  d[0].write(outfile,1,1)
  d[0].snaps[i].nselect = hold
  for one in d[1:]:
    one.tselect.one(time)
    one.write(outfile,0,1)
  if backfile: back.write(outfile,0,1)

