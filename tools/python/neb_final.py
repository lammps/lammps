#!/usr/bin/env python

# make new dump file from final snapshots from multiple NEB replica dumps
# Syntax: neb_final.py -switch arg(s) -switch arg(s) ...
#         -o outfile = new dump file
#            snapshots numbered 1 to M = # of replicas
#         -r dump1 dump2 ... = replica dump files of NEB atoms
#            must be in correct sequence
#         -b dumpfile = background atoms (optional)
#            last snapshot in this file used as static non-NEB atoms

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.insert(1,path)
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
  sys.exit("Syntax: neb_final.py -o outfile -b backfile -r dump1 dump2 ...")

if os.path.exists(outfile): os.remove(outfile)

# nback = additional atoms in each snapshot

if backfile:
  back = dump(backfile)
  t = back.time()
  back.tselect.one(t[-1])
  nback = back.snaps[-1].nselect
else: nback = 0

# write out each snapshot
# time = replica #
# natoms = ntotal, by overwriting nselect
# add background atoms if requested

n = 1
for file in rfiles:
  neb = dump(file)
  t = neb.time()
  neb.tselect.one(t[-1])
  hold = neb.snaps[-1].nselect
  neb.snaps[-1].nselect = hold + nback
  neb.snaps[-1].time = n
  neb.write(outfile,1,1)
  neb.snaps[-1].nselect = hold
  if backfile: back.write(outfile,0,1)
  n += 1
