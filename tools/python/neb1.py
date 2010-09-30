#!/usr/local/bin/python

# combine multiple NEB dump files into 1 dump file
# Syntax: neb1.py dfinal itype dfile1 dfile2 ...
#         dfinal = new combined dump file
#         itype = atoms of this type are added to each snapshot
#         dfile1, dfile2, ... = NEB output dump files from each replica

# assumes all N input files have same # of snapshots
# for each snapshot time:
#   all dfile1 atoms are added to dfinal
#   only atoms of itype in dfile2,...,dfileN are added to dfinal

import sys,os
from dump import dump
if not globals().has_key("argv"): argv = sys.argv

if len(argv) < 5:
  print "Syntax: neb1.py dfinal itype dfile1 dfile2 ..."
  sys.exit()

dfinal = argv[1]
itype = int(argv[2])
files = argv[3:]

tstr = "$type == %d" % itype

ntotal = 0
d = []
for file in files:
  one = dump(file)
  one.sort()
  if d:
    one.aselect.test(tstr)
    nnew = one.snaps[0].nselect
    idvec = range(ntotal+1,ntotal+nnew+1)
    one.setv("id",idvec)
  ntotal += one.snaps[0].nselect
  d.append(one)

if os.path.exists(dfinal): os.remove(dfinal)

times = d[0].time()
for time in times:
  d[0].tselect.one(time)
  i = d[0].findtime(time)
  hold = d[0].snaps[i].nselect
  d[0].snaps[i].nselect = ntotal
  d[0].write(dfinal,1,1)
  d[0].snaps[i].nselec = hold
  for one in d[1:]:
    one.tselect.one(time)
    one.aselect.test(tstr)
    one.write(dfinal,0,1)
