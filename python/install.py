#!/usr/local/bin/python

# copy lammps.py and LAMMPS shared library src/liblmp.so to
# Python site-packages dir

import sys,commands

# find this Python's site-packages dir in sys.path list

paths = sys.path
for i,path in enumerate(paths):
  index = path.rfind("site-packages")
  if index < 0: continue
  if index == len(path) - len("site-packages"): break
sitedir = paths[i]

str = "cp ../python/lammps.py %s" % sitedir
print str
outstr = commands.getoutput(str)
if len(outstr.strip()): print outstr

str = "cp ../src/liblmp.so %s" % sitedir
print str
outstr = commands.getoutput(str)
if len(outstr.strip()): print outstr
