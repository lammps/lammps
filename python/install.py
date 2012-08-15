#!/usr/local/bin/python

# copy LAMMPS shared library src/liblammps.so and lammps.py to system dirs
# Syntax: python install.py [libdir] [pydir]
#         libdir = target dir for src/liblammps.so, default = /usr/local/lib
#         pydir = target dir for lammps.py, default = Python site-packages dir

import sys,commands

if len(sys.argv) > 3:
  print "Syntax: python install.py [libdir] [pydir]"
  sys.exit()

if len(sys.argv) >= 2: libdir = sys.argv[1]
else: libdir = "/usr/local/lib"

if len(sys.argv) == 3: pydir = sys.argv[2]
else:
  paths = sys.path
  for i,path in enumerate(paths):
    index = path.rfind("site-packages")
    if index < 0: continue
    if index == len(path) - len("site-packages"): break
  pydir = paths[i]

str = "cp ../src/liblammps.so %s" % libdir
print str
outstr = commands.getoutput(str)
if len(outstr.strip()): print outstr

str = "cp ../python/lammps.py %s" % pydir
print str
outstr = commands.getoutput(str)
if len(outstr.strip()): print outstr

