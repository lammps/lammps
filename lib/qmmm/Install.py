#!/usr/bin/env python

# install.py tool to do a generic build of a library
# soft linked to by many of the lib/Install.py files
# used to automate the steps described in the corresponding lib/README

import sys,commands,os

# help message

help = """
Syntax: python Install.py -m machine -e suffix
  specify -m and optionally -e, order does not matter
  -m = peform a clean followed by "make -f Makefile.machine"
       machine = suffix of a lib/Makefile.* file
  -e = set EXTRAMAKE variable in Makefile.machine to Makefile.lammps.suffix
       does not alter existing Makefile.machine
"""

# print error message or help

def error(str=None):
  if not str: print help
  else: print "ERROR",str
  sys.exit()

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

machine = None
extraflag = 0

iarg = 0
while iarg < nargs:
  if args[iarg] == "-m":
    if iarg+2 > nargs: error()
    machine = args[iarg+1]
    iarg += 2  
  elif args[iarg] == "-e":
    if iarg+2 > nargs: error()
    extraflag = 1
    suffix = args[iarg+1]
    iarg += 2  
  else: error()

# set lib from working dir

cwd = os.getcwd()
lib = os.path.basename(cwd)

# create Makefile.auto as copy of Makefile.machine
# reset EXTRAMAKE if requested
  
if not os.path.exists("Makefile.%s" % machine):
  error("lib/%s/Makefile.%s does not exist" % (lib,machine))

lines = open("Makefile.%s" % machine,'r').readlines()
fp = open("Makefile.auto",'w')

for line in lines:
  words = line.split()
  if len(words) == 3 and extraflag and \
        words[0] == "EXTRAMAKE" and words[1] == '=':
    line = line.replace(words[2],"Makefile.lammps.%s" % suffix)
  print >>fp,line,

fp.close()

# make the library via Makefile.auto

print "Building lib%s.a ..." % lib
cmd = "make -f Makefile.auto clean; make -f Makefile.auto"
txt = commands.getoutput(cmd)
print txt

if os.path.exists("lib%s.a" % lib): print "Build was successful"
else: error("Build of lib/%s/lib%s.a was NOT successful" % (lib,lib))
if not os.path.exists("Makefile.lammps"):
  print "lib/%s/Makefile.lammps was NOT created" % lib
