#!/usr/bin/env python

# install.py tool to do a generic build of a library
# soft linked to by many of the lib/Install.py files
# used to automate the steps described in the corresponding lib/README

from __future__ import print_function
import sys,os,subprocess

# help message

help = """
Syntax from src dir: make lib-libname args="-m machine -e suffix"
Syntax from lib dir: python Install.py -m machine -e suffix

libname = name of lib dir (e.g. atc, h5md, meam, poems, etc)
specify -m and optionally -e, order does not matter

  -m = peform a clean followed by "make -f Makefile.machine"
       machine = suffix of a lib/Makefile.* file
  -e = set EXTRAMAKE variable in Makefile.machine to Makefile.lammps.suffix
       does not alter existing Makefile.machine

Examples:

make lib-poems args="-m g++"    # build COLVARS lib with GNU g++ compiler
make lib-meam args="-m ifort"   # build MEAM lib with Intel ifort compiler
"""

# print error message or help

def error(str=None):
  if not str: print(help)
  else: print("ERROR",str)
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

has_extramake = False
for line in lines:
  words = line.split()
  if len(words) == 3 and words[0] == "EXTRAMAKE" and words[1] == '=':
    has_extramake = True
    if extraflag:
      line = line.replace(words[2],"Makefile.lammps.%s" % suffix)
  fp.write(line)

fp.close()

# make the library via Makefile.auto optionally with parallel make

try:
  import multiprocessing
  n_cpus = multiprocessing.cpu_count()
except:
  n_cpus = 1

print("Building lib%s.a ..." % lib)
cmd = "make -f Makefile.auto clean; make -f Makefile.auto -j%d" % n_cpus
txt = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
print(txt.decode('UTF-8'))

if os.path.exists("lib%s.a" % lib): print("Build was successful")
else: error("Build of lib/%s/lib%s.a was NOT successful" % (lib,lib))
if has_extramake and not os.path.exists("Makefile.lammps"):
  print("lib/%s/Makefile.lammps was NOT created" % lib)
