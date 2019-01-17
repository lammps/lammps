#!/usr/bin/env python

# install.py tool to do a generic build of a library
# soft linked to by many of the lib/Install.py files
# used to automate the steps described in the corresponding lib/README

from __future__ import print_function
import sys,os,subprocess
sys.path.append('..')
from install_helpers import error,get_cpus

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
n_cpus = get_cpus()

print("Building lib%s.a ..." % lib)
cmd = "make -f Makefile.auto clean; make -f Makefile.auto -j%d" % n_cpus
try:
  txt = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

if os.path.exists("lib%s.a" % lib): print("Build was successful")
else: error("Build of lib/%s/lib%s.a was NOT successful" % (lib,lib))
if has_extramake and not os.path.exists("Makefile.lammps"):
  print("lib/%s/Makefile.lammps was NOT created" % lib)
