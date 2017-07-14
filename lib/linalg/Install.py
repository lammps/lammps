#!/usr/bin/env python

# install.py tool to do build of the linear algebra library
# used to automate the steps described in the README file in this dir

import sys,commands,os

# help message

help = """
Syntax from src dir: make lib-linalg args="-m machine"
Syntax from lib dir: python Install.py -m machine

  -m = peform a clean followed by "make -f Makefile.machine"
       machine = suffix of a lib/Makefile.* file

Example:

make lib-linalg args="-m gfortran"   # build with GNU Fortran compiler
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

iarg = 0
while iarg < nargs:
  if args[iarg] == "-m":
    if iarg+2 > nargs: error()
    machine = args[iarg+1]
    iarg += 2  
  else: error()

# set lib from working dir

cwd = os.getcwd()
lib = os.path.basename(cwd)

# make the library

print "Building lib%s.a ..." % lib
cmd = "make -f Makefile.%s clean; make -f Makefile.%s" % (machine,machine)
txt = commands.getoutput(cmd)
print txt

if os.path.exists("lib%s.a" % lib): print "Build was successful"
else: error("Build of lib/%s/lib%s.a was NOT successful" % (lib,lib))
