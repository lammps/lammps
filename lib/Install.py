#!/usr/bin/env python

"""
Install.py tool to do a generic build of a library
soft linked to by many of the lib/Install.py files
used to automate the steps described in the corresponding lib/README
"""

from __future__ import print_function
import sys, os, subprocess
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import get_cpus, fullpath

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

HELP = """
Syntax from src dir: make lib-libname args="-m machine -e suffix"
Syntax from lib dir: python Install.py -m machine -e suffix

libname = name of lib dir (e.g. atc, h5md, meam, poems, etc)
specify -m and optionally -e, order does not matter

Examples:

make lib-poems args="-m serial" # build POEMS lib with same settings as in the serial Makefile in src
make lib-colvars args="-m mpi"  # build COLVARS lib with same settings as in the mpi Makefile in src
make lib-meam args="-m ifort"   # build MEAM lib with custom Makefile.ifort (using Intel Fortran)
"""

# parse and process arguments

parser.add_argument("-m", "--machine",
                    help="suffix of a <libname>/Makefile.* file used for compiling this library")
parser.add_argument("-e", "--extramake",
                    help="set EXTRAMAKE variable in <libname>/Makefile.<machine> to Makefile.lammps.<extramake>")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.machine and not args.extramake:
  parser.print_help()
  sys.exit(HELP)

machine = args.machine
extraflag = args.extramake
if extraflag:
    suffix = args.extramake
else:
    suffix = 'empty'

# set lib from working dir

cwd = fullpath('.')
lib = os.path.basename(cwd)

# create Makefile.auto as copy of Makefile.machine
# reset EXTRAMAKE if requested

if not os.path.exists("Makefile.%s" % machine):
  sys.exit("lib/%s/Makefile.%s does not exist" % (lib, machine))

lines = open("Makefile.%s" % machine, 'r').readlines()
fp = open("Makefile.auto", 'w')

has_extramake = False
for line in lines:
  words = line.split()
  if len(words) == 3 and words[0] == "EXTRAMAKE" and words[1] == '=':
    has_extramake = True
    if extraflag:
      line = line.replace(words[2], "Makefile.lammps.%s" % suffix)
  fp.write(line)

fp.close()

# make the library via Makefile.auto optionally with parallel make
n_cpus = get_cpus()

print("Building lib%s.a ..." % lib)
cmd = "make -f Makefile.auto clean; make -f Makefile.auto -j%d" % n_cpus
try:
  txt = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
  print("Make failed with:\n %s" % e.output.decode('UTF-8'))
  sys.exit(1)

if os.path.exists("lib%s.a" % lib):
  print("Build was successful")
else:
  sys.exit("Build of lib/%s/lib%s.a was NOT successful" % (lib, lib))

if has_extramake and not os.path.exists("Makefile.lammps"):
  print("WARNING: lib/%s/Makefile.lammps was NOT created" % lib)
