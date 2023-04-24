#!/usr/bin/env python

"""
Install.py tool to build the Lepton library
"""

from __future__ import print_function
import sys, os, subprocess
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import get_cpus, fullpath

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS Lepton library build wrapper script")

HELP = """
Syntax from src dir: make lib-lepton args="-m machine"
Syntax from lib dir: python Install.py -m machine

specify -m

Examples:

make lib-lepton args="-m serial" # build Lepton lib with same settings as in the serial Makefile in src
python Install.py -m mpi         # build Lepton lib with same settings as in the mpi Makefile in src
"""

# parse and process arguments

parser.add_argument("-m", "--machine",
                    help="suffix of a <libname>/Makefile.* file used for compiling this library")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.machine:
  parser.print_help()
  sys.exit(HELP)

machine = args.machine

# set lib from working dir

cwd = fullpath('.')
lib = os.path.basename(cwd)

if not os.path.exists("Makefile.%s" % machine):
  sys.exit("lib/%s/Makefile.%s does not exist" % (lib, machine))

# make the library with parallel make
n_cpus = get_cpus()

print("Building lib%s.a ..." % lib)
cmd = "make -f Makefile.%s clean; make -f Makefile.%s -j%d" % (machine, machine, n_cpus)
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

if not os.path.exists("Makefile.lammps"):
  print("WARNING: lib/%s/Makefile.lammps was NOT created" % lib)
