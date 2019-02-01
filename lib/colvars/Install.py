#!/usr/bin/env python

# Install.py tool to do automate build of Colvars

from __future__ import print_function
import sys,os,subprocess
sys.path.append('..')
from install_helpers import get_cpus

from argparse import ArgumentParser

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# help message

help = """
Syntax from src dir: make lib-colvars args="-m machine -e suffix"
Syntax from lib/colvars dir: python Install.py -m machine -e suffix

specify -m and optionally -e, order does not matter

  -m = delete all existing objects, followed by "make -f Makefile.machine"
       machine = suffix of a lib/colvars/Makefile.* or of a
         src/MAKE/MACHINES/Makefile.* file
  -e = set EXTRAMAKE variable in Makefile.machine to Makefile.lammps.suffix
       does not alter existing Makefile.machine

Examples:

make lib-colvars args="-m mpi"     # build COLVARS lib with default mpi compiler wrapper
"""

# set lib from working dir

cwd = os.getcwd()
lib = os.path.basename(cwd)

# parse and process arguments

parser.add_argument("-m", "--machine",
                    help="suffix of a <libname>/Makefile.* or of a src/MAKE/MACHINES/Makefile.* file used for compiling this library")
parser.add_argument("-e", "--extramake",
                    help="set EXTRAMAKE variable in <libname>/Makefile.<machine> to Makefile.lammps.<extramake>")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.machine and not args.extramake:
  parser.print_help()
  sys.exit(help)

machine = args.machine
extraflag = args.extramake != None
suffix = args.extramake

def get_lammps_machine_flags(machine):
  """Parse Makefile.machine from LAMMPS, return dictionary of compiler flags"""
  if not os.path.exists("../../src/MAKE/MACHINES/Makefile.%s" % machine):
    sys.exit("Cannot locate src/MAKE/MACHINES/Makefile.%s" % machine)
  lines = open("../../src/MAKE/MACHINES/Makefile.%s" % machine,
               'r').readlines()
  machine_flags = {}
  for line in lines:
    line = line.partition('#')[0]
    line = line.rstrip()
    words = line.split()
    if (len(words) > 2):
      if ((words[0] == 'CC') or (words[0] == 'CCFLAGS') or
          (words[0] == 'SHFLAGS') or (words[0] == 'ARCHIVE') or
          (words[0] == 'ARFLAGS') or (words[0] == 'SHELL')):
        machine_flags[words[0]] = ' '.join(words[2:])
  return machine_flags

def gen_colvars_makefile_machine(machine, machine_flags):
  """Generate Makefile.machine for Colvars given the compiler flags"""
  machine_makefile = open("Makefile.%s" % machine, 'w')
  machine_makefile.write('''# -*- makefile -*- to build Colvars module with %s

COLVARS_LIB = libcolvars.a
COLVARS_OBJ_DIR =

CXX =		%s
CXXFLAGS =	%s %s
AR =		%s
ARFLAGS =	%s
SHELL =		%s

include Makefile.common

.PHONY: default clean

default: $(COLVARS_LIB) Makefile.lammps

clean:
	-rm -f $(COLVARS_OBJS) $(COLVARS_LIB)
''' % (machine, machine_flags['CC'],
       machine_flags['CCFLAGS'], machine_flags['SHFLAGS'] ,
       machine_flags['ARCHIVE'], machine_flags['ARFLAGS'],
       machine_flags['SHELL']))

if not os.path.exists("Makefile.%s" % machine):
  machine_flags = get_lammps_machine_flags(machine)
  gen_colvars_makefile_machine(machine, machine_flags)
if not os.path.exists("Makefile.%s" % machine):
  sys.exit("lib/%s/Makefile.%s does not exist" % (lib,machine))

# create Makefile.auto as copy of Makefile.machine
# reset EXTRAMAKE if requested

lines = open("Makefile.%s" % machine,'r').readlines()
fp = open("Makefile.auto",'w')
for line in lines:
  words = line.split()
  if len(words) == 3 and extraflag and \
        words[0] == "EXTRAMAKE" and words[1] == '=':
    line = line.replace(words[2],"Makefile.lammps.%s" % suffix)
  fp.write(line)
fp.close()

# make the library via Makefile.auto optionally with parallel make

n_cpus = get_cpus()

print("Building lib%s.a ..." % lib)
cmd = ["make -f Makefile.auto clean; make -f Makefile.auto -j%d" % n_cpus]
try:
  txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True);
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
  print("Make failed with:\n %s" % e.output.decode('UTF-8'))
  sys.exit(1)

if os.path.exists("lib%s.a" % lib): print("Build was successful")
else: sys.exit("Build of lib/%s/lib%s.a was NOT successful" % (lib,lib))
if not os.path.exists("Makefile.lammps"):
  print("lib/%s/Makefile.lammps was NOT created" % lib)
