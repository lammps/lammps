#!/usr/bin/env python

"""
Install.py tool to build the GPU library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import get_cpus

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# help message

HELP = """
Syntax from src dir: make lib-gpu args="-m machine -h hdir -a arch -p precision -e esuffix -b -o osuffix"
Syntax from lib dir: python Install.py -m machine -h hdir -a arch -p precision -e esuffix -b -o osuffix

specify one or more options, order does not matter

copies an existing Makefile.machine in lib/gpu to Makefile.auto
optionally edits these variables in Makefile.auto:
  CUDA_HOME, CUDA_ARCH, CUDA_PRECISION, EXTRAMAKE
optionally uses Makefile.auto to build the GPU library -> libgpu.a
  and to copy a Makefile.lammps.esuffix -> Makefile.lammps
optionally copies Makefile.auto to a new Makefile.osuffix

See lib/gpu/README and the LAMMPS manual for more information
on which settings to use and how to build.

Examples:

make lib-gpu args="-b"      # build GPU lib with default Makefile.linux
make lib-gpu args="-m xk7 -p single -o xk7.single"      # create new Makefile.xk7.single, altered for single-precision
make lib-gpu args="-m mpi -a sm_35 -p single -o mpi.mixed -b" # create new Makefile.mpi.mixed, also build GPU lib with these settings
"""

# parse and process arguments

parser.add_argument("-b", "--build", action="store_true",
                    help="build the GPU library from scratch from a customized Makefile.auto")
parser.add_argument("-m", "--machine", default='linux',
                    help="suffix of Makefile.machine used as base for customizing Makefile.auto")
parser.add_argument("-a", "--arch", default='sm_30',
                    choices=['sm_12', 'sm_13', 'sm_20', 'sm_21', 'sm_30', 'sm_35', 'sm_37',
                             'sm_50', 'sm_52', 'sm_60', 'sm_61', 'sm_70', 'sm_75'],
                    help="set GPU architecture and instruction set (default: 'sm_30')")
parser.add_argument("-p", "--precision", default='mixed', choices=['single', 'mixed', 'double'],
                    help="set GPU kernel precision mode (default: mixed)")
parser.add_argument("-e", "--extramake", default='standard',
                    help="set EXTRAMAKE variable in Makefile.auto to Makefile.lammps.<extramake>")
parser.add_argument("-c", "--cuda",
                    help="set CUDA_HOME variable in Makefile.auto. Will be used if $CUDA_HOME environment variable is not set")
parser.add_argument("-o", "--output",
                    help="if set, copy final Makefile.auto to Makefile.<output> for later re-use")

args = parser.parse_args()

# print help message and exit, if neither build nor output options are given
if not args.build and not args.output:
  parser.print_help()
  sys.exit(HELP)

hflag = 0
eflag = 0
makeflag = 0
outflag = 0

if args.build:
  makeflag = 1

isuffix = args.machine
arch = args.arch

if args.precision == "double":
  precstr = "-D_DOUBLE_DOUBLE"
elif args.precision == "mixed":
  precstr = "-D_SINGLE_DOUBLE"
else:
  precstr = "-D_SINGLE_SINGLE"

lmpsuffix = args.extramake

if args.cuda:
  hflag = 1
  hdir = args.cuda

if args.output:
  outflag = 1
  osuffix = args.output

# create Makefile.auto
# reset EXTRAMAKE, CUDA_HOME, CUDA_ARCH, CUDA_PRECISION if requested

if not os.path.exists("Makefile.%s" % isuffix):
  sys.exit("lib/gpu/Makefile.%s does not exist" % isuffix)

lines = open("Makefile.%s" % isuffix, 'r').readlines()
fp = open("Makefile.auto", 'w')

for line in lines:
  words = line.split()
  if len(words) != 3:
    fp.write(line)
    continue

  if hflag and words[0] == "CUDA_HOME" and words[1] == '=':
    line = line.replace(words[2], hdir)
  if words[0] == "CUDA_ARCH" and words[1] == '=':
    line = line.replace(words[2], "-arch=%s" % arch)
  if words[0] == "CUDA_PRECISION" and words[1] == '=':
    line = line.replace(words[2], precstr)
  if eflag and words[0] == "EXTRAMAKE" and words[1] == '=':
    line = line.replace(words[2], "Makefile.lammps.%s" % lmpsuffix)

  fp.write(line)
fp.close()

# perform make
# make operations copies EXTRAMAKE file to Makefile.lammps

if makeflag:
  print("Building libgpu.a ...")
  if os.path.exists("libgpu.a"):
    os.remove("libgpu.a")
  n_cpus = get_cpus()
  cmd = "make -f Makefile.auto clean; make -f Makefile.auto -j%d" % n_cpus
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

  if not os.path.exists("libgpu.a"):
    sys.exit("Build of lib/gpu/libgpu.a was NOT successful")
  if not os.path.exists("Makefile.lammps"):
    sys.exit("lib/gpu/Makefile.lammps was NOT created")

# copy new Makefile.auto to Makefile.osuffix

if outflag:
  print("Creating new Makefile.%s" % osuffix)
  shutil.copyfile("Makefile.auto", "Makefile.%s" % osuffix)
