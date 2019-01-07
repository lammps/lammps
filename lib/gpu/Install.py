#!/usr/bin/env python

# Install.py tool to build the GPU library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,subprocess
sys.path.append('..')
from install_helpers import error,get_cpus

# help message

help = """
Syntax from src dir: make lib-gpu args="-m machine -h hdir -a arch -p precision -e esuffix -b -o osuffix"
Syntax from lib dir: python Install.py -m machine -h hdir -a arch -p precision -e esuffix -b -o osuffix

specify one or more options, order does not matter

copies an existing Makefile.machine in lib/gpu to Makefile.auto
optionally edits these variables in Makefile.auto:
  CUDA_HOME, CUDA_ARCH, CUDA_PRECISION, EXTRAMAKE
optionally uses Makefile.auto to build the GPU library -> libgpu.a
  and to copy a Makefile.lammps.esuffix -> Makefile.lammps
optionally copies Makefile.auto to a new Makefile.osuffix

  -m = use Makefile.machine as starting point, copy to Makefile.auto
       default machine = linux
       default for -h, -a, -p, -e settings are those in -m Makefile
  -h = set CUDA_HOME variable in Makefile.auto to hdir
       hdir = path to NVIDIA Cuda software, e.g. /usr/local/cuda
  -a = set CUDA_ARCH variable in Makefile.auto to arch
       use arch = sm_20 for Fermi (C2050/C2070, deprecated as of CUDA 8.0)
                        or GeForce GTX 580 or similar
       use arch = sm_30 for Kepler (K10)
       use arch = sm_35 for Kepler (K40) or GeForce GTX Titan or similar
       use arch = sm_37 for Kepler (dual K80)
       use arch = sm_60 for Pascal (P100)
       use arch = sm_70 for Volta
  -p = set CUDA_PRECISION variable in Makefile.auto to precision
       use precision = double or mixed or single
  -e = set EXTRAMAKE variable in Makefile.auto to Makefile.lammps.esuffix
  -b = make the GPU library using Makefile.auto
       first performs a "make clean"
       then produces libgpu.a if successful
       also copies EXTRAMAKE file -> Makefile.lammps
         -e can set which Makefile.lammps.esuffix file is copied
  -o = copy final Makefile.auto to Makefile.osuffix

Examples:

make lib-gpu args="-b"      # build GPU lib with default Makefile.linux
make lib-gpu args="-m xk7 -p single -o xk7.single"      # create new Makefile.xk7.single, altered for single-precision
make lib-gpu args="-m mpi -a sm_35 -p single -o mpi.mixed -b" # create new Makefile.mpi.mixed, also build GPU lib with these settings
"""

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error(help=help)

isuffix = "linux"
hflag = aflag = pflag = eflag = 0
makeflag = 0
outflag = 0

iarg = 0
while iarg < nargs:
  if args[iarg] == "-m":
    if iarg+2 > nargs: error(help=help)
    isuffix = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-h":
    if iarg+2 > nargs: error(help=help)
    hflag = 1
    hdir = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-a":
    if iarg+2 > nargs: error(help=help)
    aflag = 1
    arch = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error(help=help)
    pflag = 1
    precision = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-e":
    if iarg+2 > nargs: error(help=help)
    eflag = 1
    lmpsuffix = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    makeflag = 1
    iarg += 1
  elif args[iarg] == "-o":
    if iarg+2 > nargs: error(help=help)
    outflag = 1
    osuffix = args[iarg+1]
    iarg += 2
  else: error(help=help)

if pflag:
  if precision == "double": precstr = "-D_DOUBLE_DOUBLE"
  elif precision == "mixed": precstr = "-D_SINGLE_DOUBLE"
  elif precision == "single": precstr = "-D_SINGLE_SINGLE"
  else: error("Invalid precision setting")

# create Makefile.auto
# reset EXTRAMAKE, CUDA_HOME, CUDA_ARCH, CUDA_PRECISION if requested

if not os.path.exists("Makefile.%s" % isuffix):
  error("lib/gpu/Makefile.%s does not exist" % isuffix)

lines = open("Makefile.%s" % isuffix,'r').readlines()
fp = open("Makefile.auto",'w')

for line in lines:
  words = line.split()
  if len(words) != 3:
    fp.write(line)
    continue

  if hflag and words[0] == "CUDA_HOME" and words[1] == '=':
    line = line.replace(words[2],hdir)
  if aflag and words[0] == "CUDA_ARCH" and words[1] == '=':
    line = line.replace(words[2],"-arch=%s" % arch)
  if pflag and words[0] == "CUDA_PRECISION" and words[1] == '=':
    line = line.replace(words[2],precstr)
  if eflag and words[0] == "EXTRAMAKE" and words[1] == '=':
    line = line.replace(words[2],"Makefile.lammps.%s" % lmpsuffix)

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
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

  if not os.path.exists("libgpu.a"):
    error("Build of lib/gpu/libgpu.a was NOT successful")
  if not os.path.exists("Makefile.lammps"):
    error("lib/gpu/Makefile.lammps was NOT created")

# copy new Makefile.auto to Makefile.osuffix

if outflag:
  print("Creating new Makefile.%s" % osuffix)
  cmd = "cp Makefile.auto Makefile.%s" % osuffix
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
