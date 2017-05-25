#!/usr/bin/env python

# Install.py tool to build the GPU library
# used to automate the steps described in the README file in this dir

import sys,os,re,commands

# help message

help = """
Syntax: python Install.py -i isuffix -h hdir -a arch -p precision -e esuffix -m -o osuffix
  specify one or more options, order does not matter
  copies an existing Makefile.isuffix in lib/gpu to Makefile.auto 
  optionally edits these variables in Makefile.auto:
    CUDA_HOME, CUDA_ARCH, CUDA_PRECISION, EXTRAMAKE
  optionally uses Makefile.auto to build the GPU library -> libgpu.a
    and to copy a Makefile.lammps.esuffix -> Makefile.lammps
  optionally copies Makefile.auto to a new Makefile.osuffix

  -i = use Makefile.isuffix as starting point, copy to Makefile.auto
       default isuffix = linux
  -h = set CUDA_HOME variable in Makefile.auto to hdir
       hdir = path to NVIDIA Cuda software, e.g. /usr/local/cuda
  -a = set CUDA_ARCH variable in Makefile.auto to arch
       use arch = ?? for K40 (Tesla)
       use arch = 37 for dual K80 (Tesla)
       use arch = 60 for P100 (Pascal)
  -p = set CUDA_PRECISION variable in Makefile.auto to precision
       use precision = double or mixed or single
  -e = set EXTRAMAKE variable in Makefile.auto to Makefile.lammps.esuffix
  -m = make the GPU library using Makefile.auto
       first performs a "make clean"
       produces libgpu.a if successful
       also copies EXTRAMAKE file -> Makefile.lammps
         -e can set which Makefile.lammps.esuffix file is copied
  -o = copy final Makefile.auto to Makefile.osuffix
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

isuffix = "linux"
hflag = aflag = pflag = eflag = 0
makeflag = 0
outflag = 0

iarg = 0
while iarg < nargs:
  if args[iarg] == "-i":
    if iarg+2 > nargs: error()
    isuffix = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-h":
    if iarg+2 > nargs: error()
    hflag = 1
    hdir = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-a":
    if iarg+2 > nargs: error()
    aflag = 1
    arch = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error()
    pflag = 1
    precision = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-e":
    if iarg+2 > nargs: error()
    eflag = 1
    lmpsuffix = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-m":
    makeflag = 1
    iarg += 1
  elif args[iarg] == "-o":
    if iarg+2 > nargs: error()
    outflag = 1
    osuffix = args[iarg+1]
    iarg += 2
  else: error()

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
    print >>fp,line,
    continue
  
  if hflag and words[0] == "CUDA_HOME" and words[1] == '=':
    line = line.replace(words[2],hdir)
  if aflag and words[0] == "CUDA_ARCH" and words[1] == '=':
    line = line.replace(words[2],"-arch=sm_%s" % arch)
  if pflag and words[0] == "CUDA_PRECISION" and words[1] == '=':
    line = line.replace(words[2],precstr)
  if eflag and words[0] == "EXTRAMAKE" and words[1] == '=':
    line = line.replace(words[2],"Makefile.lammps.%s" % lmpsuffix)
    
  print >>fp,line,

fp.close()

# perform make
# make operations copies EXTRAMAKE file to Makefile.lammps

if makeflag:
  print "Building libgpu.a ..."
  cmd = "rm -f libgpu.a"
  commands.getoutput(cmd)
  cmd = "make -f Makefile.auto clean; make -f Makefile.auto"
  commands.getoutput(cmd)
  if not os.path.exists("libgpu.a"):
    error("Build of lib/gpu/libgpu.a was NOT successful")
  if not os.path.exists("Makefile.lammps"):
    error("lib/gpu/Makefile.lammps was NOT created")

# copy new Makefile.auto to Makefile.osuffix

if outflag:
  print "Creating new Makefile.%s" % osuffix
  cmd = "cp Makefile.auto Makefile.%s" % osuffix
  commands.getoutput(cmd)
