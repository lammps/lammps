#!/usr/bin/env python

# Install.py tool to build the CSlib library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess
sys.path.append('..')
from install_helpers import error,get_cpus,fullpath,which

# help message

help = """
Syntax from src dir: make lib-message args="-m"
                 or: make lib-message args="-s -z"
Syntax from lib dir: python Install.py -m
                 or: python Install.py -s -z 

specify zero or more options, order does not matter

  -m = parallel build of CSlib library
  -s = serial build of CSlib library
  -z = build CSlib library with ZMQ socket support, default = no ZMQ support

Example:

make lib-message args="-m -z"   # build parallel CSlib with ZMQ support
make lib-message args="-s"   # build serial CSlib with no ZMQ support
"""

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error(help=help)

mpiflag = False
serialflag = False
zmqflag = False

iarg = 0
while iarg < nargs:
  if args[iarg] == "-m":
    mpiflag = True
    iarg += 1
  elif args[iarg] == "-s":
    serialflag = True
    iarg += 1
  elif args[iarg] == "-z":
    zmqflag = True
    iarg += 1
  else: error(help=help)

if (not mpiflag and not serialflag):
  error("Must use either -m or -s flag")

if (mpiflag and serialflag):
  error("Cannot use -m and -s flag at the same time")

# build CSlib
# copy resulting lib to cslib/src/libmessage.a
# copy appropriate Makefile.lammps.* to Makefile.lammps

print("Building CSlib ...")
srcdir = fullpath("./cslib/src")

if mpiflag and zmqflag:
  cmd = "cd %s; make lib_parallel" % srcdir
elif mpiflag and not zmqflag:
  cmd = "cd %s; make lib_parallel zmq=no" % srcdir
elif not mpiflag and zmqflag:
  cmd = "cd %s; make lib_serial" % srcdir
elif not mpiflag and not zmqflag:
  cmd = "cd %s; make lib_serial zmq=no" % srcdir
  
print(cmd)
try:
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

if mpiflag: cmd = "cd %s; cp libcsmpi.a libmessage.a" % srcdir
else: cmd = "cd %s; cp libcsnompi.a libmessage.a" % srcdir
print(cmd)
txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
print(txt.decode('UTF-8'))

if zmqflag: cmd = "cp Makefile.lammps.zmq Makefile.lammps"
else: cmd = "cp Makefile.lammps.nozmq Makefile.lammps"
print(cmd)
txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
print(txt.decode('UTF-8'))
