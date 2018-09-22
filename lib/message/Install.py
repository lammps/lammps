#!/usr/bin/env python

# Install.py tool to build the CSlib library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

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

# print error message or help

def error(str=None):
  if not str: print(help)
  else: print("ERROR",str)
  sys.exit()

# expand to full path name
# process leading '~' or relative path

def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))

def which(program):
  def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath:
    if is_exe(program):
      return program
  else:
    for path in os.environ["PATH"].split(os.pathsep):
      path = path.strip('"')
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file

  return None

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

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
  else: error()

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
txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
print(txt.decode('UTF-8'))

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
