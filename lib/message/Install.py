#!/usr/bin/env python

"""
Install.py tool to build the CSlib library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# help message

HELP = """
Syntax from src dir: make lib-message args="-m"
                 or: make lib-message args="-s -z"
Syntax from lib dir: python Install.py -m
                 or: python Install.py -s -z

Example:

make lib-message args="-m -z"   # build parallel CSlib with ZMQ support
make lib-message args="-s"   # build serial CSlib with no ZMQ support
"""

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-m", "--mpi", action="store_true",
                    help="parallel build of CSlib with MPI")
pgroup.add_argument("-s", "--serial", action="store_true",
                    help="serial build of CSlib")
parser.add_argument("-z", "--zmq", default=False, action="store_true",
                    help="build CSlib with ZMQ socket support, default ()")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.mpi and not args.serial:
  parser.print_help()
  sys.exit(HELP)

mpiflag = args.mpi
serialflag = args.serial
zmqflag = args.zmq

# build CSlib
# copy resulting lib to cslib/src/libmessage.a
# copy appropriate Makefile.lammps.* to Makefile.lammps

print("Building CSlib ...")
srcdir = fullpath(os.path.join("cslib", "src"))

if mpiflag and zmqflag:
  cmd = "make -C %s lib_parallel" % srcdir
elif mpiflag and not zmqflag:
  cmd = "make -C %s lib_parallel zmq=no" % srcdir
elif not mpiflag and zmqflag:
  cmd = "make -C %s lib_serial" % srcdir
elif not mpiflag and not zmqflag:
  cmd = "make -C %s lib_serial zmq=no" % srcdir

print(cmd)
try:
  txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

slb = os.path.join(srcdir, "libcsnompi.a")
if mpiflag:
  slb = os.path.join(srcdir, "libcsmpi.a")
shutil.copyfile(slb, os.path.join(srcdir, "libmessage.a"))

smk = "Makefile.lammps.nozmq"
if zmqflag:
  smk = "Makefile.lammps.zmq"
shutil.copyfile(smk, "Makefile.lammps")
print("Using %s for Makefile.lammps" % smk)
