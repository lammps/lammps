#!/usr/bin/env python

"""
Installer script to install the LAMMPS python module and the corresponding
shared library into either the system-wide site-packages tree, or - failing
that - into the corresponding user tree. Called from the 'install-python'
build target in the conventional and CMake based build systems
"""

# copy LAMMPS shared library and lammps.py to system dirs

from __future__ import print_function
import sys,os,shutil
from argparse import ArgumentParser

parser = ArgumentParser(prog='install.py',
                        description='LAMMPS python module installer script')

parser.add_argument("-m", "--module", required=True,
                    help="path to the source of the LAMMPS Python module")
parser.add_argument("-l", "--lib", required=True,
                    help="path to the compiled LAMMPS shared library")
parser.add_argument("-v", "--version", required=True,
                    help="path to the LAMMPS version.h header file")

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-d","--dir",
                    help="Legacy custom installation folder for module and library")
pgroup.add_argument("-p","--prefix",
                    help="Installation prefix for module and library")

args = parser.parse_args()

# validate arguments and make paths absolute

if args.module:
  if not os.path.exists(args.module):
    print( "ERROR: LAMMPS module file %s does not exist" % args.module)
    parser.print_help()
    sys.exit(1)
  else:
    args.module = os.path.abspath(args.module)

if args.lib:
  if not os.path.exists(args.lib):
    print( "ERROR: LAMMPS shared library %s does not exist" % args.lib)
    parser.print_help()
    sys.exit(1)
  else:
    args.lib = os.path.abspath(args.lib)

if args.version:
  if not os.path.exists(args.version):
    print( "ERROR: LAMMPS version header file %s does not exist" % args.version)
    parser.print_help()
    sys.exit(1)
  else:
    args.version = os.path.abspath(args.version)

if args.dir:
  if not os.path.isdir(args.dir):
    print( "ERROR: Installation folder %s does not exist" % args.dir)
    parser.print_help()
    sys.exit(1)
  else:
    args.dir = os.path.abspath(args.dir)

if args.prefix:
  if not os.path.isdir(args.prefix):
    print( "ERROR: Installation prefix folder %s does not exist" % args.prefix)
    parser.print_help()
    sys.exit(1)
  else:
    args.prefix = os.path.abspath(args.prefix)

# if a custom directory is given, we copy the files directly
# without any special processing or additional steps to that folder

if args.dir:
  print("Copying LAMMPS Python module to custom folder %s" % args.dir)
  try:
    shutil.copyfile(args.module, os.path.join(args.dir,'lammps.py'))
  except shutil.Error:
    pass # fail silently

  print("Copying LAMMPS shared library to custom folder %s" % args.dir)
  try:
    shutil.copyfile(args.lib, os.path.join(args.dir,os.path.basename(args.lib)))
  except shutil.Error:
    pass # fail silently

  sys.exit()

# extract version string from header
fp = open(args.version,'r')
txt=fp.read().split('"')[1].split()
verstr=txt[0]+txt[1]+txt[2]
fp.close()

print("Installing LAMMPS Python module version %s into site-packages folder" % verstr)

# we need to switch to the folder of the python module
os.chdir(os.path.dirname(args.module))

from distutils.core import setup
from distutils.sysconfig import get_python_lib
import site
tryuser=False

try:
  if args.prefix:
    sys.argv = ["setup.py","install","--prefix=%s" % args.prefix]    # as if had run "python setup.py install --prefix=XXX"
  else:
    sys.argv = ["setup.py","install"]    # as if had run "python setup.py install"
  setup(name = "lammps",
        version = verstr,
        author = "Steve Plimpton",
        author_email = "sjplimp@sandia.gov",
        url = "https://lammps.sandia.gov",
        description = "LAMMPS Molecular Dynamics Python module",
        license = "GPL",
        py_modules = ["lammps"],
        data_files = [(get_python_lib(prefix=args.prefix), [args.lib])])
except:
  tryuser=True
  print ("Installation into global site-packages folder failed.\nTrying user folder %s now." % site.USER_SITE)

if tryuser:
  try:
    sys.argv = ["setup.py","install","--user"]    # as if had run "python setup.py install --user"
    setup(name = "lammps",
          version = verstr,
          author = "Steve Plimpton",
          author_email = "sjplimp@sandia.gov",
          url = "https://lammps.sandia.gov",
          description = "LAMMPS Molecular Dynamics Python module",
          license = "GPL",
          py_modules = ["lammps"],
          data_files = [(site.USER_SITE, [args.lib])])
  except:
    print("Installation into user site package folder failed.")
