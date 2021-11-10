#!/usr/bin/env python

"""
Installer script to install the LAMMPS python package and the corresponding
shared library into either the system-wide site-packages tree, or - failing
that - into the corresponding user tree. Called from the 'install-python'
build target in the conventional and CMake based build systems
"""

# copy LAMMPS shared library and lammps package to system dirs

from __future__ import print_function
import sys,os,shutil,time
from argparse import ArgumentParser

parser = ArgumentParser(prog='install.py',
                        description='LAMMPS python package installer script')

parser.add_argument("-p", "--package", required=True,
                    help="path to the LAMMPS Python package")
parser.add_argument("-l", "--lib", required=True,
                    help="path to the compiled LAMMPS shared library")
parser.add_argument("-v", "--version", required=True,
                    help="path to the LAMMPS version.h header file")

parser.add_argument("-d","--dir",
                    help="Legacy custom installation folder selection for package and library")

args = parser.parse_args()

# validate arguments and make paths absolute

if args.package:
  if not os.path.exists(args.package):
    print( "ERROR: LAMMPS package %s does not exist" % args.package)
    parser.print_help()
    sys.exit(1)
  else:
    args.package = os.path.abspath(args.package)

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

# if a custom directory is given, we copy the files directly
# without any special processing or additional steps to that folder

if args.dir:
  print("Copying LAMMPS Python package to custom folder %s" % args.dir)
  try:
    shutil.copytree(args.package, os.path.join(args.dir,'lammps'))
  except shutil.Error:
    pass # fail silently

  print("Copying LAMMPS shared library to custom folder %s" % args.dir)
  try:
    shutil.copyfile(args.lib, os.path.join(args.dir,os.path.basename(args.lib)))
  except shutil.Error:
    pass # fail silently

  sys.exit()

# extract LAMMPS version string from header
# and convert to python packaging compatible version
def get_lammps_version(header):
    with open(header, 'r') as f:
        line = f.readline()
        start_pos = line.find('"')+1
        end_pos = line.find('"', start_pos)
        t = time.strptime("".join(line[start_pos:end_pos].split()), "%d%b%Y")
        return "{}.{}.{}".format(t.tm_year,t.tm_mon,t.tm_mday)

verstr = get_lammps_version(args.version)

print("Installing LAMMPS Python package version %s into site-packages folder" % verstr)

# we need to switch to the folder of the python package
os.chdir(os.path.dirname(args.package))

from distutils.core import setup
from distutils.sysconfig import get_python_lib
import site
from sys import version_info

if version_info.major >= 3:
    pkgs = ['lammps', 'lammps.mliap']
else:
    pkgs = ['lammps']

#Arguments common to global or user install -- everything but data_files
setup_kwargs= dict(name="lammps",
        version=verstr,
        author="Steve Plimpton",
        author_email="sjplimp@sandia.gov",
        url="https://www.lammps.org",
        description="LAMMPS Molecular Dynamics Python package",
        license="GPL",
        packages=pkgs,
        )

tryuser=False
try:
  sys.argv = ["setup.py","install"]    # as if had run "python setup.py install"
  setup_kwargs['data_files']=[(os.path.join(get_python_lib(), 'lammps'), [args.lib])]
  setup(**setup_kwargs)
except:                                # lgtm [py/catch-base-exception]
  tryuser=True
  print ("Installation into global site-packages folder failed.\nTrying user folder %s now." % site.USER_SITE)

if tryuser:
  try:
    sys.argv = ["setup.py","install","--user"]    # as if had run "python setup.py install --user"
    setup_kwargs['data_files']=[(os.path.join(site.USER_SITE, 'lammps'), [args.lib])]
    setup(**setup_kwargs)
  except:                              # lgtm [py/catch-base-exception]
    print("Installation into user site package folder failed.")
