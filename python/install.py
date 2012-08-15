#!/usr/bin/env python

instructions = """copy LAMMPS shared library src/liblammps.so and lammps.py to system dirs
Syntax: python install.py [libdir] [pydir]
    libdir = target dir for src/liblammps.so, default = /usr/local/lib, or the first
        item in LD_LIBRARY_PATH if it doesn't exist.
    pydir = target dir for lammps.py, default = Python site-packages, via distutils."""

import sys, shutil, os

if len(sys.argv) > 3:
  print instructions
  sys.exit()

# verify that our user-specified path is in LD_LIBRARY_PATH
# since if not, the install won't work
  
libdir = "/usr/local/lib"
libpaths = os.environ['LD_LIBRARY_PATH'].split(':')
if not libdir in libpaths:
    libdir = libpaths[0]

pydir = False
try:
    libdir = sys.argv[1]
    pydir = sys.argv[2]
except IndexError:
    pass

# copy the C library into place

shutil.copy('../src/liblammps.so', libdir)

# if user-specified, copy lammps.py into directory
# else invoke setup from Distutils to add to site-packages

if pydir:
    shutil.copy('../python/lammps.py', pydir)
    sys.exit()

from distutils.core import setup

os.chdir('../python')

setup(name = "lammps",
    version = "15Aug12",
    author = "Steve Plimpton",
    author_email = "sjplimp@sandia.gov",
    url = "http://lammps.sandia.gov",
    description = """LAMMPS molecular dynamics library""",
    py_modules = ["lammps"]
)
