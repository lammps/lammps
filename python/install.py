#!/usr/bin/env python

# copy LAMMPS src/liblammps.so and lammps.py to system dirs

instructions = """
Syntax: python install.py [-h] [pydir]
        pydir = target dir for lammps.py and liblammps.so
                default = Python site-packages dir
"""

import sys,os,commands

if (len(sys.argv) > 1 and sys.argv[1] == "-h") or len(sys.argv) > 2:
  print instructions
  sys.exit()

if len(sys.argv) == 2: pydir = sys.argv[1]
else: pydir = ""

# copy lammps.py to pydir if it exists
# if pydir not specified, install in site-packages via distutils setup()

if pydir:
  if not os.path.isdir(pydir):
    print "ERROR: pydir %s does not exist" % pydir
    sys.exit()
  str = "cp ../python/lammps.py %s" % pydir
  print str
  outstr = commands.getoutput(str)
  if len(outstr.strip()): print outstr
  str = "cp ../src/liblammps.so %s" % pydir
  print str
  outstr = commands.getoutput(str)
  sys.exit()
  
print "installing lammps.py in Python site-packages dir"

os.chdir('../python')                # in case invoked via make in src dir

from distutils.core import setup
sys.argv = ["setup.py","install"]    # as if had run "python setup.py install"
setup(name = "lammps",
      version = "15May15",
      author = "Steve Plimpton",
      author_email = "sjplimp@sandia.gov",
      url = "http://lammps.sandia.gov",
      description = "LAMMPS molecular dynamics library",
      py_modules = ["lammps"])
