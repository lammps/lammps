#!/usr/bin/env python

"""
Installer script to install the LAMMPS python package and the corresponding
shared library into either the system-wide site-packages tree, or - failing
that - into the corresponding user tree. Called from the 'install-python'
build target in the conventional and CMake based build systems
"""

# copy LAMMPS shared library and lammps package to system dirs

from __future__ import print_function
import sys,os,shutil,time,glob
from argparse import ArgumentParser

parser = ArgumentParser(prog='install.py',
                        description='LAMMPS python package installer script')

parser.add_argument("-p", "--package", required=True,
                    help="path to the LAMMPS Python package")
parser.add_argument("-l", "--lib", required=True,
                    help="path to the compiled LAMMPS shared library")

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

# we need to switch to the folder of the python package
olddir = os.path.abspath('.')
os.chdir(os.path.dirname(args.package))

print("Purging existing wheels...")
for wheel in glob.glob('lammps-*.whl'):
    print("deleting " + wheel)
    os.remove(wheel)

# create virtual environment for building the wheel
shutil.rmtree('buildwheel',True)
os.putenv('LAMMPS_SHARED_LIB',args.lib)
#os.environ['LAMMPS_SHARED_LIB'] = args.lib
shutil.copy(args.lib,'lammps')
os.system(sys.executable + ' -m virtualenv buildwheel')
os.system(sys.executable + ' makewheel.py')

# remove temporary folders and files
shutil.rmtree('buildwheel',True)
shutil.rmtree('build',True)
shutil.rmtree('lammps.egg-info',True)
os.remove(os.path.join('lammps',os.path.basename(args.lib)))

print("Installing wheel")
for wheel in glob.glob('lammps-*.whl'):
  os.system(sys.executable + " -m pip install --force-reinstall " + wheel)
  shutil.copy(wheel, olddir)
  os.remove(wheel)
