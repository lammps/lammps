#!/usr/bin/env python

"""
Script to build a "binary wheel" for the 'pip' Python package manager for
the LAMMPS python module which includes the shared library file. After a
successful build the script attempts to install the wheel into a system
specific site-packages folder or - failing that - into the corresponding
user site-packages folder.  Called from the 'install-python' build target
in the GNU make and CMake based build systems.  Can also be called
independently and used to build the wheel without installing it.
"""

# copy LAMMPS shared library and lammps package to system dirs

from __future__ import print_function
import sys,os,shutil,time,glob,subprocess
from argparse import ArgumentParser

parser = ArgumentParser(prog='install.py',
                        description='LAMMPS python package installer script')

parser.add_argument("-p", "--package", required=True,
                    help="path to the LAMMPS Python package")
parser.add_argument("-l", "--lib", required=True,
                    help="path to the compiled LAMMPS shared library")
parser.add_argument("-n", "--noinstall", action="store_true", default=False,
                    help="only build a binary wheel. Don't attempt to install it")

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
try:
  txt = subprocess.check_output([sys.executable, '-m', 'virtualenv', 'buildwheel', '-p', sys.executable], stderr=subprocess.STDOUT, shell=False)
  print(txt.decode('UTF-8'))
except subprocess.CalledProcessError as err:
  sys.exit("Failed to create a virtualenv: {0}".format(err.output.decode('UTF-8')))

os.system(sys.executable + ' makewheel.py')

# remove temporary folders and files
shutil.rmtree('buildwheel',True)
shutil.rmtree('build',True)
shutil.rmtree('lammps.egg-info',True)
os.remove(os.path.join('lammps',os.path.basename(args.lib)))

if args.noinstall:
    exit(0)

print("Installing wheel")
for wheel in glob.glob('lammps-*.whl'):
  try:
    txt = subprocess.check_output([sys.executable, '-m', 'pip', 'install', '--force-reinstall', wheel], stderr=subprocess.STDOUT, shell=False)
    print(txt.decode('UTF-8'))
    continue
  except:
    pass
  try:
    print('Installing wheel into standard site-packages folder failed. Trying user folder now')
    txt = subprocess.check_output([sys.executable, '-m', 'pip', 'install', '--user', '--force-reinstall', wheel], stderr=subprocess.STDOUT, shell=False)
    print(txt.decode('UTF-8'))
  except:
    sys.exit('Failed to install wheel ' + wheel)
  shutil.copy(wheel, olddir)
  os.remove(wheel)
