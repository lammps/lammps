#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the LATTE library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil, tarfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum, getfallback

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = '1.2.2'
suffix = 'gfortran'

# known checksums for different LATTE versions. used to validate the download.
checksums = { \
        '1.1.0' : '533635721ee222d0ed2925a18fb5b294', \
        '1.2.0' : '68bf0db879da5e068a71281020239ae7', \
        '1.2.1' : '85ac414fdada2d04619c8f936344df14', \
        '1.2.2' : '820e73a457ced178c08c71389a385de7', \
        }

# help message

HELP = """
Syntax from src dir: make lib-latte args="-b"
                 or: make lib-latte args="-p /usr/local/latte"
                 or: make lib-latte args="-m gfortran"
                 or: make lib-latte args="-b -v 1.2.2"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/latte
                 or: python Install.py -m gfortran
                 or: python Install.py -v 1.2.2 -b

Example:

make lib-latte args="-b -m gfortran"   # download/build in lib/latte
make lib-latte args="-p $HOME/latte"   # use existing LATTE installation
"""


pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the LATTE library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing LATTE installation")
parser.add_argument("-m", "--machine", choices=['gfortran', 'ifort', 'linalg', 'serial', 'mpi'],
                    help="suffix of a Makefile.lammps.* file used for linking LAMMPS with this library")
parser.add_argument("-v", "--version", default=version,
                    help="set version of LATTE to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

homepath = fullpath(".")

buildflag = args.build
pathflag = args.path is not None
version = args.version
suffixflag = args.machine is not None
if suffixflag:
  suffix = args.machine

if pathflag:
  lattedir = args.path
  if not os.path.isdir(lattedir):
    sys.exit("LATTE path %s does not exist" % lattedir)
  lattedir = fullpath(lattedir)

homedir = "LATTE-%s" % version

if buildflag:
  url = "https://github.com/lanl/LATTE/archive/v%s.tar.gz" % version
  lattepath = fullpath(homepath)
  lattedir = os.path.join(lattepath, homedir)
  fallback = getfallback('latte', url)
  filename = 'LATTE.tar.gz'

# download and unpack LATTE tarball

if buildflag:
  print("Downloading LATTE ...")
  try:
    geturl(url, filename)
  except:
    geturl(fallback, filename)

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], filename):
      print("Checksum did not match. Trying fallback URL", fallback)
      geturl(fallback, filename)
      if not checkmd5sum(checksums[version], filename):
        sys.exit("Checksum for LATTE library does not match for fallback, too.")

  print("Unpacking LATTE ...")
  if os.path.exists(lattedir):
    shutil.rmtree(lattedir)
  if tarfile.is_tarfile('LATTE.tar.gz'):
    tgz = tarfile.open('LATTE.tar.gz')
    tgz.extractall()
    os.remove('LATTE.tar.gz')
  else:
    sys.exit("File LATTE.tar.gz is not a supported archive")

  # build LATTE
  print("Building LATTE ...")
  cmd = 'cd "%s"; make' % lattedir
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    sys.exit("Make failed with:\n %s" % e.output.decode('UTF-8'))

# create 3 links in lib/latte to LATTE dirs

print("Creating links to LATTE files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
if os.path.isfile("filelink.o") or os.path.islink("filelink.o"):
  os.remove("filelink.o")
os.symlink(os.path.join(lattedir, 'src'), 'includelink')
os.symlink(lattedir, 'liblink')
os.symlink(os.path.join(lattedir, 'src', 'latte_c_bind.o'), 'filelink.o')

# copy Makefile.lammps.suffix to Makefile.lammps

if suffixflag or not os.path.exists("Makefile.lammps"):
  print("Creating Makefile.lammps")
  if os.path.exists("Makefile.lammps.%s" % suffix):
    shutil.copyfile("Makefile.lammps.%s" % suffix, 'Makefile.lammps')
