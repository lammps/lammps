#!/usr/bin/env python

"""
Install.py tool to download, unpack, and point to the Eigen library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, glob, shutil, tarfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = '3.4.0'
tarball = "eigen.tar.gz"

# known checksums for different Eigen versions. used to validate the download.
checksums = { \
              '3.4.0' : '4c527a9171d71a72a9d4186e65bea559', \
              '3.3.9' : '609286804b0f79be622ccf7f9ff2b660', \
              '3.3.7' : '9e30f67e8531477de4117506fe44669b' \
}


# help message

HELP = """
Syntax from src dir: make lib-machdyn args="-b"
                 or: make lib-machdyn args="-p /usr/include/eigen3"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/include/eigen3"
                 or: python Install.py -v 3.4.0 -b

Example:

make lib-machdyn args="-b"   # download/build in default lib/machdyn/eigen-eigen-*
make lib-machdyn args="-p /usr/include/eigen3" # use existing Eigen installation in /usr/include/eigen3
"""

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the Eigen3 library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing Eigen installation")
parser.add_argument("-v", "--version", default=version,
                    help="set version of Eigen to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

homepath = fullpath(".")
eigenpath = os.path.join(homepath, "eigen3")

buildflag = args.build
pathflag = args.path is not None
version = args.version

if pathflag:
  eigenpath = args.path
  if not os.path.isdir(eigenpath):
    sys.exit("Eigen path %s does not exist" % eigenpath)
  eigenpath = fullpath(eigenpath)

# download and unpack Eigen tarball
# use glob to find name of dir it unpacks to

if buildflag:
  print("Downloading Eigen ...")
  eigentar = os.path.join(homepath, tarball)
  url = "https://download.lammps.org/thirdparty/eigen-%s.tar.gz" %  version
  geturl(url, eigentar)

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    print("checking version %s\n" % version)
    if not checkmd5sum(checksums[version], eigentar):
      sys.exit("Checksum for Eigen library does not match")


  print("Cleaning up old folders ...")
  edir = glob.glob(os.path.join(homepath, "eigen-*"))
  edir.append(eigenpath)
  for one in edir:
    if os.path.isdir(one):
      shutil.rmtree(one)

  print("Unpacking Eigen tarball ...")
  if tarfile.is_tarfile(eigentar):
    tgz = tarfile.open(eigentar)
    tgz.extractall(path=homepath)
    os.remove(eigentar)
  else:
    sys.exit("File %s is not a supported archive" % eigentar)
  edir = os.path.join(homepath, "eigen-%s" % version)
  os.rename(edir, eigenpath)

# create link in lib/machdyn to Eigen src dir

print("Creating link to Eigen include folder")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
linkdir = eigenpath
os.symlink(linkdir, 'includelink')
