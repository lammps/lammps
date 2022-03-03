#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the Voro++ library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil, tarfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = "voro++-0.4.6"
url = "https://download.lammps.org/thirdparty/%s.tar.gz" % version

# known checksums for different Voro++ versions. used to validate the download.
checksums = { \
        'voro++-0.4.6' : '2338b824c3b7b25590e18e8df5d68af9' \
        }

# extra help message

HELP = """
Syntax from src dir: make lib-voronoi args="-b"
                 or: make lib-voronoi args="-p /usr/local/voro++-0.4.6"
                 or: make lib-voronoi args="-b -v voro++-0.4.6"

Syntax from lib dir: python Install.py -b -v voro++-0.4.6
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/voro++-0.4.6

Example:

make lib-voronoi args="-b"   # download/build in lib/voronoi/voro++-0.4.6
make lib-voronoi args="-p $HOME/voro++-0.4.6" # use existing Voro++ installation in $HOME/voro++-0.4.6
"""

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the Voro++ library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing Voro++ installation")
parser.add_argument("-v", "--version", default=version,
                    help="set version of Voro++ to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
voropath = args.path

homepath = fullpath(".")
homedir = os.path.join(homepath, version)

if pathflag:
    if not os.path.isdir(voropath):
      sys.exit("Voro++ path %s does not exist" % voropath)
    homedir = fullpath(voropath)

# download and unpack Voro++ tarball

if buildflag:
  print("Downloading Voro++ ...")
  vorotar = os.path.join(homepath, version) + '.tar.gz'
  geturl(url, vorotar)

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], vorotar):
      sys.exit("Checksum for Voro++ library does not match")

  print("Unpacking Voro++ tarball ...")
  srcpath = os.path.join(homepath, version)
  if os.path.exists(srcpath):
    shutil.rmtree(srcpath)
  if tarfile.is_tarfile(vorotar):
    tgz = tarfile.open(vorotar)
    tgz.extractall(path=homepath)
    os.remove(vorotar)
  else:
    sys.exit("File %s is not a supported archive" % vorotar)
  if os.path.basename(homedir) != version:
    if os.path.exists(homedir):
      shutil.rmtree(homedir)
    os.rename(srcpath, homedir)

# build Voro++

if buildflag:
  print("Building Voro++ ...")
  cmd = 'cd "%s"; make CXX=g++ CFLAGS="-fPIC -O3"' % homedir
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

# create 2 links in lib/voronoi to Voro++ src dir

print("Creating links to Voro++ include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
os.symlink(os.path.join(homedir, 'src'), 'includelink')
os.symlink(os.path.join(homedir, 'src'), 'liblink')
