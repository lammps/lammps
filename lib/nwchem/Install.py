#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the NWChem PWDFT
library used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil, zipfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = "PWDFT"
url = "https://github.com/ebylaska/PWDFT/archive/master.zip"

# known checksums for different PWDFT versions. used to validate the download.
checksums = { \
        'PWDFT' : '???' \
        }

# extra help message

HELP = """
Syntax from src dir: make lib-nwchem args="-b"
                 or: make lib-nwchem args="-p /usr/local/PWDFT"
                 or: make lib-nwchem args="-b -v PWDFT"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/PWDFT

Example:

make lib-nwchem args="-b"   # download/build in lib/nwchem/PWDFT
make lib-nwchem args="-p $HOME/PWDFT" # use existing PWDFT installation in $HOME/PWDFT
"""

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the PWDFT library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing PWDFT installation")
parser.add_argument("-v", "--version", default=version,
                    help="set version of PWDFT to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
PWDFTpath = args.path

homepath = fullpath(".")
homedir = os.path.join(homepath, version)

if pathflag:
    if not os.path.isdir(PWDFTpath):
      sys.exit("PWDFT path %s does not exist" % PWDFTpath)
    homedir = fullpath(PWDFTpath)

# download and unpack PWDFT zipfile from GitHub

if buildflag:
  print("Downloading NWChem PWDFT ...")
  PWDFTzip = os.path.join(homepath, version) + '.zip'
  geturl(url, PWDFTzip)

  # verify downloaded archive integrity via md5 checksum, if known.
  #if version in checksums:
  #  if not checkmd5sum(checksums[version], vorotar):
  #    sys.exit("Checksum for Voro++ library does not match")

  print("Unpacking PWDFT zipfile ...")
  srcpath = os.path.join(homepath, version)
  if os.path.exists(srcpath):
    shutil.rmtree(srcpath)
  if zipfile.is_zipfile(PWDFTzip):
    zip = zipfile.ZipFile(PWDFTzip)
    zip.extractall(path=homepath)
    os.rename("PWDFT-master","PWDFT")
    os.remove(PWDFTzip)
  else:
    sys.exit("File %s is not a supported archive" % vorotar)

  #if os.path.basename(homedir) != version:
  #  if os.path.exists(homedir):
  #    shutil.rmtree(homedir)
  #  os.rename(srcpath, homedir)

# build PWDFT

if buildflag:
  print("Building PWDFT ...")
  cmd = 'cd %s; mkdir build_library; cd build_library; cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON; make -j' % homedir
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

# create 2 links in lib/voronoi to Voro++ src dir

print("Creating links to PWDFT lib file")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
os.symlink(os.path.join(homedir, 'build_library'), 'liblink')
