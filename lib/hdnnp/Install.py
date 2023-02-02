#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the n2p2 library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, platform, subprocess, shutil
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import get_cpus, fullpath, geturl, checkmd5sum, getfallback

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = "2.2.0"

# help message

HELP = """
Syntax from src dir: make lib-hdnnp args="-b"
                 or: make lib-hdnnp args="-b -v 2.1.4"
                 or: make lib-hdnnp args="-p /usr/local/n2p2"

Syntax from lib dir: python Install.py -b -v 2.1.4
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/n2p2

Example:

make lib-hdnnp args="-b"   # download/build in lib/hdnnp/n2p2
make lib-hdnnp args="-p $HOME/n2p2" # use existing n2p2 installation in $HOME/n2p2
"""

# known checksums for different n2p2 versions. used to validate the download.
checksums = { \
        '2.2.0' : 'a2d9ab7f676b3a74a324fc1eda0a911d', \
        '2.1.4' : '9595b066636cd6b90b0fef93398297a5', \
        }

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the n2p2 library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing n2p2 installation")
parser.add_argument("-v", "--version", default=version, choices=checksums.keys(),
                    help="set version of n2p2 to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
n2p2path = args.path

homepath = fullpath('.')
homedir = "%s/n2p2" % (homepath)

if pathflag:
    if not os.path.isdir(n2p2path):
      sys.exit("n2p2 path %s does not exist" % n2p2path)
    homedir = fullpath(n2p2path)
    if not os.path.isfile(os.path.join(homedir, 'include', 'InterfaceLammps.h')):
      sys.exit("No n2p2 installation found at %s" % n2p2path)

# download and unpack n2p2 tarball

if buildflag:
  url = "https://github.com/CompPhysVienna/n2p2/archive/v%s.tar.gz" % (version)
  filename = "n2p2-%s.tar.gz" % version
  fallback = getfallback('n2p2', url)
  print("Downloading n2p2 from", url)
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
        sys.exit("Checksum for n2p2 library does not match for fallback, too.")

  print("Unpacking n2p2 source tarball ...")
  if os.path.exists("%s/n2p2-%s" % (homepath, version)):
    shutil.rmtree("%s/n2p2-%s" % (homepath, version))
  if os.path.exists(homedir):
    shutil.rmtree(homedir)
  cmd = 'cd "%s"; tar -xzvf %s' % (homepath, filename)
  subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  os.remove(os.path.join(homepath, filename))

  # build n2p2
  print("Building n2p2 ...")
  n_cpus = get_cpus()
  cmd = 'unset MAKEFLAGS MAKELEVEL MAKEOVERRIDES MFLAGS && cd %s/n2p2-%s/src && make -j%d libnnpif' % (homepath, version, n_cpus)
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

  # set correct homedir for linking step
  homedir = "%s/n2p2-%s" % (homepath, version)

# create 2 links in lib/hdnnp to n2p2 installation dir

print("Creating links to n2p2 include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
if os.path.isfile("Makefile.lammps") or os.path.islink("Makefile.lammps"):
  os.remove("Makefile.lammps")
os.symlink(os.path.join(homedir, 'include'), 'includelink')
os.symlink(os.path.join(homedir, 'lib'), 'liblink')
os.symlink(os.path.join(homedir, 'lib', 'Makefile.lammps-extra'), 'Makefile.lammps')
