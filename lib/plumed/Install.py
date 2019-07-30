#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the plumed2 library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, platform, subprocess, shutil
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import get_cpus, fullpath, geturl, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = "2.5.2"
mode = "static"

# help message

HELP = """
Syntax from src dir: make lib-plumed args="-b"
                 or: make lib-plumed args="-b -v 2.4.3"
                 or: make lib-plumed args="-p /usr/local/plumed2 -m shared"

Syntax from lib dir: python Install.py -b -v 2.4.3
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/plumed2 -m shared

Example:

make lib-plumed args="-b"   # download/build in lib/plumed/plumed2
make lib-plumed args="-p $HOME/plumed2 -m shared" # use existing Plumed2 installation in $HOME/plumed2
"""

# known checksums for different PLUMED versions. used to validate the download.
checksums = { \
        '2.4.2' : '88188743a6e03ef076e5377d03ebb0e7', \
        '2.4.3' : 'b1be7c48971627febc11c61b70767fc5', \
        '2.4.4' : '71ed465bdc7c2059e282dbda8d564e71', \
        '2.5.0' : '6224cd089493661e19ceacccd35cf911', \
        '2.5.1' : 'c2a7b519e32197a120cdf47e0f194f81', \
        '2.5.2' : 'bd2f18346c788eb54e1e52f4f6acf41a', \
        }

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the plumed2 library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing plumed2 installation")
parser.add_argument("-v", "--version", default=version, choices=checksums.keys(),
                    help="set version of plumed to download and build (default: %s)" % version)
parser.add_argument("-m", "--mode", default=mode, choices=['static', 'shared', 'runtime'],
                    help="set plumed linkage mode: static (default), shared, or runtime")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
plumedpath = args.path
mode = args.mode

homepath = fullpath('.')
homedir = "%s/plumed2" % (homepath)

if pathflag:
    if not os.path.isdir(plumedpath):
      sys.exit("Plumed2 path %s does not exist" % plumedpath)
    homedir = fullpath(plumedpath)
    if not os.path.isdir(os.path.join(homedir, 'include', 'plumed', 'core')):
      sys.exit("No Plumed2 installation found at %s" % plumedpath)

# download and unpack plumed2 tarball

if buildflag:
  url = "https://github.com/plumed/plumed2/releases/download/v%s/plumed-src-%s.tgz" % (version, version)
  filename = "plumed-src-%s.tar.gz" %version
  print("Downloading plumed  ...")
  geturl(url, filename)

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], filename):
      sys.exit("Checksum for plumed2 library does not match")

  print("Unpacking plumed2 source tarball ...")
  if os.path.exists("%s/plumed-%s" % (homepath, version)):
    shutil.rmtree("%s/plumed-%s" % (homepath, version))
  if os.path.exists(homedir):
    shutil.rmtree(homedir)
  cmd = 'cd "%s"; tar -xzvf %s' % (homepath, filename)
  subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  os.remove(os.path.join(homepath, filename))

  # build plumed
  print("Building plumed ...")
  n_cpus = get_cpus()
  cmd = 'cd %s/plumed-%s; ./configure --prefix=%s --enable-modules=all --enable-static-patch ; make -j%d ; make install' % (homepath, version, homedir, n_cpus)
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

# create 2 links in lib/plumed to plumed2 installation dir

print("Creating links to plumed2 include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
os.symlink(os.path.join(homedir, 'include'), 'includelink')
libpath = os.path.join(homedir, 'lib64')
if not os.path.exists(libpath):
  libpath = os.path.join(homedir, 'lib')
os.symlink(libpath, 'liblink')
if os.path.isfile("Makefile.lammps.%s" % mode):
  print("Creating Makefile.lammps")
  plumedinc = os.path.join('liblink', 'plumed', 'src', 'lib', 'Plumed.inc.' + mode)
  lines1 = open(plumedinc, 'r').readlines()
  if (platform.system() == 'Darwin' and os.path.isfile("Makefile.lammps.%s.macosx")):
    lines2 = open("Makefile.lammps.%s.macosx" % mode, 'r').readlines()
  else:
    lines2 = open("Makefile.lammps.%s" % mode, 'r').readlines()
  fp = open("Makefile.lammps", 'w')
  fp.write("PLUMED_LIBDIR=" + os.path.join(homedir, "lib\n"))
  for line in lines1:
    fp.write(line)
  for line in lines2:
    fp.write(line)
  fp.close()
