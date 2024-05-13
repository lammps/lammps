#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the plumed2 library
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

version = "2.8.3"
mode = "static"

# help message

HELP = """
Syntax from src dir: make lib-plumed args="-b"
                 or: make lib-plumed args="-b -v 2.8.3"
                 or: make lib-plumed args="-p /usr/local/plumed2 -m shared"

Syntax from lib dir: python Install.py -b -v 2.8.3
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/plumed2 -m shared

Example:

make lib-plumed args="-b"   # download/build in lib/plumed/plumed2
make lib-plumed args="-p $HOME/plumed2 -m shared" # use existing Plumed2 installation in $HOME/plumed2
"""

# known checksums for different PLUMED versions. used to validate downloads.
checksums = { \
        '2.4.4' : '71ed465bdc7c2059e282dbda8d564e71', \
        '2.5.7' : '1ca36226fdb8110b1009aa61d615d4e5', \
        '2.6.6' : '6b470dcdce04c221ea42d8500b03c49b', \
        '2.7.6' : 'fb8c0ec10f97a9353eb123a5c4c35aa6', \
        '2.8.1' : '6bfe72ebdae63dc38a9ca27d9b0e08f8', \
        '2.8.2' : '599092b6a0aa6fff992612537ad98994', \
        '2.8.3' : '76d23cd394eba9e6530316ed1184e219', \
        '2.9.0' : '661eabeebee05cf84bbf9dc23d7d5f46', \
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
version = args.version

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
  filename = "plumed-src-%s.tar.gz" % version
  fallback = getfallback('plumed', url)
  print("Downloading plumed  ...")
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
        sys.exit("Checksum for plumed2 library does not match for fallback, too.")

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
  if (platform.system() == 'Darwin' and os.path.isfile("Makefile.lammps.%s.macosx" % mode)):
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
