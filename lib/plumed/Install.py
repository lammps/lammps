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

version = "2.8.1"
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
        '2.5.3' : 'de30d6e7c2dcc0973298e24a6da24286', \
        '2.5.4' : 'f31b7d16a4be2e30aa7d5c19c3d37853', \
        '2.5.7' : '1ca36226fdb8110b1009aa61d615d4e5', \
        '2.6.0' : '204d2edae58d9b10ba3ad460cad64191', \
        '2.6.1' : '89a9a450fc6025299fe16af235957163', \
        '2.6.3' : 'a9f8028fd74528c2024781ea1fdefeee', \
        '2.6.5' : 'b67356f027e5c2747823b0422c3b0ec2', \
        '2.6.6' : '6b470dcdce04c221ea42d8500b03c49b', \
        '2.7.0' : '95f29dd0c067577f11972ff90dfc7d12', \
        '2.7.1' : '4eac6a462ec84dfe0cec96c82421b8e8', \
        '2.7.2' : 'cfa0b4dd90a81c25d3302e8d97bfeaea', \
        '2.7.3' : 'f00cc82edfefe6bb3df934911dbe32fb', \
        '2.7.4' : 'f858e0b6aed173748fc85b6bc8a9dcb3', \
        '2.7.5' : '2aca1986d6c1ca3ba7e9eb51b1102792', \
        '2.8.0' : '489b23daba70da78cf0506cbc31689c6', \
        '2.8.1' : '6bfe72ebdae63dc38a9ca27d9b0e08f8', \
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
