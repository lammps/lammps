#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the LATTE library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,shutil
sys.path.append('..')
from install_helpers import error,get_cpus,fullpath,which,geturl,checkmd5sum

# help message

help = """
Syntax from src dir: make lib-latte args="-b"
                 or: make lib-latte args="-p /usr/local/latte"
                 or: make lib-latte args="-m gfortran"
                 or: make lib-latte args="-b -v 1.2.1"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/latte
                 or: python Install.py -m gfortran
                 or: python Install.py -v 1.2.1 -b

specify one or more options, order does not matter

  -b = download and build the LATTE library
  -p = specify folder of existing LATTE installation
  -m = copy Makefile.lammps.suffix to Makefile.lammps
  -v = set version of LATTE library to download and set up (default = 1.2.1)

Example:

make lib-latte args="-b -m gfortran"   # download/build in lib/latte
make lib-latte args="-p $HOME/latte"   # use existing LATTE installation
"""

# settings

version = '1.2.1'

# known checksums for different LATTE versions. used to validate the download.
checksums = { \
        '1.1.0' : '533635721ee222d0ed2925a18fb5b294', \
        '1.2.0' : '68bf0db879da5e068a71281020239ae7', \
        '1.2.1' : '85ac414fdada2d04619c8f936344df14', \
        }

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error(help=help)

homepath = "."

buildflag = False
pathflag = False
suffixflag = False
linkflag = True

iarg = 0
while iarg < nargs:
  if args[iarg] == "-p":
    if iarg+2 > nargs: error(help=help)
    lattedir = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  elif args[iarg] == "-m":
    if iarg+2 > nargs: error(help=help)
    suffix = args[iarg+1]
    suffixflag = True
    iarg += 2
  elif args[iarg] == "-v":
    if iarg+2 > nargs: error(help=help)
    version = args[iarg+1]
    iarg += 2
  else: error(help=help)

homedir = "LATTE-%s" % version

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if buildflag:
  url = "https://github.com/lanl/LATTE/archive/v%s.tar.gz" % version
  lattepath = fullpath(homepath)
  lattedir = "%s/%s" % (lattepath,homedir)

if pathflag:
  if not os.path.isdir(lattedir): error("LATTE path does not exist")

# download and unpack LATTE tarball

if buildflag:
  print("Downloading LATTE ...")
  geturl(url,"LATTE.tar.gz")

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version],'LATTE.tar.gz'):
      error("Checksum for LATTE library does not match")

  print("Unpacking LATTE ...")
  if os.path.exists(lattedir):
    shutil.rmtree(lattedir)
  cmd = 'cd "%s"; tar zxvf LATTE.tar.gz' % lattepath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/LATTE.tar.gz" % lattepath)

# build LATTE

if buildflag:
  print("Building LATTE ...")
  cmd = 'cd "%s"; make' % lattedir
  try:
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

# create 3 links in lib/latte to LATTE dirs
# do this -b or -p is set

if buildflag or pathflag:
  print("Creating links to LATTE files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  if os.path.isfile("filelink.o") or os.path.islink("filelink.o"):
    os.remove("filelink.o")
  cmd = 'ln -s "%s/src" includelink' % lattedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s" liblink' % lattedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/src/latte_c_bind.o" filelink.o' % lattedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

# copy Makefile.lammps.suffix to Makefile.lammps

if suffixflag:
  print("Creating Makefile.lammps")
  if os.path.exists("Makefile.lammps.%s" % suffix):
    cmd = 'cp Makefile.lammps.%s Makefile.lammps' % suffix
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
