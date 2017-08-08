#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the MS-CG library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-mscg args="-p [path] -m [suffix]"
                 or: make lib-mscg args="-b -m [suffix]"
Syntax from lib dir: python Install.py -p [path]  -m [suffix]
Syntax from lib dir: python Install.py -b -m [suffix]

specify one or more options, order does not matter

  -b = download and build MS-CG library
  -p = specify folder of existing MS-CG installation
  -m = machine suffix specifies which src/Make/Makefile.suffix to use
       default suffix = g++_simple

Example:

make lib-mscg args="-b -m serial " # download/build in lib/mscg/MSCG-release-master with settings compatible with "make serial"
make lib-mscg args="-b -m mpi " # download/build in lib/mscg/MSCG-release-master with settings compatible with "make mpi"
make lib-mscg args="-p /usr/local/mscg-release " # use existing MS-CG installation in /usr/local/mscg-release
"""

# settings

url = "http://github.com/uchicago-voth/MSCG-release/archive/master.tar.gz"
tarfile = "MS-CG-master.tar.gz"
tardir = "MSCG-release-master"

# print error message or help

def error(str=None):
  if not str: print(help)
  else: print("ERROR",str)
  sys.exit()

# expand to full path name
# process leading '~' or relative path

def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))

def geturl(url,fname):
  cmd = 'curl -L -o "%s" %s' % (fname,url)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  return txt

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

homepath = "."
homedir = tardir

buildflag = False
pathflag = False
linkflag = True
msuffix = "g++_simple"

iarg = 0
while iarg < nargs:
  if args[iarg] == "-p":
    if iarg+2 > nargs: error()
    mscgpath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-m":
    if iarg+2 > nargs: error()
    msuffix = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,homedir)

if (pathflag):
    if not os.path.isdir(mscgpath): error("MS-CG path does not exist")
    homedir = mscgpath

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if (not buildflag and not pathflag):
    error("Have to use either -b or -p flag")

# download and unpack MS-CG tarfile

if buildflag:
  print("Downloading MS-CG ...")
  geturl(url,"%s/%s" % (homepath,tarfile))

  print("Unpacking MS-CG tarfile ...")
  if os.path.exists("%s/%s" % (homepath,tardir)):
    cmd = 'rm -rf "%s/%s"' % (homepath,tardir)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar -xzvf %s' % (homepath,tarfile)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s" % (homepath,tarfile))
  if os.path.basename(homedir) != tardir:
    if os.path.exists(homedir):
      cmd = 'rm -rf "%s"' % homedir
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    os.rename("%s/%s" % (homepath,tardir),homedir)

# build MS-CG

if buildflag:
  print("Building MS-CG ...")
  if os.path.exists("%s/src/Make/Makefile.%s" % (homedir,msuffix)):
    cmd = 'cd "%s/src"; cp Make/Makefile.%s .; make -f Makefile.%s' % \
        (homedir,msuffix,msuffix)
  elif os.path.exists("Makefile.%s" % msuffix):
    cmd = 'cd "%s/src"; cp ../../Makefile.%s .; make -f Makefile.%s' % \
        (homedir,msuffix,msuffix)
  else:
    error("Cannot find Makefile.%s" % msuffix)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))
  if not os.path.exists("Makefile.lammps"):
    print("Creating Makefile.lammps")
    if os.path.exists("Makefile.lammps.%s" % msuffix):
      cmd = 'cp Makefile.lammps.%s Makefile.lammps' % msuffix
    else:
      cmd = 'cp Makefile.lammps.default Makefile.lammps'
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  else: print("Makefile.lammps exists. Please check its settings")

# create 2 links in lib/mscg to MS-CG src dir

if linkflag:
  print("Creating links to MS-CG include and lib files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = 'ln -s "%s/src" includelink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/src" liblink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
