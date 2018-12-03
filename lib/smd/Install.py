#!/usr/bin/env python

# Install.py tool to download, unpack, and point to the Eigen library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,glob,subprocess,shutil
sys.path.append('..')
from install_helpers import error,get_cpus,fullpath,which,geturl

# help message

help = """
Syntax from src dir: make lib-smd args="-b"
                 or: make lib-smd args="-p /usr/include/eigen3"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/include/eigen3"
                 or: python Install.py -v 3.3.4 -b

specify one or more options, order does not matter

  -b = download and unpack/configure the Eigen library
  -p = specify folder holding an existing installation of Eigen
  -v = set version of Eigen library to download and set up (default = 3.3.4)


Example:

make lib-smd args="-b"   # download/build in default lib/smd/eigen-eigen-*
make lib-smd args="-p /usr/include/eigen3" # use existing Eigen installation in /usr/include/eigen3
"""

# settings

version = '3.3.4'
tarball = "eigen.tar.gz"


# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error(help=help)

homepath = "."
homedir = "eigen3"

buildflag = False
pathflag = False
linkflag = True

iarg = 0
while iarg < nargs:
  if args[iarg] == "-v":
    if iarg+2 > nargs: error(help=help)
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error(help=help)
    eigenpath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error(help=help)

homepath = fullpath(homepath)

if (pathflag):
  if not os.path.isdir(eigenpath): error("Eigen path does not exist")

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if (not buildflag and not pathflag):
    error("Have to use either -b or -p flag")

# download and unpack Eigen tarball
# use glob to find name of dir it unpacks to

if buildflag:
  print("Downloading Eigen ...")
  url = "http://bitbucket.org/eigen/eigen/get/%s.tar.gz" % version
  geturl(url,"%s/%s" % (homepath,tarball))

  print("Unpacking Eigen tarball ...")
  edir = glob.glob("%s/eigen-eigen-*" % homepath)
  for one in edir:
    if os.path.isdir(one):
      shutil.rmtree(one)
  cmd = 'cd "%s"; tar -xzvf %s' % (homepath,tarball)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  edir = glob.glob("%s/eigen-eigen-*" % homepath)
  os.rename(edir[0],"%s/%s" % (homepath,homedir))
  os.remove(tarball)

# create link in lib/smd to Eigen src dir

if linkflag:
  print("Creating link to Eigen files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if pathflag: linkdir = eigenpath
  else: linkdir = "%s/%s" % (homepath,homedir)
  cmd = 'ln -s "%s" includelink' % linkdir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
