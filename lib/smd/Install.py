#!/usr/bin/env python

# Install.py tool to download, unpack, and point to the Eigen library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,glob,subprocess
try: from urllib.request import urlretrieve as geturl
except: from urllib import urlretrieve as geturl

# help message

help = """
Syntax from src dir: make lib-smd
                 or: make lib-smd args="-p /usr/include/eigen3"

Syntax from lib dir: python Install.py
                 or: python Install.py -p /usr/include/eigen3"
                 or: python Install.py -v 3.3.4 -b

specify one or more options, order does not matter

  -b = download and unpack/configure the Eigen library (default)
  -p = specify folder holding an existing installation of Eigen
  -v = set version of Eigen library to download and set up (default = 3.3.4)


Example:

make lib-smd args="-b"   # download/build in default lib/smd/eigen-eigen-*
"""

# settings

version = '3.3.4'
tarball = "eigen.tar.gz"

# print error message or help

def error(str=None):
  if not str: print(help)
  else: print("ERROR",str)
  sys.exit()

# expand to full path name
# process leading '~' or relative path

def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))

# parse args

args = sys.argv[1:]
nargs = len(args)

homepath = "."
homedir = "eigen3"

grabflag = True
buildflag = True
pathflag = False
linkflag = True

iarg = 0
while iarg < nargs:
  if args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error()
    eigenpath = fullpath(args[iarg+1])
    pathflag = True
    buildflag = False
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)

if (pathflag):
  if not os.path.isdir(eigenpath): error("Eigen path does not exist")

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

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
      subprocess.check_output('rm -rf "%s"' % one,stderr=subprocess.STDOUT,shell=True)
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
