#!/usr/bin/env python

# Install.py tool to download, unpack, and point to the Eigen library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,glob,subprocess

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

# print error message or help

def error(str=None):
  if not str: print(help)
  else: print("ERROR",str)
  sys.exit()

# expand to full path name
# process leading '~' or relative path

def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))

def which(program):
  def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath:
    if is_exe(program):
      return program
  else:
    for path in os.environ["PATH"].split(os.pathsep):
      path = path.strip('"')
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file

  return None

def geturl(url,fname):
  success = False

  if which('curl') != None:
    cmd = 'curl -L -o "%s" %s' % (fname,url)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling curl failed with: %s" % e.output.decode('UTF-8'))

  if not success and which('wget') != None:
    cmd = 'wget -O "%s" %s' % (fname,url)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling wget failed with: %s" % e.output.decode('UTF-8'))

  if not success:
    error("Failed to download source code with 'curl' or 'wget'")
  return

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

homepath = "."
homedir = "eigen3"

buildflag = False
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
