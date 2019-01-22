#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,shutil
sys.path.append('..')
from install_helpers import error,get_cpus,fullpath,which,geturl

# help message

help = """
Syntax from src dir: make lib-voronoi args="-b"
                 or: make lib-voronoi args="-p /usr/local/voro++-0.4.6"
                 or: make lib-voronoi args="-b -v voro++-0.4.6"
Syntax from lib dir: python Install.py -b -v voro++-0.4.6
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/voro++-0.4.6

specify one or more options, order does not matter

  -b = download and build the Voro++ library
  -p = specify folder of existing Voro++ installation
  -v = set version of Voro++ to download and build (default voro++-0.4.6)

Example:

make lib-voronoi args="-b"   # download/build in lib/voronoi/voro++-0.4.6
make lib-voronoi args="-p $HOME/voro++-0.4.6" # use existing Voro++ installation in $HOME/voro++-0.4.6
"""

# settings

version = "voro++-0.4.6"
url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version


# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error(help=help)

homepath = "."
homedir = version

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
    voropath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error(help=help)

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,version)

if (pathflag):
    if not os.path.isdir(voropath): error("Voro++ path does not exist")
    homedir = voropath

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if (not buildflag and not pathflag):
    error("Have to use either -b or -p flag")

# download and unpack Voro++ tarball

if buildflag:
  print("Downloading Voro++ ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,version))

  print("Unpacking Voro++ tarball ...")
  if os.path.exists("%s/%s" % (homepath,version)):
    shutil.rmtree("%s/%s" % (homepath,version))
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s.tar.gz" % (homepath,version))
  if os.path.basename(homedir) != version:
    if os.path.exists(homedir):
      shutil.rmtree(homedir)
    os.rename("%s/%s" % (homepath,version),homedir)

# build Voro++

if buildflag:
  print("Building Voro++ ...")
  cmd = 'cd "%s"; make CXX=g++ CFLAGS="-fPIC -O3"' % homedir
  try:
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

# create 2 links in lib/voronoi to Voro++ src dir

if linkflag:
  print("Creating links to Voro++ include and lib files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = 'ln -s "%s/src" includelink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/src" liblink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
