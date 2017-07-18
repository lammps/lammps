#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir
from __future__ import print_function
import sys,os,re,subprocess
try: from urllib.request import urlretrieve as geturl
except: from urllib import urlretrieve as geturl

# help message

help = """
Syntax from src dir: make lib-voronoi args="-v version -h hpath hdir -g -b -l"
Syntax from lib dir: python Install.py -v version -h hpath hdir -g -b -l

specify one or more options, order does not matter

  -v = version of Voro++ to download and build
       default version = voro++-0.4.6 (current as of Jan 2015)
  -h = set home dir of Voro++ to be hpath/hdir
       hpath can be full path, contain '~' or '.' chars
       default hpath = . = lib/voronoi
       default hdir = voro++-0.4.6 = what tarball unpacks to
  -g = grab (download) tarball from math.lbl.gov/voro++ website
       unpack it to hpath/hdir
       hpath must already exist
       if hdir already exists, it will be deleted before unpack
  -b = build Voro++ library in its src dir
  -l = create 2 softlinks (includelink,liblink) in lib/voronoi to Voro++ src dir

Example:

make lib-voronoi args="-g -b -l"   # download/build in lib/voronoi/voro++-0.4.6
"""

# settings

version = "voro++-0.4.6"
url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version

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
if nargs == 0: error()

homepath = "."
homedir = version

grabflag = False
buildflag = False
linkflag = False

iarg = 0
while iarg < nargs:
  if args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-h":
    if iarg+3 > nargs: error()
    homepath = args[iarg+1]
    homedir = args[iarg+2]
    iarg += 3
  elif args[iarg] == "-g":
    grabflag = True
    iarg += 1
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  elif args[iarg] == "-l":
    linkflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)
if not os.path.isdir(homepath): error("Voro++ path does not exist")
homedir = "%s/%s" % (homepath,homedir)

# download and unpack Voro++ tarball

if grabflag:
  print("Downloading Voro++ ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,version))
  
  print("Unpacking Voro++ tarball ...")
  if os.path.exists("%s/%s" % (homepath,version)):
    cmd = 'rm -rf "%s/%s"' % (homepath,version)
    subprocess.check_output(cmd,shell=True)
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,shell=True)
  if os.path.basename(homedir) != version:
    if os.path.exists(homedir):
      cmd = 'rm -rf "%s"' % homedir
      subprocess.check_output(cmd,shell=True)
    os.rename("%s/%s" % (homepath,version),homedir)

# build Voro++

if buildflag:
  print("Building Voro++ ...")
  cmd = 'cd "%s"; make' % homedir
  txt = subprocess.check_output(cmd,shell=True)
  print(txt)

# create 2 links in lib/voronoi to Voro++ src dir

if linkflag:
  print("Creating links to Voro++ include and lib files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = ['ln -s "%s/src" includelink' % homedir, 'includelink']
  subprocess.check_output(cmd,shell=True)
  cmd = ['ln -s "%s/src" liblink' % homedir]
  subprocess.check_output(cmd,shell=True)
