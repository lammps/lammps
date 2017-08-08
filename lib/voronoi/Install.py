#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

try:
  import ssl
  try: from urllib.request import urlretrieve as geturl
  except: from urllib import urlretrieve as geturl
except:
  def geturl(url,fname):
    cmd = 'curl -L -o "%s" %s' % (fname,url)
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    return txt

# help message

help = """
Syntax from src dir: make lib-voronoi
                 or: make lib-voronoi args="-p /usr/local/voro++-0.4.6"
                 or: make lib-voronoi args="-v voro++-0.4.6 -b"
Syntax from lib dir: python Install.py -v voro++-0.4.6 -b
                 or: python Install.py
                 or: python Install.py -p /usr/local/voro++-0.4.6

specify one or more options, order does not matter

  -b = download and build the Voro++ library (default)
  -p = specify folder of existing Voro++ installation 
  -v = set version of Voro++ to download and build (default voro++-0.4.6)

Example:

make lib-voronoi args="-b"   # download/build in lib/voronoi/voro++-0.4.6
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

homepath = "."
homedir = version

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
    voropath = fullpath(args[iarg+1])
    pathflag = True
    buildflag = False
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,version)

if (pathflag):
    if not os.path.isdir(voropath): error("Voro++ path does not exist")
    homedir = voropath

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

# download and unpack Voro++ tarball

if buildflag:
  print("Downloading Voro++ ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,version))

  print("Unpacking Voro++ tarball ...")
  if os.path.exists("%s/%s" % (homepath,version)):
    cmd = 'rm -rf "%s/%s"' % (homepath,version)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s.tar.gz" % (homepath,version))
  if os.path.basename(homedir) != version:
    if os.path.exists(homedir):
      cmd = 'rm -rf "%s"' % homedir
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    os.rename("%s/%s" % (homepath,version),homedir)

# build Voro++

if buildflag:
  print("Building Voro++ ...")
  cmd = 'cd "%s"; make' % homedir
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))

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
