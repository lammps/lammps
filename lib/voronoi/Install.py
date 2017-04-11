#!/usr/bin/env python

# install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

import sys,os,re,urllib,commands

# help message

help = """
Syntax: install.py -v version -g gdir [gname] -b bdir -l ldir
  specify one or more options, order does not matter
  gdir,bdir,ldir can be paths relative to lib/latte, full paths, or contain ~
  -v = version of Voro++ to download and build
       default = voro++-0.4.6 (current as of Jan 2015)
  -g = grab (download) from math.lbl.gov/voro++ website
       unpack tarfile in gdir to produce version dir (e.g. voro++-0.4.6)
       if optional gname specified, rename version dir to gname within gdir
  -b = build Voro++, bdir = Voro++ home directory
       note that bdir must include the version suffix unless renamed
  -l = create 2 softlinks (includelink,liblink)
         in lib/voronoi to src dir of ldir = Voro++ home directory
       note that ldir must include the version suffix unless renamed
"""

# settings

version = "voro++-0.4.6"
url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version

# print error message or help

def error(str=None):
  if not str: print help
  else: print "ERROR",str
  sys.exit()

# expand to full path name
# process leading '~' or relative path
  
def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))
  
# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

grabflag = 0
buildflag = 0
linkflag = 0

iarg = 0
while iarg < nargs:
  if args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2  
  elif args[iarg] == "-g":
    if iarg+2 > nargs: error()
    grabflag = 1
    grabdir = args[iarg+1]
    grabname = None
    if iarg+2 < nargs and args[iarg+2][0] != '-':
      grabname = args[iarg+2]
      iarg += 1
    iarg += 2
  elif args[iarg] == "-b":
    if iarg+2 > nargs: error()
    buildflag = 1
    builddir = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-l":
    if iarg+2 > nargs: error()
    linkflag = 1
    linkdir = args[iarg+1]
    iarg += 2
  else: error()

# download and unpack Voro++ tarball

if grabflag:
  print "Downloading Voro++ ..."
  grabdir = fullpath(grabdir)
  if not os.path.isdir(grabdir): error("Grab directory does not exist")
  urllib.urlretrieve(url,"%s/%s.tar.gz" % (grabdir,version))
  
  print "Unpacking Voro++ tarball ..."
  tardir = "%s/%s" % (grabdir,version)
  if os.path.exists(tardir): commands.getoutput("rm -rf %s" % tardir)
  cmd = "cd %s; tar zxvf %s.tar.gz" % (grabdir,version)
  txt = commands.getoutput(cmd)
  print tardir,grabdir,grabname
  if grabname: os.rename(tardir,"%s/%s" % (grabdir,grabname))

# build Voro++

if buildflag:
  print "Building Voro++ ..."
  cmd = "cd %s; make" % builddir
  txt = commands.getoutput(cmd)
  print txt

# create 2 links in lib/voronoi to Voro++ src dir

if linkflag:
  print "Creating links to Voro++ include and lib files"
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = "ln -s %s/src includelink" % linkdir
  commands.getoutput(cmd)
  cmd = "ln -s %s/src liblink" % linkdir
  commands.getoutput(cmd)
