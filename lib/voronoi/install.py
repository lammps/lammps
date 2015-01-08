#!usr/local/python

# install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

import sys,os,urllib,commands

help = """
Syntax: install.py -d dir -v version -g -b -i installdir -l incdir libdir
        specify one or more options, order does not matter
        -d = dir to download tarball to, unpack tarball in, perform build in
             dir will be created if it doesn't exist (only last level)
             default = this dir
        -v = version of Voro++ to download and work with
             default = voro++-0.4.6 (current as of Jan 2015)
        -g = download (grab) tarball from
             http://math.lbl.gov/voro++/download/dir/version
        -b = build Voro++ by invoking "make" in its home dir
             no default
        -i = install Voro++ by invoking "make install" in its home dir
             installdir arg is optional:
               if not specified, installs at PREFIX defined in config.mk file
               if specified, will overwrite PREFIX and install there
             if PREFIX starts with /usr, will invoke "sudo make install"
        -l = create two links to incdir and libdir
             incdir and libdir are optional (specify neither or both):
               if specified, includelink and liblink are to those two dirs
               if not specified and no install, links are to Voro++ src dir
               if not specified and install performed,
                 links are to include and lib dirs under PREFIX
"""

def error():
  print help
  sys.exit()
  
# parse args

args = sys.argv

if len(args) == 1: error()

dir = "."
version = "voro++-0.4.6"
grabflag = 0
buildflag = 0
installflag = 0
linkflag = 0

iarg = 1
while iarg < len(args):
  if args[iarg] == "-d":
    if iarg+2 > len(args): error()
    dir = args[iarg+1]
    iarg += 2  
  elif args[iarg] == "-v":
    if iarg+2 > len(args): error()
    version = args[iarg+1]
    iarg += 2  
  elif args[iarg] == "-g":
    grabflag = 1
    iarg += 1
  elif args[iarg] == "-b":
    buildflag = 1
    iarg += 1
  elif args[iarg] == "-i":
    if iarg+2 > len(args): error()
    if iarg+1 == len(args) or args[iarg+1][0] == '-':
      installdir = ""
      iarg += 1
    else:
      installdir = args[iarg+1]
      iarg += 2  
  elif args[iarg] == "-l":
    if iarg+2 > len(args): error()
    if iarg+1 == len(args) or args[iarg+1][0] == '-':
      incdir = libdir = ""
      iarg += 1
    else:
      incdir = args[iarg+1]
      libdir = args[iarg+2]
      iarg += 3
  else: error()

url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version

# create dir if does not exist

if not os.path.isdir(dir):
  if os.path.isfile(dir):
    print "ERROR: Dir already exists as file"
    sys.exit()
  os.mkdir(dir)
  if not os.path.isdir(dir):
    print "ERROR: Unable to create dir"
    sys.exit()

# download and unpack tarball

if grabflag:
  print "Downloading Voro++ tarball ..."
  urllib.urlretrieve(url,"%s/%s.tar.gz" % (dir,version))
  print "Unpacking Voro++ tarball ..."
  cmd = "cd %s; tar zxvf %s.tar.gz" % (dir,version)
  txt = commands.getoutput(cmd)

# build Voro++

if buildflag:
  print "Building Voro++ ..."
  cmd = "cd %s/%s; make" % (dir,version)
  txt = commands.getoutput(cmd)
  print txt

# install Voro++



  
# create links in this dir to Voro++ include and lib files

if linkflag:
  print "Creating links to Voro++ include and lib files"
  cmd = "ln -s %s/src includelink" % version
  txt = commands.getoutput(cmd)
  cmd = "ln -s %s/src liblink" % version
  txt = commands.getoutput(cmd)
