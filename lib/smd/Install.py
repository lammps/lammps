#!/usr/bin/env python

# Install.py tool to download, unpack, and point to the Eigen library
# used to automate the steps described in the README file in this dir

import sys,os,re,glob,commands

# help message

help = """
Syntax from src dir: make lib-smd args="-h hpath hdir -g -l"
Syntax from lib dir: python Install.py -h hpath hdir -g -l

specify one or more options, order does not matter

  -h = set home dir of Eigen to be hpath/hdir
       hpath can be full path, contain '~' or '.' chars
       default hpath = . = lib/smd
       default hdir = "ee" = what tarball unpacks to (eigen-eigen-*)
  -g = grab (download) tarball from http://eigen.tuxfamily.org website
       unpack it to hpath/hdir
       hpath must already exist
       if hdir already exists, it will be deleted before unpack
  -l = create softlink (includelink) in lib/smd to Eigen src dir

Example:

make lib-smd args="-g -l"   # download/build in default lib/smd/eigen-eigen-*
"""

# settings

url = "http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz"
tarball = "eigen.tar.gz"

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

homepath = "."
homedir = "ee"

grabflag = 0
linkflag = 0

iarg = 0
while iarg < nargs:
  if args[iarg] == "-h":
    if iarg+3 > nargs: error()
    homepath = args[iarg+1]
    homedir = args[iarg+2]
    iarg += 3
  elif args[iarg] == "-g":
    grabflag = 1
    iarg += 1
  elif args[iarg] == "-l":
    linkflag = 1
    iarg += 1
  else: error()

homepath = fullpath(homepath)
if not os.path.isdir(homepath): error("Eigen path does not exist")

# download and unpack Eigen tarball
# glob to find name of dir it unpacks to

if grabflag:
  print "Downloading Eigen ..."
  cmd = "curl -L %s > %s/%s" % (url,homepath,tarball)
  print cmd
  print commands.getoutput(cmd)

  print "Unpacking Eigen tarball ..."
  edir = glob.glob("%s/eigen-eigen-*" % homepath)
  for one in edir:
    if os.path.isdir(one): commands.getoutput("rm -rf %s" % one)
  cmd = "cd %s; tar zxvf %s" % (homepath,tarball)
  commands.getoutput(cmd)
  if homedir != "ee":
    if os.path.exists(homedir): commands.getoutput("rm -rf %s" % homedir)
    edir = glob.glob("%s/eigen-eigen-*" % homepath)
    os.rename(edir[0],"%s/%s" % (homepath,homedir))

# create link in lib/smd to Eigen src dir

if linkflag:
  print "Creating link to Eigen files"
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if homedir == "ee":
    edir = glob.glob("%s/eigen-eigen-*" % homepath)
    linkdir = edir[0]
  else: linkdir = "%s/%s" % (homepath,homedir)
  cmd = "ln -s %s includelink" % linkdir
  commands.getoutput(cmd)
