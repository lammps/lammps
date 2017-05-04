#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the MS-CG library
# used to automate the steps described in the README file in this dir

import sys,os,re,commands

# help message

help = """
Syntax: python Install.py -h hpath hdir -g -b [suffix] -l
  specify one or more options, order does not matter
  -h = set home dir of MS-CG to be hpath/hdir
       hpath can be full path, contain '~' or '.' chars
       default hpath = . = lib/mscg
       default hdir = MSCG-release-master = what GitHub zipfile unpacks to
  -g = grab (download) zipfile from MS-CG GitHub website
       unpack it to hpath/hdir
       hpath must already exist
       if hdir already exists, it will be deleted before unpack
  -b = build MS-CG library in its src dir
       optional suffix specifies which src/Make/Makefile.suffix to use
       default suffix = g++_simple
  -l = create 2 softlinks (includelink,liblink) in lib/mscg to MS-CG src dir
"""

# settings

url = "https://github.com/uchicago-voth/MSCG-release/archive/master.zip"
zipfile = "MS-CG-master.zip"
zipdir = "MSCG-release-master"

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
homedir = zipdir

grabflag = 0
buildflag = 0
msuffix = "g++_simple"
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
  elif args[iarg] == "-b":
    buildflag = 1
    if iarg+1 < nargs and args[iarg+1][0] != '-':
      msuffix = args[iarg+1]
      iarg += 1
    iarg += 1
  elif args[iarg] == "-l":
    linkflag = 1
    iarg += 1
  else: error()

homepath = fullpath(homepath)
if not os.path.isdir(homepath): error("MS-CG path does not exist")
homedir = "%s/%s" % (homepath,homedir)

# download and unpack MS-CG zipfile

if grabflag:
  print "Downloading MS-CG ..."
  cmd = "curl -L %s > %s/%s" % (url,homepath,zipfile)
  print cmd
  print commands.getoutput(cmd)

  print "Unpacking MS-CG zipfile ..."
  if os.path.exists("%s/%s" % (homepath,zipdir)):
    commands.getoutput("rm -rf %s/%s" % (homepath,zipdir))
  cmd = "cd %s; unzip %s" % (homepath,zipfile)
  commands.getoutput(cmd)
  if os.path.basename(homedir) != zipdir:
    if os.path.exists(homedir): commands.getoutput("rm -rf %s" % homedir)
    os.rename("%s/%s" % (homepath,zipdir),homedir)

# build MS-CG

if buildflag:
  print "Building MS-CG ..."
  cmd = "cd %s/src; cp Make/Makefile.%s .; make -f Makefile.%s" % \
      (homedir,msuffix,msuffix)
  txt = commands.getoutput(cmd)
  print txt

# create 2 links in lib/mscg to MS-CG src dir

if linkflag:
  print "Creating links to MS-CG include and lib files"
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = "ln -s %s/src includelink" % homedir
  commands.getoutput(cmd)
  cmd = "ln -s %s/src liblink" % homedir
  commands.getoutput(cmd)
