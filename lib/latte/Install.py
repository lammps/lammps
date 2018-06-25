#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the LATTE library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-latte args="-b"
                 or: make lib-latte args="-p /usr/local/latte"
                 or: make lib-latte args="-m gfortran"
                 or: make lib-latte args="-b -v 1.2.1"

Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/latte
                 or: python Install.py -m gfortran
                 or: python Install.py -v 1.2.1 -b

specify one or more options, order does not matter

  -b = download and build the LATTE library
  -p = specify folder of existing LATTE installation
  -m = copy Makefile.lammps.suffix to Makefile.lammps
  -v = set version of LATTE library to download and set up (default = 1.1.1)

Example:

make lib-latte args="-b -m gfortran"   # download/build in lib/latte
make lib-latte args="-p $HOME/latte"   # use existing LATTE installation
"""

# settings

version = '1.2.1'

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
    print(cmd)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling curl failed with: %s" % e.output.decode('UTF-8'))

  if not success and which('wget') != None:
    cmd = 'wget -O "%s" %s' % (fname,url)
    print(cmd)
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

buildflag = False
pathflag = False
suffixflag = False
linkflag = True

iarg = 0
while iarg < nargs:
  if args[iarg] == "-p":
    if iarg+2 > nargs: error()
    lattedir = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  elif args[iarg] == "-m":
    if iarg+2 > nargs: error()
    suffix = args[iarg+1]
    suffixflag = True
    iarg += 2
  elif args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2
  else: error()

homedir = "LATTE-%s" % version

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if buildflag:
  url = "https://github.com/lanl/LATTE/archive/v%s.tar.gz" % version
  lattepath = fullpath(homepath)
  lattedir = "%s/%s" % (lattepath,homedir)

if pathflag:
  if not os.path.isdir(lattedir): error("LATTE path does not exist")

# download and unpack LATTE tarball

if buildflag:
  print("Downloading LATTE ...")
  geturl(url,"LATTE.tar.gz")

  print("Unpacking LATTE zipfile ...")
  if os.path.exists(lattedir):
    cmd = 'rm -rf "%s"' % lattedir
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar zxvf LATTE.tar.gz' % lattepath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/LATTE.tar.gz" % lattepath)

# build LATTE

if buildflag:
  print("Building LATTE ...")
  cmd = 'cd "%s"; make' % lattedir
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))

# create 3 links in lib/latte to LATTE dirs
# do this -b or -p is set
  
if buildflag or pathflag:
  print("Creating links to LATTE files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  if os.path.isfile("filelink.o") or os.path.islink("filelink.o"):
    os.remove("filelink.o")
  cmd = 'ln -s "%s/src" includelink' % lattedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s" liblink' % lattedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/src/latte_c_bind.o" filelink.o' % lattedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

# copy Makefile.lammps.suffix to Makefile.lammps

if suffixflag:
  print("Creating Makefile.lammps")
  if os.path.exists("Makefile.lammps.%s" % suffix):
    cmd = 'cp Makefile.lammps.%s Makefile.lammps' % suffix
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
