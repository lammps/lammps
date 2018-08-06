#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Scafacos library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-scafacos args="-b"
                 or: make lib-scafacos args="-p /usr/local/scafacos"
Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/scafacos

specify zero or more options, order does not matter

  -b = download and build the Scafacos library
  -p = specify folder of existing Scafacos installation

   always creates includelink, liblink to Scafacos dirs

Example:

make lib-scafacos args="-b"   # download/build in lib/scafacos/scafacos
make lib-scafacos args="-p $HOME/scafacos" # use existing Scafacos installation in $HOME
"""

# settings

version = "scafacos-1.0.0"
url = "https://github.com/scafacos/scafacos/releases/download/v1.0.0/scafacos-1.0.0.tar.gz"
#url = "https://gigamove.rz.rwth-aachen.de/d/id/CTzyApN76MXMJ6/dd/100" % version

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
    print("Wget command: %s" % cmd)
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

homepath = "."

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
    scafacospath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,version)

if (pathflag):
    if not os.path.isdir(scafacospath): error("Scafacos path does not exist")
    homedir =scafacospath

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

# download and unpack Scafacos tarball

if buildflag:
  print("Downloading Scafacos ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,version))

  print("Unpacking Scafacos tarball ...")
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

# build Scafacos

if buildflag:
  print("Building Scafacos ...")
  cmd = 'cd "%s"; CC=mpicc FC=mpif90 CXX=mpicxx ./configure --prefix="`pwd`/build" --disable-doc --enable-fcs-solvers=fmm,p2nfft,direct,ewald > log.txt; make -j; make install' % homedir
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))

# create 2 links in lib/scafacos to Scafacos include/lib dirs

if linkflag:
  print("Creating links to Scafacos include and lib files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = 'ln -s "%s/build/include" includelink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/build/lib" liblink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
