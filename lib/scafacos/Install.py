#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Scafacos library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,shutil
sys.path.append('..')
from install_helpers import error,fullpath,which,geturl

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

version = "scafacos-1.0.1"
url = "https://github.com/scafacos/scafacos/releases/download/v1.0.1/scafacos-1.0.1.tar.gz"
#url = "https://gigamove.rz.rwth-aachen.de/d/id/CTzyApN76MXMJ6/dd/100" % version

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
    if iarg+2 > nargs: error(help=help)
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error(help=help)
    scafacospath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error(help=help)

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
    shutil.rmtree("%s/%s" % (homepath,version))
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s.tar.gz" % (homepath,version))
  if os.path.basename(homedir) != version:
    if os.path.exists(homedir):
      shutil.rmtree(homedir)
    os.rename("%s/%s" % (homepath,version),homedir)

# build Scafacos

if buildflag:
  print("Building Scafacos ...")
  n_cpu = get_gpus()
  cmd = 'cd "%s"; ./configure --prefix="%s/build" --disable-doc --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m --with-internal-fftw --with-internal-pfft --with-internal-pnfft CC=mpicc FC=mpif90 CXX=mpicxx F77=; make -j%d; make install' % (homedir,homedir,n_cpu)
  try:
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

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
