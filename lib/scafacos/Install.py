#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Scafacos library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess
from argparse import ArgumentParser

# help message

help_text = """
Syntax from src dir: make lib-scafacos args="-b"
                 or: make lib-scafacos args="-p /usr/local/scafacos"
Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/scafacos

specify zero or more options, order does not matter

  -b = download and build the Scafacos library
  -p = specify folder of existing Scafacos installation
  -v = set version of Scafacos to download and build (default scafacos-1.0.1)

   always creates includelink, liblink to Scafacos dirs

Example:

make lib-scafacos args="-b"   # download/build in lib/scafacos/scafacos
make lib-scafacos args="-p $HOME/scafacos" # use existing Scafacos installation in $HOME
"""

# settings

version = "scafacos-1.0.1"
url = "https://github.com/scafacos/scafacos/releases/download/v1.0.1/scafacos-1.0.1.tar.gz"
#url = "https://gigamove.rz.rwth-aachen.de/d/id/CTzyApN76MXMJ6/dd/100" % version

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
    sys.exit("ERROR Failed to download source code with 'curl' or 'wget'")
  return

# parse args
parser = ArgumentParser(description = help_text)

parser.add_argument("--build", help="download and build the Scafacos library.", default=True)
parser.add_argument("--path", help="specify folder of existing Scafacos installation.")
parser.add_argument("--version", help="set version of Voro++ to download and build (default scafacos-1.0.1).", default=version)

homepath = "."

args = parser.parse_args()

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,args.version)

if args.path is not None:
    if not os.path.isdir(args.path): sys.exit("ERROR Scafacos path does not exist")
    homedir = args.path

if args.build is None and args.path is None:
    sys.exit("ERROR Cannot use -b and -p flag at the same time")

# download and unpack Scafacos tarball

if args.build:
  print("Downloading Scafacos ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,args.version))

  print("Unpacking Scafacos tarball ...")
  if os.path.exists("%s/%s" % (homepath,args.version)):
    cmd = 'rm -rf "%s/%s"' % (homepath,args.version)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,args.version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s.tar.gz" % (homepath,args.version))
  if os.path.basename(homedir) != args.version:
    if os.path.exists(homedir):
      cmd = 'rm -rf "%s"' % homedir
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    os.rename("%s/%s" % (homepath,args.version),homedir)

# build Scafacos

  print("Building Scafacos ...")
  cmd = 'cd "%s"; ./configure --prefix="`pwd`/build" --disable-doc --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m --with-internal-fftw --with-internal-pfft --with-internal-pnfft CC=mpicc FC=mpif90 CXX=mpicxx F77= > log.txt; make -j; make install' % homedir
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))

# create 2 links in lib/scafacos to Scafacos include/lib dirs

print("Creating links to Scafacos include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
cmd = 'ln -s "%s/build/include" includelink' % homedir
subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
cmd = 'ln -s "%s/build/lib" liblink' % homedir
subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
