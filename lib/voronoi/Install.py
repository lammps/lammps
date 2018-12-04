#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess
from argparse import ArgumentParser

# help message

help_text = """
Syntax from src dir: make lib-voronoi args="-b"
                 or: make lib-voronoi args="-p /usr/local/voro++-0.4.6"
                 or: make lib-voronoi args="-b -v voro++-0.4.6"
Syntax from lib dir: python Install.py -b -v voro++-0.4.6
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/voro++-0.4.6

specify one or more options, order does not matter

  -b = download and build the Voro++ library
  -p = specify folder of existing Voro++ installation
  -v = set version of Voro++ to download and build (default voro++-0.4.6)

Example:

make lib-voronoi args="-b"   # download/build in lib/voronoi/voro++-0.4.6
make lib-voronoi args="-p $HOME/voro++-0.4.6" # use existing Voro++ installation in $HOME/voro++-0.4.6
"""

# settings

version = "voro++-0.4.6"
url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version

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
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling wget failed with: %s" % e.output.decode('UTF-8'))

  if not success:
    sys.exit("Failed to download source code with 'curl' or 'wget'")
  return

# parse args
parser = ArgumentParser(description=help_text)

parser.add_argument("--build", help="download and build the Voro++ library.")
parser.add_argument("--path", help="specify folder of existing Voro++ installation.")
parser.add_argument("--version", help="set version of Voro++ to download and build (default voro++-0.4.6).", default=version)

args = parser.parse_args()

homepath = "."
homedir = args.version

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,args.version)

if args.path is not None:
    if not os.path.isdir(args.path): sys.exit("Voro++ path does not exist")
    homedir = args.path

if args.build is not None and args.path is not None:
    sys.exit("Cannot use -b and -p flag at the same time")

if args.build is None and args.path is None:
    sys.exit("Have to use either -b or -p flag")

# download and unpack Voro++ tarball

if args.build:
  print("Downloading Voro++ ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,args.version))

  print("Unpacking Voro++ tarball ...")
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

# build Voro++

  print("Building Voro++ ...")
  cmd = 'cd "%s"; make CXX=g++ CFLAGS="-fPIC -O3"' % homedir
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))

# create 2 links in lib/voronoi to Voro++ src dir

print("Creating links to Voro++ include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
cmd = 'ln -s "%s/src" includelink' % homedir
subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
cmd = 'ln -s "%s/src" liblink' % homedir
subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
