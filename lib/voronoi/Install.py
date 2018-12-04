#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,shutil
sys.path.append('..')
from install_helpers import error,get_cpus,fullpath,which,geturl

from argparse import ArgumentParser

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = "voro++-0.4.6"
url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version

# extra help message

help = """
Syntax from src dir: make lib-voronoi args="-b"
                 or: make lib-voronoi args="-p /usr/local/voro++-0.4.6"
                 or: make lib-voronoi args="-b -v voro++-0.4.6"
Syntax from lib dir: python Install.py -b -v voro++-0.4.6
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/voro++-0.4.6

Example:

make lib-voronoi args="-b"   # download/build in lib/voronoi/voro++-0.4.6
make lib-voronoi args="-p $HOME/voro++-0.4.6" # use existing Voro++ installation in $HOME/voro++-0.4.6
"""

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the Voro++ library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing Voro++ installation")
parser.add_argument("-v", "--version", default=version,
                    help="set Voro++ version of Voro++ to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if args.build == False and not args.path:
  parser.print_help()
  sys.exit(help)

buildflag = args.build
pathflag = args.path != None
voropath = args.path
linkflag = True

homepath = fullpath(".")
homedir = "%s/%s" % (homepath,version)

if (pathflag):
    if not os.path.isdir(voropath): error("Voro++ path does not exist")
    homedir = voropath

# download and unpack Voro++ tarball

if buildflag:
  print("Downloading Voro++ ...")
  geturl(url,"%s/%s.tar.gz" % (homepath,version))

  print("Unpacking Voro++ tarball ...")
  if os.path.exists("%s/%s" % (homepath,version)):
    shutil.rmtree("%s/%s" % (homepath,version))
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s.tar.gz" % (homepath,version))
  if os.path.basename(homedir) != version:
    if os.path.exists(homedir):
      shutil.rmtree(homedir)
    os.rename("%s/%s" % (homepath,version),homedir)

# build Voro++

if buildflag:
  print("Building Voro++ ...")
  cmd = 'cd "%s"; make CXX=g++ CFLAGS="-fPIC -O3"' % homedir
  try:
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    print("Make failed with:\n %s" % e.output.decode('UTF-8'))
    sys.exit(1)

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
