#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-plumed args="-b"
                 or: make lib-plumed args="-b -v 2.4.2"
Syntax from lib dir: python Install.py -b -v 2.4.2
                 or: python Install.py -b

specify one or more options, order does not matter

  -b = download and build the plumed2 library
  -v = set version of Voro++ to download and build (default: latest stable version)

Example:

make lib-plumed args="-b"   # download/build in lib/plumed/plumed2
"""

# settings

version = "2.4.2"

# Add known checksums for different PLUMED versions and use them to validate the download

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
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling wget failed with: %s" % e.output.decode('UTF-8'))

  if not success:
    error("Failed to download source code with 'curl' or 'wget'")
  return

# Here add function to check fsum

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

homepath = "."

buildflag = False
suffixflag = False
linkflag = True

iarg = 0
while iarg < nargs:
  if args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)

# download and unpack plumed  tarball

if buildflag:
  url = "https://github.com/plumed/plumed2/archive/v%s.tar.gz" % version
  filename = "v%s.tar.gz" %version
  print("Downloading plumed  ...")
  geturl(url,filename) 

  print("Unpacking plumed tarball ...")
  if os.path.exists("%s/%s" % (homepath,version)):
    cmd = 'rm -rf "%s/%s"' % (homepath,version)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar -xzvf v%s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/v%s.tar.gz" % (homepath,version))

# build plumed
 
if buildflag:
   print("Building plumed ...")
   cmd = 'cd %s/plumed2-%s; ./configure ; make' % (homepath,version)
   txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
   print(txt.decode('UTF-8'))
# 
# # create 2 links in lib/voronoi to Voro++ src dir
# 
# if linkflag:
#   print("Creating links to Voro++ include and lib files")
#   if os.path.isfile("includelink") or os.path.islink("includelink"):
#     os.remove("includelink")
#   if os.path.isfile("liblink") or os.path.islink("liblink"):
#     os.remove("liblink")
#   cmd = 'ln -s "%s/src" includelink' % homedir
#   subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
#   cmd = 'ln -s "%s/src" liblink' % homedir
#   subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
