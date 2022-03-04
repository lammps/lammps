#!/usr/bin/env python

# install.py tool to do a generic build of a library
# soft linked to by many of the lib/Install.py files
# used to automate the steps described in the corresponding lib/README

from __future__ import print_function
import sys,os,subprocess
import glob

sys.path.append('..')
from install_helpers import checkmd5sum

# help message

help = """
Syntax from src dir: make lib-libname args="-m machine -e suffix"
Syntax from lib dir: python Install.py -m machine -e suffix

libname = name of lib dir (e.g. atc, h5md, meam, poems, etc)
specify -m and optionally -e, order does not matter

  -m = peform a clean followed by "make -f Makefile.machine"
       machine = suffix of a lib/Makefile.* file
  -e = set EXTRAMAKE variable in Makefile.machine to Makefile.lammps.suffix
       does not alter existing Makefile.machine

Examples:

make lib-poems args="-m serial" # build POEMS lib with same settings as in the serial Makefile in src
make lib-colvars args="-m mpi"  # build COLVARS lib with same settings as in the mpi Makefile in src
make lib-meam args="-m ifort"   # build MEAM lib with custom Makefile.ifort (using Intel Fortran)
"""

# settings

version = "1.2.9"
url = "https://github.com/MolSSI-MDI/MDI_Library/archive/v%s.tar.gz" % version

# known checksums for different MDI versions. used to validate the download.
checksums = { \
              '1.2.7' : '2f3177b30ccdbd6ae28ea3bdd5fed0db', \
              '1.2.9' : 'ddfa46d6ee15b4e59cfd527ec7212184', \
              }

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

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

machine = None
extraflag = 0

iarg = 0
while iarg < nargs:
  if args[iarg] == "-m":
    if iarg+2 > nargs: error()
    machine = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-e":
    if iarg+2 > nargs: error()
    extraflag = 1
    suffix = args[iarg+1]
    iarg += 2
  else: error()

# set lib from working dir

cwd = os.getcwd()
lib = os.path.basename(cwd)

# download and unpack MDI_Library tarball

homepath = "."
homedir = "%s/MDI_Library" % homepath

print("Downloading MDI_Library ...")
mditar = "%s/v%s.tar.gz" % (homepath,version)
geturl(url, mditar)

# verify downloaded archive integrity via md5 checksum, if known.
if version in checksums:
  if not checkmd5sum(checksums[version], mditar):
    sys.exit("Checksum for MDI library does not match")

print("Unpacking MDI_Library tarball ...")
if os.path.exists("%s/v%s" % (homepath,version)):
  cmd = 'rm -rf "%s/v%s"' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
cmd = 'cd "%s"; tar -xzvf v%s.tar.gz' % (homepath,version)
subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
os.remove("%s/v%s.tar.gz" % (homepath,version))
if os.path.basename(homedir) != version:
  if os.path.exists(homedir):
    cmd = 'rm -rf "%s"' % homedir
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.rename("%s/MDI_Library-%s" % (homepath,version),homedir)

# create Makefile.auto as copy of Makefile.machine
# reset EXTRAMAKE if requested

if not os.path.exists("Makefile.%s" % machine):
  error("lib/%s/Makefile.%s does not exist" % (lib,machine))

lines = open("Makefile.%s" % machine,'r').readlines()
fp = open("Makefile.auto",'w')

has_extramake = False
for line in lines:
  words = line.split()
  if len(words) == 3 and words[0] == "EXTRAMAKE" and words[1] == '=':
    has_extramake = True
    if extraflag:
      line = line.replace(words[2],"Makefile.lammps.%s" % suffix)
  fp.write(line)

fp.close()

# make the library via Makefile.auto optionally with parallel make

try:
  import multiprocessing
  n_cpus = multiprocessing.cpu_count()
except:
  n_cpus = 1

print("Building lib%s.so ..." % lib)
cmd = "make -f Makefile.auto clean; make -f Makefile.auto -j%d" % n_cpus
txt = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
print(txt.decode('UTF-8'))

# create 2 links in lib/mdi to MDI Library src dir

print("Creating links to MDI Library include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
os.symlink(os.path.join(homedir, 'MDI_Library'), 'includelink')
os.symlink(os.path.join(homepath, 'build', 'MDI_Library'), 'liblink')

# Append the -rpath option to Makefile.lammps

dir_path = os.path.dirname(os.path.realpath(__file__))
rpath_option = "-Wl,-rpath=" + str(dir_path) + "/liblink"
makefile_lammps = open(str(dir_path) + "/Makefile.lammps", "a")
makefile_lammps.write(str(rpath_option) + "\n")
makefile_lammps.close()


shared_files = glob.glob( os.path.join( homepath, "liblink", "lib%s.so*" % lib) )
if len(shared_files) > 0:
  print("Build was successful")
else:
  error("Build of lib/%s/lib%s.so was NOT successful" % (lib,lib))
if has_extramake and not os.path.exists("Makefile.lammps"):
  print("lib/%s/Makefile.lammps was NOT created" % lib)
