#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the plumed2 library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,hashlib

# help message

help = """
Syntax from src dir: make lib-plumed args="-b"
                 or: make lib-plumed args="-b -v 2.4.3"
                 or: make lib-plumed args="-p /usr/local/plumed2-2.4.3"

Syntax from lib dir: python Install.py -b -v 2.4.3
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/plumed2-2.4.3

specify one or more options, order does not matter

  -b = download and build the plumed2 library
  -p = specify folder of existing plumed2 installation
  -v = set version of plumed2 to download and build (default: 2.4.3)

Example:

make lib-plumed args="-b"   # download/build in lib/plumed/plumed2
make lib-plumed args="-p $HOME/plumed-2.4.3" # use existing Plumed2 installation in $HOME/plumed-2.4.3
"""

# settings

version = "2.4.3"

# known checksums for different PLUMED versions. used to validate the download.
checksums = { \
        '2.4.2' : '0f66f24b4c763ae8b2f39574113e9935', \
        '2.4.3' : 'dc38de0ffd59d13950d8f1ef1ce05574', \
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

def checkmd5sum(md5sum,fname):
    with open(fname,'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(81920)
            if not data:
                break
            m.update(data)
    fh.close()
    return m.hexdigest() == md5sum

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
  if args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error()
    plumedpath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)

if (pathflag):
    if not os.path.isdir(plumedpath): error("Plumed2 path does not exist")
    homedir = plumedpath

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if (not buildflag and not pathflag):
    error("Have to use either -b or -p flag")

# download and unpack plumed2 tarball

if buildflag:
  url = "https://github.com/plumed/plumed2/archive/v%s.tar.gz" % version
  filename = "v%s.tar.gz" %version
  print("Downloading plumed  ...")
  geturl(url,filename)

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version],filename):
      error("Checksum for plumed2 library does not match")

  print("Unpacking plumed2 tarball ...")
  if os.path.exists("%s/plumed2-%s" % (homepath,version)):
    cmd = 'rm -rf "%s/plumed2-%s"' % (homepath,version)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if os.path.exists("%s/plumed2" % (homepath)):
    cmd = 'rm -rf "%s/plumed2"' % (homepath)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar -xzvf v%s.tar.gz' % (homepath,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/v%s.tar.gz" % (homepath,version))

# build plumed
 
if buildflag:
   print("Building plumed ...")
   cmd = 'cd %s/plumed2-%s; ./configure --prefix=%s/plumed2 ; make ; make install' % (homepath,version,homepath)
   txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
   print(txt.decode('UTF-8'))
# 
# create 2 links in lib/plumed to plumed2 installation dir

if linkflag:
  print("Creating links to plumed2 include and lib files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = 'ln -s "%s/plumed2/include" includelink' % homepath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/plumed2/lib" liblink' % homepath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if os.path.isfile("Makefile.lammps.static"):
    print("Creating Makefile.lammps")
    cmd = 'cat liblink/plumed/src/lib/Plumed.inc.static Makefile.lammps.static > Makefile.lammps'

