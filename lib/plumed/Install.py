#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the plumed2 library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,hashlib

# help message

help = """
Syntax from src dir: make lib-plumed args="-b"
                 or: make lib-plumed args="-b -v 2.4.3"
                 or: make lib-plumed args="-p /usr/local/plumed2 -m shared"

Syntax from lib dir: python Install.py -b -v 2.4.3
                 or: python Install.py -b
                 or: python Install.py -p /usr/local/plumed2 -m shared

specify one or more options, order does not matter

  -b = download and build the plumed2 library
  -v = set version of plumed2 to download and build (default: 2.4.3)
  -p = specify folder of existing plumed2 installation
  -m = set plumed linkage mode: static (default), shared, or runtime

Example:

make lib-plumed args="-b"   # download/build in lib/plumed/plumed2
make lib-plumed args="-p $HOME/plumed2 -m shared" # use existing Plumed2 installation in $HOME/plumed2
"""

# settings

version = "2.4.3"
mode = "static"

# known checksums for different PLUMED versions. used to validate the download.
checksums = { \
        '2.4.2' : '88188743a6e03ef076e5377d03ebb0e7', \
        '2.4.3' : 'b1be7c48971627febc11c61b70767fc5', \
        '2.5b'  : 'e341bdef469be1da058b8a0b97a3db22', \
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
  elif args[iarg] == "-m":
    if iarg+2 > nargs: error()
    mode = args[iarg+1]
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

if ((mode != "static") and (mode != "shared") and (mode != "runtime")):
    error("Unknown linkage mode '%s' for Plumed" % mode)

# download and unpack plumed2 tarball

if buildflag:
  url = "https://github.com/plumed/plumed2/releases/download/v%s/plumed-src-%s.tgz" % (version,version)
  filename = "plumed-src-%s.tar.gz" %version
  print("Downloading plumed  ...")
  geturl(url,filename)

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version],filename):
      error("Checksum for plumed2 library does not match")

  print("Unpacking plumed2 source tarball ...")
  if os.path.exists("%s/plumed-%s" % (homepath,version)):
    cmd = 'rm -rf "%s/plumed-%s"' % (homepath,version)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if os.path.exists("%s/plumed2" % (homepath)):
    cmd = 'rm -rf "%s/plumed2"' % (homepath)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'cd "%s"; tar -xzvf %s' % (homepath,filename)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s/%s" % (homepath,filename))

# build plumed
 
if buildflag:
   print("Building plumed ...")
   cmd = 'cd %s/plumed-%s; ./configure --prefix=%s/plumed2 --enable-static-patch ; make ; make install' % (homepath,version,homepath)
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
  cmd = 'ln -s "%s/include" includelink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/lib" liblink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if os.path.isfile("Makefile.lammps.%s" % mode):
    print("Creating Makefile.lammps")
    cmd = 'echo PLUMED_LIBDIR="%s/lib" > Makefile.lammps; cat liblink/plumed/src/lib/Plumed.inc.%s Makefile.lammps.%s >> Makefile.lammps' % (homedir,mode,mode)
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

