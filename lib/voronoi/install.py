#!usr/local/python

# install.py tool to download, unpack, build, and link to the Voro++ library
# used to automate the steps described in the README file in this dir

import sys,os,re,urllib,commands

help = """
Syntax: install.py -d dir -v version -g -b -i installdir -l incdir libdir
        specify one or more options, order does not matter
        -d = dir to download tarball to, unpack tarball in, perform build in
             dir will be created if it doesn't exist (only last level)
             default = this dir
        -v = version of Voro++ to download and work with
             default = voro++-0.4.6 (current as of Jan 2015)
        -g = download (grab) tarball from
             http://math.lbl.gov/voro++/download/dir/version
        -b = build Voro++ by invoking "make" in its home dir
             no default
        -i = install Voro++ by invoking "make install" in its home dir
             installdir arg is optional:
               if not specified, installs at PREFIX defined in config.mk file
               if specified, will overwrite PREFIX and install there
             if PREFIX starts with /usr, will invoke "sudo make install"
        -l = create two links to incdir and libdir
             incdir and libdir are optional (specify neither or both):
               if specified, includelink and liblink are to those two dirs
                 these are dirs where Voro++ include files and lib file are
               if not specified and no install, links are to Voro++ src dir
               if not specified and install performed,
                 links are to include and lib dirs under PREFIX
"""

def error():
  print help
  sys.exit()
  
# parse args

args = sys.argv

if len(args) == 1: error()

dir = "."
version = "voro++-0.4.6"
grabflag = 0
buildflag = 0
installflag = 0
linkflag = 0

iarg = 1
while iarg < len(args):
  if args[iarg] == "-d":
    if iarg+2 > len(args): error()
    dir = args[iarg+1]
    iarg += 2  
  elif args[iarg] == "-v":
    if iarg+2 > len(args): error()
    version = args[iarg+1]
    iarg += 2  
  elif args[iarg] == "-g":
    grabflag = 1
    iarg += 1
  elif args[iarg] == "-b":
    buildflag = 1
    iarg += 1
  elif args[iarg] == "-i":
    installflag = 1
    if iarg+1 == len(args) or args[iarg+1][0] == '-':
      installdir = ""
      iarg += 1
    else:
      if iarg+2 > len(args): error()
      installdir = args[iarg+1]
      iarg += 2  
  elif args[iarg] == "-l":
    linkflag = 1
    if iarg+1 == len(args) or args[iarg+1][0] == '-' or \
          iarg+2 == len(args) or args[iarg+2][0] == '-':
      includedir = libdir = ""
      iarg += 1
    else:
      if iarg+3 > len(args): error()
      includedir = args[iarg+1]
      libdir = args[iarg+2]
      iarg += 3
  else: error()

dir = os.path.abspath(dir)
url = "http://math.lbl.gov/voro++/download/dir/%s.tar.gz" % version

# create dir if does not exist

if not os.path.isdir(dir):
  if os.path.isfile(dir):
    print "ERROR: Dir already exists as file"
    sys.exit()
  os.mkdir(dir)
  if not os.path.isdir(dir):
    print "ERROR: Unable to create dir"
    sys.exit()

# download and unpack tarball

if grabflag:
  print "Downloading Voro++ tarball ..."
  urllib.urlretrieve(url,"%s/%s.tar.gz" % (dir,version))
  print "Unpacking Voro++ tarball ..."
  cmd = "cd %s; tar zxvf %s.tar.gz" % (dir,version)
  txt = commands.getoutput(cmd)

# build Voro++ in its dir

if buildflag:
  print "Building Voro++ ..."
  cmd = "cd %s/%s; make" % (dir,version)
  txt = commands.getoutput(cmd)
  print txt

# install Voro++
# if installdir set, overwrite PREFIX var in its config.mk file
# if PREFIX var starts with /usr, invoke sudo make install, else make install
  
if installflag:
  print "Installing Voro++ ..."
  if installdir:
    txt = open("%s/%s/config.mk" % (dir,version),'r').read()
    txt = re.sub("PREFIX=.*?\n","PREFIX=%s\n" % installdir,txt)
    open("%s/%s/config.mk" % (dir,version),'w').write(txt)
    print "TXT:",txt
  txt = open("%s/%s/config.mk" % (dir,version),'r').read()
  var = re.findall("PREFIX=.*?\n",txt)
  prefix = var[0].split('=')[1].strip()
  if prefix.startswith("/usr"):
    cmd = "cd %s/%s; sudo make install" % (dir,version)
  else:
    cmd = "cd %s/%s; make install" % (dir,version)
  txt = commands.getoutput(cmd)
  print txt
  
# create links in this dir to Voro++ include and lib files

if linkflag:
  print "Creating links to Voro++ include and lib files"
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  if includedir:
    cmd = "ln -s %s includelink" % includedir
    txt = commands.getoutput(cmd)
    cmd = "ln -s %s liblink" % linkdir
    txt = commands.getoutput(cmd)
  elif not installflag:
    cmd = "ln -s %s/%s/src includelink" % (dir,version)
    txt = commands.getoutput(cmd)
    cmd = "ln -s %s/%s/src liblink" % (dir,version)
    txt = commands.getoutput(cmd)
  else:
    cmd = "ln -s %s/include includelink" % prefix
    txt = commands.getoutput(cmd)
    cmd = "ln -s %s/lib liblink" % prefix
    txt = commands.getoutput(cmd)
