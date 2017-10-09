#!/usr/bin/env python

# install.py tool to download, compile, and setup the kim-api library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-kim args="-b -v version  -a kim-name"
                 or: make lib-kim args="-b -a everything"
                 or: make lib-kim args="-n -a kim-name"
                 or: make lib-kim args="-p /usr/local/open-kim -a kim-name"
Syntax from lib dir: python Install.py -b -v version  -a kim-name
                 or: python Install.py -b -a everything
                 or: python Install.py -n -a kim-name
                 or: python Install.py -p /usr/local/open-kim -a kim-name

specify one or more options, order does not matter

  -v = version of KIM API library to use
       default = kim-api-v1.9.2 (current as of Oct 2017)
  -b = download and build base KIM API library with example Models
       this will delete any previous installation in the current folder
  -n = do NOT download and build base KIM API library.
       Use an existing installation
  -p = specify location of KIM API installation (implies -n)
  -a = add single KIM model or model driver with kim-name
       to existing KIM API lib (see example below).
       If kim-name = everything, then rebuild KIM API library with
       *all* available OpenKIM Models (make take a long time).
  -vv = be more verbose about what is happening while the script runs

Examples:

make lib-kim args="-b" # install KIM API lib with only example models
make lib-kim args="-a Glue_Ercolessi_Adams_Al__MO_324507536345_001"  # Ditto plus one model
make lib-kim args="-b -a everything"   # install KIM API lib with all models
make lib-kim args="-n -a EAM_Dynamo_Ackland_W__MO_141627196590_002"   # only add one model or model driver

See the list of KIM model drivers here:
https://openkim.org/kim-items/model-drivers/alphabetical

See the list of all KIM models here:
https://openkim.org/kim-items/models/by-model-drivers

See the list of example KIM models included by default here:
https://openkim.org/kim-api
in the "What is in the KIM API source package?" section
"""

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

thisdir = os.environ['PWD']
version = "kim-api-v1.9.2"

buildflag = False
everythingflag = False
addflag = False
verboseflag = False
pathflag = False

iarg = 0
while iarg < len(args):
  if args[iarg] == "-v":
    if iarg+2 > len(args): error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  elif args[iarg] == "-n":
    buildflag = False
    iarg += 1
  elif args[iarg] == "-p":
    if iarg+2 > len(args): error()
    kimdir = fullpath(args[iarg+1])
    pathflag = True
    buildflag = False
    iarg += 2
  elif args[iarg] == "-a":
    addflag = True
    if iarg+2 > len(args): error()
    addmodelname = args[iarg+1]
    if addmodelname == "everything":
      buildflag = True
      everythingflag = True
      addflag = False
    iarg += 2
  elif args[iarg] == "-vv":
    verboseflag = True
    iarg += 1
  else: error()

thisdir = os.path.abspath(thisdir)
url = "https://s3.openkim.org/kim-api/%s.txz" % version

# set KIM API directory

if pathflag:
  if not os.path.isdir(kimdir):
    print("\nkim-api is not installed at %s" % kimdir)
    error()

  # configure LAMMPS to use existing kim-api installation
  with open("%s/Makefile.KIM_DIR" % thisdir, 'w') as mkfile:
    mkfile.write("KIM_INSTALL_DIR=%s\n\n" % kimdir)
    mkfile.write(".DUMMY: print_dir\n\n")
    mkfile.write("print_dir:\n")
    mkfile.write("	@printf $(KIM_INSTALL_DIR)\n")

  with open("%s/Makefile.KIM_Config" % thisdir, 'w') as cfgfile:
    cfgfile.write("include %s/lib/kim-api/Makefile.KIM_Config" % kimdir)

  print("Created %s/Makefile.KIM_DIR\n  using %s" % (thisdir,kimdir))
else:
  kimdir = os.path.join(os.path.abspath(thisdir), "installed-" + version)

# download KIM tarball, unpack, build KIM
if buildflag:

  # check to see if an installed kim-api already exists and wipe it out.

  if os.path.isdir(kimdir):
    print("kim-api is already installed at %s.\nRemoving it for re-install" % kimdir)
    cmd = 'rm -rf "%s"' % kimdir
    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  # configure LAMMPS to use kim-api to be installed

  with open("%s/Makefile.KIM_DIR" % thisdir, 'w') as mkfile:
    mkfile.write("KIM_INSTALL_DIR=%s\n\n" % kimdir)
    mkfile.write(".DUMMY: print_dir\n\n")
    mkfile.write("print_dir:\n")
    mkfile.write("	@printf $(KIM_INSTALL_DIR)\n")

  with open("%s/Makefile.KIM_Config" % thisdir, 'w') as cfgfile:
    cfgfile.write("include %s/lib/kim-api/Makefile.KIM_Config" % kimdir)

  print("Created %s/Makefile.KIM_DIR\n  using %s" % (thisdir,kimdir))

  # download entire kim-api tarball

  print("Downloading kim-api tarball ...")
  geturl(url,"%s/%s.txz" % (thisdir,version))
  print("Unpacking kim-api tarball ...")
  cmd = 'cd "%s"; rm -rf "%s"; tar -xJvf %s.txz' % (thisdir,version,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  # configure kim-api

  print("Configuring kim-api ...")
  cmd = 'cd "%s/%s"; ./configure --prefix="%s"' % (thisdir,version,kimdir)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  # build kim-api
  print("Building kim-api ...")
  cmd = 'cd "%s/%s"; make' % (thisdir,version)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print(txt.decode("UTF-8"))

  # install kim-api

  print("Installing kim-api ...")
  cmd = 'cd "%s/%s"; make install' % (thisdir,version)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print(txt.decode("UTF-8"))

  # remove source files

  print("Building and installing example Models")
  cmd = 'cd "%s/%s/examples"; make model-drivers-all-system' % (thisdir,version)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print (txt.decode("UTF-8"))
  cmd = 'cd "%s/%s/examples"; make models-all-system' % (thisdir,version)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print (txt.decode("UTF-8"))

  print("Removing kim-api source and build files ...")
  cmd = 'cd "%s"; rm -rf %s; rm -rf %s.txz' % (thisdir,version,version)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  # add all OpenKIM models, if desired
  if everythingflag:
    print("Adding all OpenKIM models, this will take a while ...")
    cmd = '%s/bin/kim-api-v1-collections-management install system OpenKIM' % (kimdir)
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    if verboseflag: print(txt.decode("UTF-8"))

# add single OpenKIM model
if addflag:

  if not os.path.isdir(kimdir):
    print("\nkim-api is not installed")
    error()

  # download single model
  cmd = '%s/bin/kim-api-v1-collections-management install system %s' % (kimdir, addmodelname)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print (txt.decode("UTF-8"))
