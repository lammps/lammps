#!/usr/bin/env python

# install.pa tool to setup the kim-api library
# used to automate the steps described in the README file in this dir
from __future__ import print_function
import sys,os,re,urllib,subprocess

help = """
Syntax from src dir: make lib-kim args="-v version -b kim-install-dir kim-name -a kim-name"
Syntax from lib dir: python Install.py -v version -b kim-install-dir kim-name -a kim-name

specify one or more options, order does not matter

  -v = version of KIM API library to use
       default = kim-api-v1.8.2 (current as of June 2017)
  -b = download and build KIM API library with KIM models
       kim-dir = where to install/build the KIM API library
         use "." to install in lib/kim
       kim-name = none to install only the example KIM models
       kim-name = KIM model name (see example below) + examples
       kim-name = OpenKIM to install all models
         from the openkim.org site (this can take 30 minutes or more)
  -a = add single KIM model or model driver with kim-name
       to existing KIM API lib (see example below)

Examples:

make lib-kim args="-b . none"   # install KIM API lib with only example models
make lib-kim args="-b . Glue_Ercolessi_Adams_Al__MO_324507536345_001"  # ditto plus one model
make lib-kim args="-b . OpenKIM"   # install KIM API lib with all models
make lib-kim args="-a EAM_Dynamo_Ackland_W__MO_141627196590_002"   # add one model or model driver

See the list of KIM model drivers here:
https://openkim.org/kim-items/model-drivers/alphabetical

See the list of all KIM models here:
https://openkim.org/kim-items/models/by-model-drivers

See the list of example KIM models included by default here:
https://openkim.org/kim-api
in the "What is in the KIM API source package?" section
"""

def error():
  print(help)
  sys.exit()

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

thisdir = os.environ['PWD']
version = "kim-api-v1.8.2"

buildflag = False
addflag = False

iarg = 0
while iarg < len(args):
  if args[iarg] == "-v":
    if iarg+2 > len(args): error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    if iarg+3 > len(args): error()
    dir = args[iarg+1]
    modelname = args[iarg+2]
    iarg += 3
  elif args[iarg] == "-a":
    addflag = True
    if iarg+2 > len(args): error()
    addmodelname = args[iarg+1]
    iarg += 2
  else: error()

thisdir = os.path.abspath(thisdir)
url = "https://s3.openkim.org/kim-api/%s.tgz" % version

# download KIM tarball, unpack, build KIM
# either in lib/kim or user-requested location

if buildflag:

  # set install directory
  dir = os.path.join(os.path.abspath(dir), "installed-" + version)

  # check to see if an installed kim-api already exists

  if os.path.isdir(dir):
     print("kim-api is already installed at %s" % dir)
     print("Must remove this directory in order to resintall at this location")
     sys.exit()

  # configure LAMMPS to use kim-api to be installed

  with open("%s/Makefile.KIM_DIR" % thisdir, 'w') as mkfile:
    mkfle.write("KIM_INSTALL_DIR=%s\n\n" % dir)
    mkfle.write(".DUMMY: print_dir\n\n")
    mkfle.write("print_dir:\n")
    mkfle.write("	@printf $(KIM_INSTALL_DIR)\n")

  with open("%s/Makefile.KIM_Config" % thisdir, 'w') as cfgfile:
    cfgfile.write("include %s/lib/kim-api/Makefile.KIM_Config" % dir)

  print("Created %s/Makefile.KIM_DIR : using %s" % (thisdir,dir))

  # download entire kim-api tarball
  # try first via urllib
  # if fails (probably due to no SSL support), use wget

  print("Downloading kim-api tarball ...")

  try: urllib.urlretrieve(url,"%s/%s.tgz" % (thisdir,version))
  except:
    cmd = "wget %s %s/%s.tgz" % (url,thisdir,version)
    txt = subprocess.getstatusoutput(cmd)
    print(txt[1])
    if not os.path.isfile("%s/%s.tgz" % (thisdir,version)):
      print("Both urllib.urlretrieve() and wget command failed to download")
      sys.exit()

  print("Unpacking kim-api tarball ...")
  cmd = "cd %s; rm -rf %s; tar zxvf %s.tgz" % (thisdir,version,version)
  txt = subprocess.getstatusoutput(cmd)
  if txt[0] != 0: error()

  # configure kim-api

  print("Configuring kim-api ...")
  cmd = "cd %s/%s; ./configure --prefix='%s'" % (thisdir,version,dir)
  txt = subprocess.getstatusoutput(cmd)
  print(txt[1])
  if txt[0] != 0: error()

  # build kim-api

  print("Configuring model : %s" % modelname)
  if modelname == "none":
    cmd = "cd %s/%s; make add-examples" % (thisdir,version)
  else:
    if modelname == "OpenKIM":
      print("configuring all OpenKIM models, this will take a while ...")
    cmd = "cd %s/%s; make add-examples; make add-%s" % \
        (thisdir,version,modelname)
  txt = subprocess.getstatusoutput(cmd)
  print(txt[1])
  if txt[0] != 0: error()

  print("Building kim-api ...")
  cmd = "cd %s/%s; make" % (thisdir,version)
  txt = subprocess.getstatusoutput(cmd)
  print(txt[1])
  if txt[0] != 0: error()

  # install kim-api

  print("Installing kim-api ...")
  cmd = "cd %s/%s; make install" % (thisdir,version)
  txt = subprocess.getstatusoutput(cmd)
  print(txt[1])
  if txt[0] != 0: error()

  cmd = "cd %s/%s; make install-set-default-to-v1" %(thisdir,version)
  txt = subprocess.getstatusoutput(cmd)
  print(txt[1])
  if txt[0] != 0: error()

  # remove source files

  print("Removing kim-api source and build files ...")
  cmd = "cd %s; rm -rf %s; rm -rf %s.tgz" % (thisdir,version,version)
  txt = subprocess.getstatusoutput(cmd)
  print(txt[1])
  if txt[0] != 0: error()

# add a single model (and possibly its driver) to existing KIM installation

if addflag:

  # get location of installed kim-api

  if not os.path.isfile("%s/Makefile.KIM_DIR" % thisdir):
    print("kim-api is not installed")
    error()
  else:
    cmd = "cd %s; make -f Makefile.KIM_DIR print_dir" % thisdir
    dir = subprocess.getstatusoutput(cmd)[1]

  # download single model
  # try first via urllib
  # if fails (probably due to no SSL support), use wget

  print("Downloading item tarball ...")

  url = "https://openkim.org/download/%s.tgz" % addmodelname

  try: urllib.urlretrieve(url,"%s/%s.tgz" % (thisdir,addmodelname))
  except:
    cmd = "wget %s %s/%s.tgz" % (url,thisdir,addmodelname)
    txt = subprocess.getstatusoutput(cmd)
    print(txt[1])
    if not os.path.isfile("%s/%s.tgz" % (thisdir,addmodelname)):
      print("Both urllib.urlretrieve() and wget command failed to download")
      sys.exit()

  print("Unpacking item tarball ...")
  cmd = "cd %s; tar zxvf %s.tgz" % (thisdir,addmodelname)
  txt = subprocess.getstatusoutput(cmd)
  if txt[0] != 0: error()

  print("Building item ...")
  cmd = "cd %s/%s; make; make install" %(thisdir,addmodelname)
  txt = subprocess.getstatusoutput(cmd)
  firstRunOutput = txt[1]
  if txt[0] != 0:
    # Error: but first, check to see if it needs a driver

    cmd = "cd %s/%s; make kim-item-type" % (thisdir,addmodelname)
    txt = subprocess.getstatusoutput(cmd)
    if txt[1] == "ParameterizedModel":

      # Get and install driver

      cmd = "cd %s/%s; make model-driver-name" % (thisdir,addmodelname)
      txt = subprocess.getstatusoutput(cmd)
      adddrivername = txt[1]
      print("First Installing model driver: %s" % adddrivername)
      cmd = "cd %s; python Install.py -a %s" % (thisdir,adddrivername)
      txt = subprocess.getstatusoutput(cmd)
      if txt[0] != 0:
         print(firstRunOutput)
         print(txt[1])
         error()
      else:
        print(txt[1])
      cmd = "cd %s; python Install.py -a %s" % (thisdir,addmodelname)
      txt = subprocess.getstatusoutput(cmd)
      print(txt[1])
      if txt[0] != 0:
        error()
    else:
      print(firstRunOutput)
      error()
  else:

    # success

    print(firstRunOutput)
    print("Removing kim item source and build files ...")
    cmd = "cd %s; rm -rf %s; rm -rf %s.tgz" %(thisdir,addmodelname,addmodelname)
    txt = subprocess.getstatusoutput(cmd)
    print(txt[1])
    if txt[0] != 0: error()
