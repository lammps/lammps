#!/usr/bin/env python

# install.pa tool to setup the kim-api library
# used to automate the steps described in the README file in this dir
from __future__ import print_function
import sys,os,re,subprocess
try: from urllib.request import urlretrieve as geturl
except: from urllib import urlretrieve as geturl

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

  # check to see if an installed kim-api already exists and wipe it out.

  if os.path.isdir(dir):
    print("kim-api is already installed at %s.\nRemoving it for re-install" % dir)
    cmd = "rm -rf %s" % dir
    subprocess.check_output(cmd,shell=True)

  # configure LAMMPS to use kim-api to be installed

  with open("%s/Makefile.KIM_DIR" % thisdir, 'w') as mkfile:
    mkfile.write("KIM_INSTALL_DIR=%s\n\n" % dir)
    mkfile.write(".DUMMY: print_dir\n\n")
    mkfile.write("print_dir:\n")
    mkfile.write("	@printf $(KIM_INSTALL_DIR)\n")

  with open("%s/Makefile.KIM_Config" % thisdir, 'w') as cfgfile:
    cfgfile.write("include %s/lib/kim-api/Makefile.KIM_Config" % dir)

  print("Created %s/Makefile.KIM_DIR : using %s" % (thisdir,dir))

  # download entire kim-api tarball

  print("Downloading kim-api tarball ...")
  geturl(url,"%s/%s.tgz" % (thisdir,version))
  print("Unpacking kim-api tarball ...")
  cmd = "cd %s; rm -rf %s; tar zxvf %s.tgz" % (thisdir,version,version)
  subprocess.check_output(cmd,shell=True)

  # configure kim-api

  print("Configuring kim-api ...")
  cmd = "cd %s/%s; ./configure --prefix='%s'" % (thisdir,version,dir)
  subprocess.check_output(cmd,shell=True)

  # build kim-api

  print("Configuring model : %s" % modelname)
  if modelname == "none":
    cmd = "cd %s/%s; make add-examples" % (thisdir,version)
  else:
    if modelname == "OpenKIM":
      print("configuring all OpenKIM models, this will take a while ...")
    cmd = "cd %s/%s; make add-examples; make add-%s" % \
        (thisdir,version,modelname)
  subprocess.check_output(cmd,shell=True)

  print("Building kim-api ...")
  cmd = "cd %s/%s; make" % (thisdir,version)
  subprocess.check_output(cmd,shell=True)

  # install kim-api

  print("Installing kim-api ...")
  cmd = "cd %s/%s; make install" % (thisdir,version)
  subprocess.check_output(cmd,shell=True)

  cmd = "cd %s/%s; make install-set-default-to-v1" %(thisdir,version)
  subprocess.check_output(cmd,shell=True)

  # remove source files

  print("Removing kim-api source and build files ...")
  cmd = "cd %s; rm -rf %s; rm -rf %s.tgz" % (thisdir,version,version)
  subprocess.check_output(cmd,shell=True)

# add a single model (and possibly its driver) to existing KIM installation

if addflag:

  # get location of installed kim-api

  if not os.path.isfile("%s/Makefile.KIM_DIR" % thisdir):
    print("kim-api is not installed")
    error()
  else:
    cmd = "cd %s; make -f Makefile.KIM_DIR print_dir" % thisdir
    dir = subprocess.check_output(cmd,shell=True)[1]

  # download single model
  # try first via urllib
  # if fails (probably due to no SSL support), use wget

  print("Downloading item tarball ...")
  url = "https://openkim.org/download/%s.tgz" % addmodelname
  geturl(url,"%s/%s.tgz" % (thisdir,addmodelname))

  print("Unpacking item tarball ...")
  cmd = "cd %s; tar zxvf %s.tgz" % (thisdir,addmodelname)
  subprocess.check_output(cmd,shell=True)

  print("Building item ...")
  cmd = "cd %s/%s; make; make install" %(thisdir,addmodelname)
  try:
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  except subprocess.CalledProcessError as e:

    # Error: but first, check to see if it needs a driver
    firstRunOutput = e.output.decode()

    cmd = "cd %s/%s; make kim-item-type" % (thisdir,addmodelname)
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
    txt = txt.decode().strip()
    if txt == "ParameterizedModel":

      # Get and install driver

      cmd = "cd %s/%s; make model-driver-name" % (thisdir,addmodelname)
      txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      adddrivername = txt.decode().strip()
      print("First installing model driver: %s" % adddrivername)
      cmd = "cd %s; python Install.py -a %s" % (thisdir,adddrivername)
      try:
        txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      except subprocess.CalledProcessError as e:
        print(e.output)
        sys.exit()

      # now install the model that needed the driver

      print("Now installing model : %s" % addmodelname)
      cmd = "cd %s; python Install.py -a %s" % (thisdir,addmodelname)
      try:
        txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      except subprocess.CalledProcessError as e:
        print(e.output)
        sys.exit()
      print(txt.decode())
      sys.exit()
    else:
      print(firstRunOutput)
      print("Error, unable to build and install OpenKIM item: %s" \
            % addmodelname)
      sys.exit()

  # success the first time

  print(txt)
  print("Removing kim item source and build files ...")
  cmd = "cd %s; rm -rf %s; rm -rf %s.tgz" %(thisdir,addmodelname,addmodelname)
  subprocess.check_output(cmd,shell=True)
