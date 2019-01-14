#!/usr/bin/env python

# install.py tool to download, compile, and setup the kim-api library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,shutil
sys.path.append('..')
from install_helpers import fullpath,geturl
from argparse import ArgumentParser

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

thisdir = fullpath('.')
version = "kim-api-v1.9.5"

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

Examples:

make lib-kim args="-b" # install KIM API lib with only example models
make lib-kim args="-b -a Glue_Ercolessi_Adams_Al__MO_324507536345_001" # Ditto plus one model
make lib-kim args="-b -a everything"   # install KIM API lib with all models
make lib-kim args="-n -a EAM_Dynamo_Ackland_W__MO_141627196590_002"    # only add one model or model driver

See the list of KIM model drivers here:
https://openkim.org/kim-items/model-drivers/alphabetical

See the list of all KIM models here:
https://openkim.org/kim-items/models/by-model-drivers

See the list of example KIM models included by default here:
https://openkim.org/kim-api
in the "What is in the KIM API source package?" section
"""

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build base KIM API library with example Models.")
pgroup.add_argument("-n", "--nobuild", action="store_true",
                    help="use the previously downloaded and compiled base KIM API.")
pgroup.add_argument("-p", "--path",
                    help="specify location of existing KIM API installation.")
parser.add_argument("-v", "--version", default=version,
                    help="set version of KIM API library to download and build (default: %s)" % version)
parser.add_argument("-a", "--add",
                    help="add single KIM model or model driver. If adding 'everything', then all available OpenKIM models are added (may take a long time)")
parser.add_argument("-vv", "--verbose", action="store_true",
                    help="be more verbose about is happening while this script runs")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if args.build == False and not args.path and args.nobuild == False:
  parser.print_help()
  sys.exit(help)

buildflag = args.build
pathflag = args.path != None
addflag = args.add != None
addmodelname = args.add
everythingflag = False
if addflag and addmodelname == "everything":
  everythingflag = True
  buildflag = True
verboseflag = args.verbose

if pathflag:
  buildflag = False
  kimdir = args.path
  if not os.path.isdir(kimdir): sys.exit("KIM API path %s does not exist" % kimdir)
  kimdir = fullpath(kimdir)

url = "https://s3.openkim.org/kim-api/%s.txz" % version

# set KIM API directory

if pathflag:
  # configure LAMMPS to use existing kim-api installation
  with open("%s/Makefile.KIM_DIR" % thisdir, 'w') as mkfile:
    mkfile.write("KIM_INSTALL_DIR=%s\n\n" % kimdir)
    mkfile.write(".DUMMY: print_dir\n\n")
    mkfile.write("print_dir:\n")
    mkfile.write("	@printf $(KIM_INSTALL_DIR)\n")

  print("Created %s/Makefile.KIM_DIR\n  using %s" % (thisdir,kimdir))
else:
  kimdir = os.path.join(os.path.abspath(thisdir), "installed-" + version)
  if args.nobuild and not os.path.isdir(kimdir):
    sys.exit("Cannot use -n/--nobuild without first building the KIM API with -b")

# download KIM tarball, unpack, build KIM
if buildflag:

  # check to see if an installed kim-api already exists and wipe it out.

  if os.path.isdir(kimdir):
    print("kim-api is already installed at %s.\nRemoving it for re-install" % kimdir)
    shutil.rmtree(kimdir)

  # configure LAMMPS to use kim-api to be installed

  with open("%s/Makefile.KIM_DIR" % thisdir, 'w') as mkfile:
    mkfile.write("KIM_INSTALL_DIR=%s\n\n" % kimdir)
    mkfile.write(".DUMMY: print_dir\n\n")
    mkfile.write("print_dir:\n")
    mkfile.write("	@printf $(KIM_INSTALL_DIR)\n")

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

  makefile_path = os.path.join(thisdir, "Makefile.KIM_DIR")
  if os.path.isfile(makefile_path):
    cmd = 'make --no-print-directory -f %s print_dir' % makefile_path
    kimdir = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  if not os.path.isdir(kimdir):
    sys.exit("\nkim-api is not installed")

  # download single model
  cmd = '%s/bin/kim-api-v1-collections-management install system %s' % (kimdir.decode("UTF-8"), addmodelname)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print (txt.decode("UTF-8"))
