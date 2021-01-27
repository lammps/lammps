#!/usr/bin/env python

"""
Install.py tool to download, compile, and setup the kim-api library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

thisdir = fullpath('.')
version = "2.2.1"

# known checksums for different KIM-API versions. used to validate the download.
checksums = { \
        '2.1.2' : '6ac52e14ef52967fc7858220b208cba5', \
        '2.1.3' : '6ee829a1bbba5f8b9874c88c4c4ebff8', \
        '2.2.0' : 'e7f944e1593cffd7444679a660607f6c', \
        '2.2.1' : 'ae1ddda2ef7017ea07934e519d023dca', \
        }



# help message

HELP = """
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
make lib-kim args="-b -a EAM_ErcolessiAdams_1994_Al__MO_324507536345_002" # Ditto plus one model
make lib-kim args="-b -a everything"   # install KIM API lib with all models
make lib-kim args="-n -a EAM_Dynamo_Ackland_2003_W__MO_141627196590_005"    # only add one model or model driver

See the list of all KIM models here:
https://openkim.org/browse/models
"""

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build base KIM API library with example Models.")
pgroup.add_argument("-n", "--nobuild", action="store_true",
                    help="use the previously downloaded and compiled base KIM API.")
pgroup.add_argument("-p", "--path",
                    help="specify location of existing KIM API installation.")
parser.add_argument("-v", "--version", default=version, choices=checksums.keys(),
                    help="set version of KIM API library to download and build (default: %s)" % version)
parser.add_argument("-a", "--add",
                    help="add single KIM model or model driver. If adding 'everything', then all available OpenKIM models are added (may take a long time)")
parser.add_argument("-vv", "--verbose", action="store_true",
                    help="be more verbose about is happening while this script runs")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path and not args.nobuild:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
addflag = args.add is not None
addmodelname = args.add
everythingflag = False
if addflag and addmodelname == "everything":
  everythingflag = True
  buildflag = True
verboseflag = args.verbose
version = args.version

if pathflag:
  buildflag = False
  kimdir = args.path
  if not os.path.isdir(kimdir):
    sys.exit("KIM API path %s does not exist" % kimdir)
  kimdir = fullpath(kimdir)

url = "https://s3.openkim.org/kim-api/kim-api-%s.txz" % version

# set KIM API directory

if pathflag:
  # configure LAMMPS to use existing kim-api installation
  with open("%s/kim-prefix.txt" % thisdir, 'w') as pffile:
    pffile.write("%s" % kimdir)

  print("Created %s/kim-prefix.txt\n  using %s" % (thisdir,kimdir))
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

  with open("%s/kim-prefix.txt" % thisdir, 'w') as pffile:
    pffile.write("%s" % kimdir)

  print("Created %s/kim-prefix.txt\n  using %s" % (thisdir,kimdir))

  # download entire kim-api tarball

  print("Downloading kim-api tarball ...")
  filename = "kim-api-%s.txz" % version
  geturl(url, "%s/%s" % (thisdir, filename))

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], filename):
      sys.exit("Checksum for KIM-API library does not match")

  print("Unpacking kim-api tarball ...")
  cmd = 'cd "%s"; rm -rf "kim-api-%s"; tar -xJvf %s' % (thisdir, version, filename)
  subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

  # configure kim-api

  print("Configuring kim-api ...")
  cmd = 'cd "%s/kim-api-%s" && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX="%s" -DCMAKE_BUILD_TYPE=Release' % (thisdir,version,kimdir)
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if verboseflag: print(txt.decode("UTF-8"))

  # build kim-api
  print("Building kim-api ...")
  cmd = 'cd "%s/kim-api-%s/build" && make -j2' % (thisdir, version)
  txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  if verboseflag:
    print(txt.decode("UTF-8"))

  # install kim-api

  print("Installing kim-api ...")
  cmd = 'cd "%s/kim-api-%s/build" && make -j2 install' % (thisdir, version)
  txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  if verboseflag:
    print(txt.decode("UTF-8"))

  # remove source files

  print("Removing kim-api source and build files ...")
  cmd = 'cd "%s"; rm -rf kim-api-%s; rm -rf kim-api-%s.txz' % (thisdir, version, version)
  subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

  # add all OpenKIM models, if desired
  if everythingflag:
    print("Adding all OpenKIM models, this will take a while ...")
    cmd = '%s/bin/kim-api-collections-management install system OpenKIM' % (kimdir)
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    if verboseflag:
      print(txt.decode("UTF-8"))

# add single OpenKIM model
if addflag:

  pf_path = os.path.join(thisdir, "kim-prefix.txt")
  if os.path.isfile(pf_path):
    cmd = 'cat %s' % pf_path
    kimdir = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  if not os.path.isdir(kimdir):
    sys.exit("\nkim-api is not installed")

  # download single model
  cmd = '%s/bin/kim-api-collections-management install system %s' % (kimdir.decode("UTF-8"), addmodelname)
  txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  if verboseflag:
    print(txt.decode("UTF-8"))
