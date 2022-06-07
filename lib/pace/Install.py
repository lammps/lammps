# TODO#!/usr/bin/env python

"""
Install.py tool to download, compile, and setup the pace library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, subprocess
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

# settings

thisdir = fullpath('.')
version = 'v.2021.10.25.fix2'

# known checksums for different PACE versions. used to validate the download.
checksums = { \
        'v.2021.2.3.upd2' : '8fd1162724d349b930e474927197f20d',
        'v.2021.4.9'      : '4db54962fbd6adcf8c18d46e1798ceb5',
        'v.2021.9.28'     : 'f98363bb98adc7295ea63974738c2a1b',
        'v.2021.10.25'    : 'a2ac3315c41a1a4a5c912bcb1bc9c5cc',
        'v.2021.10.25.fix': 'e0572de57039d4afedefb25707b6ceae',
        'v.2021.10.25.fix2': '32394d799bc282bb57696c78c456e64f'
        }


parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")


# help message

HELP = """
Syntax from src dir: make lib-pace args="-b"
                 or: make lib-pace args="-b -v version"
Syntax from lib dir: python Install.py -b
                 or: python Install.py -b -v version

Examples:

make lib-pace args="-b" # install default version of PACE lib
make lib-pace args="-b -v version" # install specified version of PACE lib


"""

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build base PACE  library")
parser.add_argument("-v", "--version", default=version, choices=checksums.keys(),
                    help="set version of PACE library to download and build (default: %s)" % version)
parser.add_argument("-vv", "--verbose", action="store_true",
                    help="be more verbose about is happening while this script runs")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build

verboseflag = args.verbose
version = args.version


archive_extension = "tar.gz"
url = "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/%s.%s" % (version, archive_extension)
unarchived_folder_name = "lammps-user-pace-%s"%(version)

# download PACE tarball, unpack, build PACE
if buildflag:

  # download entire tarball

  print("Downloading pace tarball ...")
  archive_filename = "%s.%s" % (version, archive_extension)
  download_filename = "%s/%s" % (thisdir, archive_filename)
  print("Downloading from ",url," to ",download_filename, end=" ")
  geturl(url, download_filename)
  print(" done")

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], archive_filename):
      sys.exit("Checksum for pace library does not match")

  print("Unpacking pace tarball ...")
  src_folder = thisdir+"/src"
  cmd = 'cd "%s"; rm -rf "%s"; tar -xvf %s; mv %s %s' % (thisdir, src_folder, archive_filename, unarchived_folder_name, src_folder)
  subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

  # build
  print("Building libpace ...")
  cmd = 'make lib -j2'
  txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
  if verboseflag:
    print(txt.decode("UTF-8"))

#   remove source files

  print("Removing pace build files and archive ...")
  cmd = 'rm %s; make clean-build' % (download_filename)
  subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
