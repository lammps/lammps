#!/usr/bin/env python

"""
Install.py tool to download, compile, and setup the pace library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function

import shutil
import subprocess
import sys
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, checkmd5sum

# settings

thisdir = fullpath('.')
version ='v.2023.01.3'

# known checksums for different PACE versions. used to validate the download.
checksums = { \
    'v.2023.01.3': 'f418d32b60e531063ac4285bf702b468'
}

parser = ArgumentParser(prog='Install.py', description="LAMMPS library build wrapper script")

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
parser.add_argument("-l", "--local", default=None,
                    help="use local version of PACE library build")

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build:
    parser.print_help()
    sys.exit(HELP)

buildflag = args.build

verboseflag = args.verbose
version = args.version
local = args.local

archive_extension = "tar.gz"
url = "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags/%s.%s" % (version, archive_extension)
unarchived_folder_name = "lammps-user-pace-%s" % (version)

# download PACE tarball, unpack, build PACE
if buildflag:
    if not local:
        # download entire tarball
        print("Downloading pace tarball ...")
        archive_filename = "%s.%s" % (version, archive_extension)
        download_filename = "%s/%s" % (thisdir, archive_filename)
        print("Downloading from ", url, " to ", download_filename, end=" ")
        geturl(url, download_filename)
        print(" done")

        # verify downloaded archive integrity via md5 checksum, if known.
        if version in checksums:
            if not checkmd5sum(checksums[version], archive_filename):
                sys.exit("Checksum for pace library does not match")

        print("Unpacking pace tarball ...")
        src_folder = thisdir + "/src"
        cmd = 'cd "%s"; rm -rf "%s"; tar -xvf %s; mv %s %s' % (
            thisdir, src_folder, archive_filename, unarchived_folder_name, src_folder)
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    else:
        # copy from local version of library PACE
        print("Copy pace from ", local)
        src_folder = thisdir + "/src"
        shutil.copytree(local, src_folder,
                        # ignore=lambda (s1,s2): ('.git' in s1 or '.git' in s2),
                        dirs_exist_ok=True)


    # build
    print("Building libpace ...")
    cmd = 'make lib -j2'
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    if verboseflag:
        print(txt.decode("UTF-8"))

    #   remove source files

    print("Removing pace build files and archive ...")
    cmd = 'make clean-build'
    if not local:
        cmd = ('rm %s;' % (download_filename))+cmd
    subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

