#!/usr/bin/env python

"""
Install.py tool to download, unpack, build, and link to the ScaFaCoS library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil, tarfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, get_cpus, checkmd5sum

parser = ArgumentParser(prog='Install.py',
                        description="LAMMPS library build wrapper script")

# settings

version = "1.0.1"
url = "https://github.com/scafacos/scafacos/releases/download/v%s/scafacos-%s.tar.gz" % (version, version)

# known checksums for different ScaFaCoS versions. used to validate the download.
checksums = { \
        '1.0.1' : 'bd46d74e3296bd8a444d731bb10c1738' \
        }

# extra help message

HELP = """
Syntax from src dir: make lib-scafacos args="-b"
                 or: make lib-scafacos args="-p /usr/local/scafacos"
Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/scafacos

Example:

make lib-scafacos args="-b"   # download/build in lib/scafacos/scafacos
make lib-scafacos args="-p $HOME/scafacos" # use existing ScaFaCoS installation in $HOME
"""

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the ScaFaCoS library")
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing ScaFaCoS installation")
parser.add_argument("-v", "--version", default=version,
                    help="set version of ScaFaCoS to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
version = args.version

homepath = fullpath(".")
scafacospath = os.path.join(homepath, "scafacos-%s" % version)

if pathflag:
  scafacospath = args.path
  if not os.path.isdir(os.path.join(scafacospath, "include")):
    sys.exit("ScaFaCoS include path for %s does not exist" % scafacospath)
  if (not os.path.isdir(os.path.join(scafacospath, "lib64"))) \
     and (not os.path.isdir(os.path.join(scafacospath, "lib"))):
    sys.exit("ScaFaCoS lib path for %s does not exist" % scafacospath)
  scafacospath = fullpath(scafacospath)

# download and unpack ScaFaCoS tarball

if buildflag:
  print("Downloading ScaFaCoS ...")
  geturl(url, "%s/scafacos-%s.tar.gz" % (homepath, version))

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version], '%s/scafacos-%s.tar.gz' % (homepath, version)):
      sys.exit("Checksum for ScaFaCoS library does not match")

  print("Unpacking ScaFaCoS tarball ...")
  if os.path.exists(scafacospath):
    shutil.rmtree(scafacospath)
  tarname = os.path.join(homepath, "%s.tar.gz" % scafacospath)
  if tarfile.is_tarfile(tarname):
    tgz = tarfile.open(tarname)
    tgz.extractall(path=homepath)
    os.remove(tarname)
  else:
    sys.exit("File %s is not a supported archive" % tarname)

  # build ScaFaCoS
  print("Building ScaFaCoS ...")
  n_cpu = get_cpus()
  cmd = 'cd "%s"; ./configure --prefix="%s" --disable-doc --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m --with-internal-fftw --with-internal-pfft --with-internal-pnfft CC=mpicc FC=mpif90 CXX=mpicxx F77=; make -j%d; make install' % (scafacospath, os.path.join(homepath, 'build'), n_cpu)
  try:
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    print(txt.decode('UTF-8'))
  except subprocess.CalledProcessError as e:
    sys.exit("Make failed with:\n %s" % e.output.decode('UTF-8'))

# create 2 links in lib/scafacos to ScaFaCoS include/lib dirs

print("Creating links to ScaFaCoS include and lib files")
if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")
if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")
if buildflag:
  os.symlink(os.path.join(homepath, 'build', 'include'), 'includelink')
  os.symlink(os.path.join(homepath, 'build', 'lib'), 'liblink')
else:
  os.symlink(os.path.join(scafacospath, 'include'), 'includelink')
  if os.path.isdir(os.path.join(scafacospath, "lib64")):
    os.symlink(os.path.join(scafacospath, 'lib64'), 'liblink')
  else:
    os.symlink(os.path.join(scafacospath, 'lib'), 'liblink')
