#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the ScaFaCoS library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess,shutil
sys.path.append('..')
from install_helpers import fullpath,geturl,get_cpus,checkmd5sum
from argparse import ArgumentParser

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

help = """
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
if args.build == False and not args.path:
  parser.print_help()
  sys.exit(help)

buildflag = args.build
pathflag = args.path != None
version = args.version

homepath = fullpath(".")
scafacospath = "%s/scafacos-%s" % (homepath,version)

if (pathflag):
  scafacospath = args.path
  if not os.path.isdir("%s/include" % scafacospath):
    sys.exit("ScaFaCoS include path for %s does not exist" % scafacospath)
  if (not os.path.isdir("%s/lib64" % scafacospath)) \
     and (not os.path.isdir("%s/lib" % scafacospath)):
    sys.exit("ScaFaCoS lib path for %s does not exist" % scafacospath)
  scafacospath = fullpath(scafacospath)

# download and unpack ScaFaCoS tarball

if buildflag:
  print("Downloading ScaFaCoS ...")
  geturl(url,"%s/scafacos-%s.tar.gz" % (homepath,version))

  # verify downloaded archive integrity via md5 checksum, if known.
  if version in checksums:
    if not checkmd5sum(checksums[version],'%s/scafacos-%s.tar.gz' % (homepath,version)):
      sys.exit("Checksum for ScaFaCoS library does not match")

  print("Unpacking ScaFaCoS tarball ...")
  if os.path.exists(scafacospath):
    shutil.rmtree(scafacospath)
  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,scafacospath)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  os.remove("%s.tar.gz" % scafacospath)

  # build ScaFaCoS
  print("Building ScaFaCoS ...")
  n_cpu = get_cpus()
  cmd = 'cd "%s"; ./configure --prefix="%s/build" --disable-doc --enable-fcs-solvers=fmm,p2nfft,direct,ewald,p3m --with-internal-fftw --with-internal-pfft --with-internal-pnfft CC=mpicc FC=mpif90 CXX=mpicxx F77=; make -j%d; make install' % (scafacospath,homepath,n_cpu)
  try:
    txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
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
  cmd = 'ln -s "%s/build/include" includelink' % homepath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s/build/lib" liblink' % homepath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
else:
  cmd = 'ln -s "%s/include" includelink' % scafacospath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  if os.path.isdir("%s/lib64" % scafacospath):
    cmd = 'ln -s "%s/lib64" liblink' % scafacospath
  else:
    cmd = 'ln -s "%s/lib" liblink' % scafacospath
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
