#!/usr/bin/env python3

from __future__ import print_function
from glob import glob
from argparse import ArgumentParser
import os, re, sys

parser = ArgumentParser(prog='check-packages.py',
                        description="Check package table completeness")

parser.add_argument("-v", "--verbose",
                    action='store_const',
                    const=True, default=False,
                    help="Enable verbose output")

parser.add_argument("-d", "--doc",
                    help="Path to LAMMPS documentation sources")
parser.add_argument("-s", "--src",
                    help="Path to LAMMPS sources")

args = parser.parse_args()
verbose = args.verbose
src = args.src
doc = args.doc

if not args.src or not args.doc:
  parser.print_help()
  sys.exit(1)

if not os.path.isdir(src):
    sys.exit("LAMMPS source path %s does not exist" % src)

if not os.path.isdir(doc):
    sys.exit("LAMMPS documentation source path %s does not exist" % doc)

pkgdirs = glob(os.path.join(src, '[A-Z][A-Z]*'))
dirs = re.compile(".*/([0-9A-Z-]+)$")
user = re.compile("USER-.*")

stdpkg = []
usrpkg = []

# find package names and add to standard and user package lists.
# anything starting with at least two upper case characters, is a
# folder, and is not called 'MAKE' is a package 

for d in pkgdirs:
  pkg = dirs.match(d).group(1)
  if not os.path.isdir(os.path.join(src,pkg)): continue
  files = glob(os.path.join(src,pkg,"*.cpp"))
  if len(files) == 0: continue
  if pkg in ['DEPEND','MAKE','STUBS']: continue
  if user.match(pkg):
    usrpkg.append(pkg)
  else:
    stdpkg.append(pkg)

print("Found %d standard and %d user packages" % (len(stdpkg),len(usrpkg)))

counter = 0
fp = open(os.path.join(doc,'Packages_standard.rst'))
text = fp.read()
fp.close()
matches = re.findall(':ref:`([A-Z0-9-]+) <[A-Z0-9-]+>`',text,re.MULTILINE)
for p in stdpkg:
  if not p in matches:
    ++counter
    print("Standard package %s missing in Packages_standard.rst" % p)

fp = open(os.path.join(doc,'Packages_user.rst'))
text = fp.read()
fp.close()
matches = re.findall(':ref:`([A-Z0-9-]+) <[A-Z0-9-]+>`',text,re.MULTILINE)
for p in usrpkg:
  if not p in matches:
    ++counter
    print("User package %s missing in Packages_user.rst" % p)

fp = open(os.path.join(doc,'Packages_details.rst'))
text = fp.read()
fp.close()
matches = re.findall(':ref:`([A-Z0-9]+) <PKG-\\1>`',text,re.MULTILINE)
for p in stdpkg:
  if not p in matches:
    ++counter
    print("Standard package %s missing in Packages_details.rst"
          % p)
matches = re.findall(':ref:`(USER-[A-Z0-9]+) <PKG-\\1>`',text,re.MULTILINE)
for p in usrpkg:
  if not p in matches:
    ++counter
    print("User package %s missing in Packages_details.rst"
          % p)

if counter:
    print("Found %d issue(s) with package lists" % counter)

