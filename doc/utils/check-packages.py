#!/usr/bin/env python3

import os, re, sys
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser(prog='check-packages.py',
                        description="Check package table completeness")

parser.add_argument("-v", "--verbose",
                    action='store_true',
                    help="Enable verbose output")

parser.add_argument("-d", "--doc",
                    help="Path to LAMMPS documentation sources")
parser.add_argument("-s", "--src",
                    help="Path to LAMMPS sources")

args = parser.parse_args()
verbose = args.verbose
src_dir = args.src
doc_dir = args.doc

LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

if not src_dir:
    src_dir = os.path.join(LAMMPS_DIR , 'src')

if not doc_dir:
    doc_dir = os.path.join(LAMMPS_DIR, 'doc', 'src')

if not src_dir or not doc_dir:
    parser.print_help()
    sys.exit(1)

if not os.path.isdir(src_dir):
    sys.exit(f"LAMMPS source path {src_dir} does not exist")

if not os.path.isdir(doc_dir):
    sys.exit(f"LAMMPS documentation source path {doc_dir} does not exist")

pkgdirs = glob(os.path.join(src_dir, '[A-Z][A-Z]*'))
dirs = re.compile(".*/([0-9A-Z-]+)$")
user = re.compile(".*")

stdpkg = []
usrpkg = []

# find package names and add to standard and user package lists.
# anything starting with at least two upper case characters, is a
# folder, and is not called 'MAKE' is a package 

for d in pkgdirs:
    pkg = dirs.match(d).group(1)
    if not os.path.isdir(os.path.join(src_dir, pkg)): continue
    if pkg in ['DEPEND','MAKE','STUBS']: continue
    if user.match(pkg):
        usrpkg.append(pkg)
    else:
        stdpkg.append(pkg)

print(f"Found {len(stdpkg)} standard and {len(usrpkg)} user packages")

counter = 0

with open(os.path.join(doc_dir, 'Packages_standard.rst')) as fp:
    text = fp.read()

matches = set(re.findall(':ref:`([A-Z0-9-]+) <[A-Z0-9-]+>`', text, re.MULTILINE))
for p in stdpkg:
  if not p in matches:
    counter += 1
    print(f"Standard package {p} missing in Packages_standard.rst")

with open(os.path.join(doc_dir, 'Packages_user.rst')) as fp:
    text = fp.read()

matches = set(re.findall(':ref:`([A-Z0-9-]+) <[A-Z0-9-]+>`', text, re.MULTILINE))
for p in usrpkg:
  if not p in matches:
    counter += 1
    print(f"User package {p} missing in Packages_user.rst")

with open(os.path.join(doc_dir, 'Packages_details.rst')) as fp:
    text = fp.read()

matches = set(re.findall(':ref:`([A-Z0-9]+) <PKG-\\1>`', text, re.MULTILINE))
for p in stdpkg:
    if not p in matches:
        counter += 1
        print(f"Standard package {p} missing in Packages_details.rst")

matches = set(re.findall(':ref:`([A-Z0-9]+) <PKG-\\1>`', text, re.MULTILINE))
for p in usrpkg:
    if not p in matches:
        counter +=1 
        print(f"User package {p} missing in Packages_details.rst")

if counter:
    print(f"Found {counter} issue(s) with package lists")

