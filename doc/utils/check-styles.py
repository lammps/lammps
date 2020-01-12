#!/usr/bin/env python

from __future__ import print_function
from glob import glob
from argparse import ArgumentParser
import os, re, sys

parser = ArgumentParser(prog='check-styles.py',
                        description="Check style table completeness")

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

headers = glob(os.path.join(src, '*', '*.h'))
headers += glob(os.path.join(src, '*.h'))

angle = {}
atom = {}
body = {}
bond = {}
command = {}
compute = {}
dihedral = {}
dump = {}
fix = {}
improper = {}
integrate = {}
kspace = {}
minimize = {}
pair = {}
reader = {}
region = {}

upper = re.compile("[A-Z]+")
gpu = re.compile("(.+)/gpu$")
intel = re.compile("(.+)/intel$")
kokkos = re.compile("(.+)/kk$")
kokkos_skip = re.compile("(.+)/kk/(host|device)$")
omp = re.compile("(.+)/omp$")
opt = re.compile("(.+)/opt$")
removed = re.compile("(.+)Deprecated$")

def register_style(list,style,info):
    if style in list.keys():
        list[style]['gpu'] += info['gpu']
        list[style]['intel'] += info['intel']
        list[style]['kokkos'] += info['kokkos']
        list[style]['omp'] += info['omp']
        list[style]['opt'] += info['opt']
        list[style]['removed'] += info['removed']
    else:
        list[style] = info

def add_suffix(list,style):
    suffix = ""
    if list[style]['gpu']:
        suffix += 'g'
    if list[style]['intel']:
        suffix += 'i'
    if list[style]['kokkos']:
        suffix += 'k'
    if list[style]['omp']:
        suffix += 'o'
    if list[style]['opt']:
        suffix += 't'
    if suffix:
        return style + ' (' + suffix + ')'
    else:
        return style
        
for h in headers:
    if verbose: print("Checking ", h)
    fp = open(h)
    text = fp.read()
    fp.close()
    matches = re.findall("(.+)Style\((.+),(.+)\)",text,re.MULTILINE)
    for m in matches:

        # skip over internal styles w/o explicit documentation
        style = m[1]
        if upper.match(style):
            continue

        # detect, process, and flag suffix styles:
        info = { 'kokkos':  0, 'gpu':     0, 'intel':   0, \
                 'omp':     0, 'opt':     0, 'removed': 0 }
        suffix = kokkos_skip.match(style)
        if suffix:
            continue
        suffix = gpu.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['gpu'] = 1
        suffix = intel.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['intel'] = 1
        suffix = kokkos.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['kokkos'] = 1
        suffix = omp.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['omp'] = 1
        suffix = opt.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['opt'] = 1
        deprecated = removed.match(m[2])
        if deprecated:
            info['removed'] = 1

        # register style and suffix flags
        if m[0] == 'Angle':
            register_style(angle,style,info)
        elif m[0] == 'Atom':
            register_style(atom,style,info)
        elif m[0] == 'Body':
            register_style(body,style,info)
        elif m[0] == 'Bond':
            register_style(bond,style,info)
        elif m[0] == 'Command':
            register_style(command,style,info)
        elif m[0] == 'Compute':
            register_style(compute,style,info)
        elif m[0] == 'Dihedral':
            register_style(dihedral,style,info)
        elif m[0] == 'Dump':
            register_style(dump,style,info)
        elif m[0] == 'Fix':
            register_style(fix,style,info)
        elif m[0] == 'Improper':
            register_style(improper,style,info)
        elif m[0] == 'Integrate':
            register_style(integrate,style,info)
        elif m[0] == 'KSpace':
            register_style(kspace,style,info)
        elif m[0] == 'Minimize':
            register_style(minimize,style,info)
        elif m[0] == 'Pair':
            register_style(pair,style,info)
        elif m[0] == 'Reader':
            register_style(reader,style,info)
        elif m[0] == 'Region':
            register_style(region,style,info)
        else:
            print("Skipping over: ",m)


if verbose:
    print("""
Parsed styles from %s:
  Angle styles:     %3d
  Atom styles:      %3d
  Body styles:      %3d
  Bond styles:      %3d
  Command styles:   %3d
  Compute styles:   %3d
  Dihedral styles:  %3d
  Dump styles:      %3d
  Fix styles:       %3d
  Improper styles:  %3d
  Integrate styles: %3d
  Kspace styles:    %3d
  Minimize styles:  %3d
  Pair styles:      %3d
  Reader styles:    %3d
  Region styles:    %3d
""" % (src, len(angle), len(atom), len(body), len(bond),     \
       len(command), len(compute), len(dihedral), len(dump), \
       len(fix), len(improper), len(integrate), len(kspace), \
       len(minimize), len(pair), len(reader), len(region)))


# check main commands lists
f = os.path.join(doc, 'Commands_all.rst')
fp = open(f)
text = fp.read()
fp.close()
matches = re.findall(":doc:`(.+) <.+>`",text,re.MULTILINE)
for c in command.keys():
    if not c in matches:
        if not command[c]['removed']:
            print("Command %s is missing in Commands_all.rst" % c)

f = os.path.join(doc, 'Commands_compute.rst')
fp = open(f)
text = fp.read()
fp.close()
matches = re.findall(":doc:`(.+) <compute.+>`",text,re.MULTILINE)
for c in compute.keys():
    if not add_suffix(compute,c) in matches:
        if not compute[c]['removed']:
            print("Compute style entry %s is missing or" % c,
                  "incomplete in Commands_compute.rst")

f = os.path.join(doc, 'Commands_fix.rst')
fp = open(f)
text = fp.read()
fp.close()
matches = re.findall(":doc:`(.+) <fix.+>`",text,re.MULTILINE)
for c in fix.keys():
    # known undocumented aliases we need to skip
    if c in ('python'): continue
    if not add_suffix(fix,c) in matches:
        if not fix[c]['removed']:
            print("Fix style entry %s is missing or" % c,
                  "incomplete in Commands_fix.rst")

f = os.path.join(doc, 'Commands_pair.rst')
fp = open(f)
text = fp.read()
fp.close()
matches = re.findall(":doc:`(.+) <pair.+>`",text,re.MULTILINE)
for c in pair.keys():
    # known undocumented aliases we need to skip
    if c in ('meam','lj/sf'): continue
    if not add_suffix(pair,c) in matches:
        if not pair[c]['removed']:
            print("Pair style entry %s is missing or" % c,
                  "incomplete in Commands_pair.rst")

f = os.path.join(doc, 'Commands_bond.rst')
fp = open(f)
text = fp.read()
fp.close()
matches = re.findall(":doc:`(.+) <bond.+>`",text,re.MULTILINE)
for c in bond.keys():
    if not add_suffix(bond,c) in matches:
        if not bond[c]['removed']:
            print("Bond style entry %s is missing or" % c,
                  "incomplete in Commands_bond.rst")

matches = re.findall(":doc:`(.+) <angle.+>`",text,re.MULTILINE)
for c in angle.keys():
    if not add_suffix(angle,c) in matches:
        if not angle[c]['removed']:
            print("Angle style entry %s is missing or" % c,
                  "incomplete in Commands_bond.rst")

matches = re.findall(":doc:`(.+) <dihedral.+>`",text,re.MULTILINE)
for c in dihedral.keys():
    if not add_suffix(dihedral,c) in matches:
        if not dihedral[c]['removed']:
            print("Dihedral style entry %s is missing or" % c,
                  "incomplete in Commands_bond.rst")

matches = re.findall(":doc:`(.+) <improper.+>`",text,re.MULTILINE)
for c in improper.keys():
    if not add_suffix(improper,c) in matches:
        if not improper[c]['removed']:
            print("Improper style entry %s is missing or" % c,
                  "incomplete in Commands_bond.rst")

f = os.path.join(doc, 'Commands_kspace.rst')
fp = open(f)
text = fp.read()
fp.close()
matches = re.findall(":doc:`(.+) <kspace_style>`",text,re.MULTILINE)
for c in kspace.keys():
    if not add_suffix(kspace,c) in matches:
        if not kspace[c]['removed']:
            print("KSpace style entry %s is missing or" % c,
                  "incomplete in Commands_kspace.rst")
            print(kspace[c])
