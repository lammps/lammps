#!/usr/bin/env python3

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

def check_style(file,dir,pattern,list,name,suffix=False,skip=()):
    f = os.path.join(dir, file)
    fp = open(f)
    text = fp.read()
    fp.close()
    matches = re.findall(pattern,text,re.MULTILINE)
    counter = 0
    for c in list.keys():
        # known undocumented aliases we need to skip
        if c in skip: continue
        s = c
        if suffix: s = add_suffix(list,c)
        if not s in matches:
            if not list[c]['removed']:
                print("%s style entry %s" % (name,s),
                      "is missing or incomplete in %s" % file)
                counter += 1
    return counter

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


print("""Parsed style names w/o suffixes from C++ tree in %s:
   Angle styles:     %3d    Atom styles:      %3d
   Body styles:      %3d    Bond styles:      %3d
   Command styles:   %3d    Compute styles:   %3d
   Dihedral styles:  %3d    Dump styles:      %3d
   Fix styles:       %3d    Improper styles:  %3d
   Integrate styles: %3d    Kspace styles:    %3d
   Minimize styles:  %3d    Pair styles:      %3d
   Reader styles:    %3d    Region styles:    %3d""" \
      % (src, len(angle), len(atom), len(body), len(bond),   \
       len(command), len(compute), len(dihedral), len(dump), \
       len(fix), len(improper), len(integrate), len(kspace), \
       len(minimize), len(pair), len(reader), len(region)))


counter = 0

counter += check_style('Commands_all.rst',doc,":doc:`(.+) <.+>`",
                       command,'Command',suffix=False)
counter += check_style('Commands_compute.rst',doc,":doc:`(.+) <compute.+>`",
                       compute,'Compute',suffix=True)
counter += check_style('compute.rst',doc,":doc:`(.+) <compute.+>` -",
                       compute,'Compute',suffix=False)
counter += check_style('Commands_fix.rst',doc,":doc:`(.+) <fix.+>`",
                       fix,'Fix',skip=('python'),suffix=True)
counter += check_style('fix.rst',doc,":doc:`(.+) <fix.+>` -",
                       fix,'Fix',skip=('python'),suffix=False)
counter += check_style('Commands_pair.rst',doc,":doc:`(.+) <pair.+>`",
                       pair,'Pair',skip=('meam','lj/sf'),suffix=True)
counter += check_style('pair_style.rst',doc,":doc:`(.+) <pair.+>` -",
                       pair,'Pair',skip=('meam','lj/sf'),suffix=False)
counter += check_style('Commands_bond.rst',doc,":doc:`(.+) <bond.+>`",
                       bond,'Bond',suffix=True)
counter += check_style('bond_style.rst',doc,":doc:`(.+) <bond.+>` -",
                       bond,'Bond',suffix=False)
counter += check_style('Commands_bond.rst',doc,":doc:`(.+) <angle.+>`",
                       angle,'Angle',suffix=True)
counter += check_style('angle_style.rst',doc,":doc:`(.+) <angle.+>` -",
                       angle,'Angle',suffix=False)
counter += check_style('Commands_bond.rst',doc,":doc:`(.+) <dihedral.+>`",
                       dihedral,'Dihedral',suffix=True)
counter += check_style('dihedral_style.rst',doc,":doc:`(.+) <dihedral.+>` -",
                       dihedral,'Dihedral',suffix=False)
counter += check_style('Commands_bond.rst',doc,":doc:`(.+) <improper.+>`",
                       improper,'Improper',suffix=True)
counter += check_style('improper_style.rst',doc,":doc:`(.+) <improper.+>` -",
                       improper,'Improper',suffix=False)
counter += check_style('Commands_kspace.rst',doc,":doc:`(.+) <kspace_style>`",
                       kspace,'KSpace',suffix=True)

if counter:
    print("Found %d issue(s) with style lists" % counter)

