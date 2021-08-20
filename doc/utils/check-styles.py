#!/usr/bin/env python3

import os, re, sys
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser(prog='check-styles.py',
                        description="Check style table completeness")

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

headers = glob(os.path.join(src_dir, '*', '*.h'))
headers += glob(os.path.join(src_dir, '*.h'))

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
total = 0

index_pattern = re.compile(r"^.. index:: (compute|fix|pair_style|angle_style|bond_style|dihedral_style|improper_style|kspace_style)\s+([a-zA-Z0-9/_]+)$")
style_pattern = re.compile(r"(.+)Style\((.+),(.+)\)")
upper = re.compile("[A-Z]+")
gpu = re.compile("(.+)/gpu$")
intel = re.compile("(.+)/intel$")
kokkos = re.compile("(.+)/kk$")
kokkos_skip = re.compile("(.+)/kk/(host|device)$")
omp = re.compile("(.+)/omp$")
opt = re.compile("(.+)/opt$")
removed = re.compile("(.*)Deprecated$")

def load_index_entries_in_file(path):
    entries = []
    with open(path, 'r') as reader:
        for line in reader:
            m = index_pattern.match(line)
            if m:
                command_type = m.group(1)
                style = m.group(2)
                entries.append((command_type, style))
    return entries

def load_index_entries():
    index = {'compute': set(), 'fix': set(), 'pair_style': set(), 'angle_style': set(),
             'bond_style': set(), 'dihedral_style': set(), 'improper_style': set(), 'kspace_style': set()}
    rst_files = glob(os.path.join(doc_dir, '*.rst'))
    for f in rst_files:
        for command_type, style in load_index_entries_in_file(f):
            index[command_type].add(style)
    return index

def register_style(styles, name, info):
    if name in styles:
        for key, value in info.items():
            styles[name][key] += value
    else:
        styles[name] = info

def add_suffix(styles, name):
    suffix = ""
    if styles[name]['gpu']:
        suffix += 'g'
    if styles[name]['intel']:
        suffix += 'i'
    if styles[name]['kokkos']:
        suffix += 'k'
    if styles[name]['omp']:
        suffix += 'o'
    if styles[name]['opt']:
        suffix += 't'
    if suffix:
        return f"{name} ({suffix})"
    else:
        return name

def check_style(filename, dirname, pattern, styles, name, suffix=False, skip=set()):
    with open(os.path.join(dirname, filename)) as f:
        text = f.read()

    matches = re.findall(pattern, text, re.MULTILINE)

    counter = 0
    for c in styles:
        # known undocumented aliases we need to skip
        if c in skip: continue
        s = c
        if suffix: s = add_suffix(styles, c)
        if not s in matches:
            if not styles[c]['removed']:
                print(f"{name} style entry {s} is missing or incomplete in {filename}")
                counter += 1
    return counter

def check_style_index(name, styles, index, skip=[]):
    counter = 0
    for style in styles:
        if style not in index and not styles[style]['removed'] and style not in skip:
            print(f"{name} index entry {style} is missing")
            counter += 1

        for suffix in styles[style]:
            if suffix == 'removed': continue
            if suffix == 'kokkos':
                suffix_style = f"{style}/kk"
            else:
                suffix_style = f"{style}/{suffix}"
            if styles[style][suffix] and suffix_style not in index and style not in skip:
                print(f"{name} index entry {suffix_style} is missing")
                counter += 1
    return counter

for header in headers:
    if verbose: print("Checking ", header)
    with open(header) as f:
        for line in f:
            matches = style_pattern.findall(line)
            for m in matches:
                # skip over internal styles w/o explicit documentation
                style = m[1]
                total += 1
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
   Reader styles:    %3d    Region styles:    %3d
----------------------------------------------------
Total number of styles (including suffixes): %d""" \
      % (src_dir, len(angle), len(atom), len(body), len(bond),   \
       len(command), len(compute), len(dihedral), len(dump), \
       len(fix), len(improper), len(integrate), len(kspace), \
       len(minimize), len(pair), len(reader), len(region), total))

index = load_index_entries()

total_index = 0
for command_type, entries in index.items():
    total_index += len(entries)

print("Total number of style index entries:", total_index)

skip_fix = ('python', 'NEIGH_HISTORY/omp','qeq/reax','reax/c/bonds','reax/c/species')
skip_pair = ('meam/c','lj/sf','reax/c')

counter = 0

counter += check_style('Commands_all.rst', doc_dir, ":doc:`(.+) <.+>`",command,'Command',suffix=False)
counter += check_style('Commands_compute.rst', doc_dir, ":doc:`(.+) <compute.+>`",compute,'Compute',suffix=True)
counter += check_style('compute.rst', doc_dir, ":doc:`(.+) <compute.+>` -",compute,'Compute',suffix=False)
counter += check_style('Commands_fix.rst', doc_dir, ":doc:`(.+) <fix.+>`",fix,'Fix',skip=skip_fix,suffix=True)
counter += check_style('fix.rst', doc_dir, ":doc:`(.+) <fix.+>` -",fix,'Fix',skip=skip_fix,suffix=False)
counter += check_style('Commands_pair.rst', doc_dir, ":doc:`(.+) <pair.+>`",pair,'Pair',skip=skip_pair,suffix=True)
counter += check_style('pair_style.rst', doc_dir, ":doc:`(.+) <pair.+>` -",pair,'Pair',skip=skip_pair,suffix=False)
counter += check_style('Commands_bond.rst', doc_dir, ":doc:`(.+) <bond.+>`",bond,'Bond',suffix=True)
counter += check_style('bond_style.rst', doc_dir, ":doc:`(.+) <bond.+>` -",bond,'Bond',suffix=False)
counter += check_style('Commands_bond.rst', doc_dir, ":doc:`(.+) <angle.+>`",angle,'Angle',suffix=True)
counter += check_style('angle_style.rst', doc_dir, ":doc:`(.+) <angle.+>` -",angle,'Angle',suffix=False)
counter += check_style('Commands_bond.rst', doc_dir, ":doc:`(.+) <dihedral.+>`",dihedral,'Dihedral',suffix=True)
counter += check_style('dihedral_style.rst', doc_dir, ":doc:`(.+) <dihedral.+>` -",dihedral,'Dihedral',suffix=False)
counter += check_style('Commands_bond.rst', doc_dir, ":doc:`(.+) <improper.+>`",improper,'Improper',suffix=True)
counter += check_style('improper_style.rst', doc_dir, ":doc:`(.+) <improper.+>` -",improper,'Improper',suffix=False)
counter += check_style('Commands_kspace.rst', doc_dir, ":doc:`(.+) <kspace_style>`",kspace,'KSpace',suffix=True)

if counter:
    print(f"Found {counter} issue(s) with style lists")

counter = 0

counter += check_style_index("compute", compute, index["compute"])
counter += check_style_index("fix", fix, index["fix"], skip=['python','qeq/reax','reax/c/bonds','reax/c/species'])
counter += check_style_index("angle_style", angle, index["angle_style"])
counter += check_style_index("bond_style", bond, index["bond_style"])
counter += check_style_index("dihedral_style", dihedral, index["dihedral_style"])
counter += check_style_index("improper_style", improper, index["improper_style"])
counter += check_style_index("kspace_style", kspace, index["kspace_style"])
counter += check_style_index("pair_style", pair, index["pair_style"], skip=['meam/c', 'lj/sf','reax/c'])

if counter:
    print(f"Found {counter} issue(s) with style index")
