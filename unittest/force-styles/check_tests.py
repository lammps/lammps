#!/usr/bin/env python3

from __future__ import print_function
from glob import glob
from argparse import ArgumentParser
import os, re, sys

parser = ArgumentParser(prog='check_tests.py',
                        description="Check force tests for completeness")

parser.add_argument("-v", "--verbose",
                    action='store_const',
                    const=True, default=False,
                    help="Enable verbose output")

parser.add_argument("-t", "--tests",
                    help="Path to LAMMPS test YAML format input files")
parser.add_argument("-s", "--src",
                    help="Path to LAMMPS sources")

args = parser.parse_args()
verbose = args.verbose
src = args.src
tests = args.tests

if not src:
  src = os.path.join('.','..','..','src')

if not tests:
  tests = os.path.join(src,'..','unittest','force-styles','tests')

try:
  src = os.path.abspath(os.path.expanduser(src))
  tests = os.path.abspath(os.path.expanduser(tests))
except:
  parser.print_help()
  sys.exit(1)

if not os.path.isdir(src):
    sys.exit("LAMMPS source path %s does not exist" % src)

if not os.path.isdir(tests):
    sys.exit("LAMMPS test inputs path %s does not exist" % tests)

headers = glob(os.path.join(src, '*', '*.h'))
headers += glob(os.path.join(src, '*.h'))

angle = {}
bond = {}
compute = {}
dihedral = {}
dump = {}
fix = {}
improper = {}
kspace = {}
pair = {}
total = 0

upper = re.compile("[A-Z]+")
gpu = re.compile("(.+)/gpu$")
intel = re.compile("(.+)/intel$")
kokkos = re.compile("(.+)/kk$")
kokkos_skip = re.compile("(.+)/kk/(host|device)$")
omp = re.compile("(.+)/omp$")
opt = re.compile("(.+)/opt$")
removed = re.compile("(.*)Deprecated$")

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

print("Parsing style names from C++ tree in: ",src)

for h in headers:
    if verbose: print("Checking ", h)
    fp = open(h)
    text = fp.read()
    fp.close()
    matches = re.findall("(.+)Style\((.+),(.+)\)",text,re.MULTILINE)
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
        elif m[0] == 'Bond':
            register_style(bond,style,info)
        elif m[0] == 'Dihedral':
            register_style(dihedral,style,info)
        elif m[0] == 'Improper':
            register_style(improper,style,info)
        elif m[0] == 'KSpace':
            register_style(kspace,style,info)
        elif m[0] == 'Pair':
            register_style(pair,style,info)


counter = 0

def check_tests(name,list,yaml,search,skip=()):
    num = 0
    yaml_files = glob(os.path.join(tests, yaml))
    styles = []
    missing = []
    for y in yaml_files:
        if verbose: print("Checking: ",y)
        fp = open(y)
        text = fp.read()
        fp.close()
        matches = re.findall(search,text,re.MULTILINE)
        for m in matches:
            if m[1] == 'hybrid' or m[1] == 'hybrid/overlay':
                for s in m[0].split():
                    if s in list.keys():
                        styles.append(s)
            else:
                styles.append(m[1])
    for s in list.keys():
        # known undocumented aliases we need to skip
        if s in skip: continue
        if not s in styles:
            if not list[s]['removed']:
                if verbose: print("No test for %s style %s" % (name,s))
                num += 1
                missing.append(s)
    total = len(list)
    print("\nTests available for %s styles: %d of %d"
          % (name,total - num, total))
    print("No tests for: ", missing)
    return num

counter += check_tests('pair',pair,'*-pair-*.yaml',
                       '.*pair_style:\s*((\S+).*)?',skip=('meam','lj/sf'))
counter += check_tests('bond',bond,'bond-*.yaml',
                       '.*bond_style:\s*((\S+).*)?')
counter += check_tests('angle',angle,'angle-*.yaml',
                       '.*angle_style:\s*((\S+).*)?')
counter += check_tests('dihedral',dihedral,'dihedral-*.yaml',
                       '.*dihedral_style:\s*((\S+).*)?')
counter += check_tests('improper',improper,'improper-*.yaml',
                       '.*improper_style:\s*((\S+).*)?')
counter += check_tests('kspace',kspace,'kspace-*.yaml',
                       '.*kspace_style\s*((\S+).*)?')

total = len(pair)+len(bond)+len(angle)+len(dihedral)+len(improper)+len(kspace)
print("\nTotal tests missing: %d of %d" % (counter, total))
