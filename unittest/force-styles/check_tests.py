#!/usr/bin/env python3
import os, re, sys
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser(prog='check_tests.py',
                        description="Check force tests for completeness")

parser.add_argument("-v", "--verbose",
                    action='store_true',
                    help="Enable verbose output")

parser.add_argument("-t", "--tests",
                    help="Path to LAMMPS test YAML format input files")
parser.add_argument("-s", "--src",
                    help="Path to LAMMPS sources")

args = parser.parse_args()
verbose = args.verbose
src_dir = args.src
tests_dir = args.tests

LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

if not src_dir:
    src_dir = os.path.join(LAMMPS_DIR , 'src')

if not tests_dir:
    tests_dir = os.path.join(LAMMPS_DIR, 'unittest', 'force-styles', 'tests')

try:
    src_dir = os.path.abspath(os.path.expanduser(src_dir))
    tests_dir = os.path.abspath(os.path.expanduser(tests_dir))
except:                                 # lgtm [py/catch-base-exception]
    parser.print_help()
    sys.exit(1)

if not os.path.isdir(src_dir):
    sys.exit(f"LAMMPS source path {src_dir} does not exist")

if not os.path.isdir(tests_dir):
    sys.exit(f"LAMMPS test inputs path {tests_dir} does not exist")

headers = glob(os.path.join(src_dir, '*', '*.h'))
headers += glob(os.path.join(src_dir, '*.h'))

angle = {}
bond = {}
compute = {}
dihedral = {}
dump = {}
fix = {}
improper = {}
kspace = {}
pair = {}

style_pattern = re.compile("(.+)Style\((.+),(.+)\)")
upper = re.compile("[A-Z]+")
gpu = re.compile("(.+)/gpu$")
intel = re.compile("(.+)/intel$")
kokkos = re.compile("(.+)/kk$")
kokkos_skip = re.compile("(.+)/kk/(host|device)$")
omp = re.compile("(.+)/omp$")
opt = re.compile("(.+)/opt$")
removed = re.compile("(.*)Deprecated$")

def register_style(styles, name, info):
    if name in styles:
        styles[name]['gpu'] += info['gpu']
        styles[name]['intel'] += info['intel']
        styles[name]['kokkos'] += info['kokkos']
        styles[name]['omp'] += info['omp']
        styles[name]['opt'] += info['opt']
        styles[name]['removed'] += info['removed']
    else:
        styles[name] = info

print("Parsing style names from C++ tree in: ", src_dir)

for header in headers:
    if verbose: print("Checking ", header)
    with open(header) as f:
        for line in f:
            matches = style_pattern.findall(line)

            for m in matches:
                # skip over internal styles w/o explicit documentation
                style = m[1]
                if upper.match(style):
                    continue
                if style in ['reax/c', 'reax/c/omp', 'reax/c/kk',
                        'reax/c/kk/device', 'reax/c/kk/host',
                        'reax/c/species', 'reax/c/bonds',
                        'reax/c/species/kk', 'reax/c/bonds/kk', 'meam/c']:
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
                elif m[0] == 'Fix':
                    register_style(fix,style,info)



def check_tests(name,styles,yaml,search,skip=()):
    yaml_files = glob(os.path.join(tests_dir, yaml))
    tests = set()
    missing = set()
    search_pattern = re.compile(search)
    for y in yaml_files:
        if verbose: print("Checking: ",y)
        with open(y) as f:
            for line in f:
                matches = search_pattern.findall(line)
                for m in matches:
                    if m[1] == 'hybrid' or m[1] == 'hybrid/overlay':
                        for s in m[0].split():
                            if s in styles:
                                tests.add(s)
                    else:
                        tests.add(m[1])
    for s in styles:
        # known undocumented aliases we need to skip
        if s in skip: continue
        if not s in tests:
            if not styles[s]['removed']:
                if verbose: print("No test for %s style %s" % (name,s))
                missing.add(s)

    num_missing = len(missing)
    total = len(styles)
    num_tests = total - num_missing
    print(f"\nTests available for {name} styles: {num_tests} of {total}")
    print("No tests for: ", sorted(missing))
    return num_missing

counter = 0
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
counter += check_tests('fix',fix,'fix-*.yaml',
                       '  fix\s+((\S+)\s*)?')

total = len(pair)+len(bond)+len(angle)+len(dihedral)+len(improper)+len(kspace)+len(fix)
print(f"\nTotal tests missing: {counter} of {total}")
