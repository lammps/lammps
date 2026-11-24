#!/usr/bin/env python3
"""Check or update the "gpu" and "kokkos" entries on the tags: line of the YAML
files of the force-style tests.

The test programs run every YAML file also with the accelerated variants of
the tested style, as far as they are compiled in: "gpu" when the GPU package
variant of the test runs (-sf gpu), "kokkos" when one of the KOKKOS package
variants runs (-sf kk).  The tags become CTest labels, so that

    ctest --test-dir build -L gpu
    ctest --test-dir build -L kokkos

run only the tests that exercise styles of these packages.  This script
applies the rules of the test programs to the styles in the source tree:

  - pair styles (also the kspace tests): every "pair" prerequisite must have
    a variant with the suffix; hybrid pair styles only need their sub-styles
    (which are prerequisites as well) except with KOKKOS, which also needs
    hybrid/kk (or hybrid/overlay/kk, ...)
  - bond, angle, dihedral, improper, and fix-timestep tests: every
    prerequisite of that kind must have a variant with the suffix
  - minimizer tests: every "minimize" prerequisite needs the suffix
  - output tests (compute-*, fix-output-*): no suffix needed, the KOKKOS
    variants always run, with the plain style where there is no /kk variant.
    These count, too: KOKKOS keeps its data on the GPU as much as possible,
    and the transfers for styles without a /kk variant must work as well
    (the GPU package only accelerates pair styles and ewald/pppm, so this
    does not apply there)
  - only the pair style tests have a GPU package variant
  - a variant listed in skip_tests does not count; entries qualified with a
    precision or "_devicerng" apply only to some builds and are ignored, so
    a file can have a tag although some builds skip all its variants

Without options it lists the files whose tags differ and exits with status 1
if there are any; with -u it rewrites the tags: lines of those files.
"""

import os
import re
import sys
from argparse import ArgumentParser

parser = ArgumentParser(prog='check_accelerator_tags.py',
                        description="Check or update the gpu and kokkos tags of the force-style tests")
parser.add_argument("-u", "--update", action='store_true',
                    help="Rewrite the tags: lines of the YAML files that differ")
parser.add_argument("-v", "--verbose", action='store_true',
                    help="List the tags of every file")
parser.add_argument("-s", "--src", help="Path to the LAMMPS sources")
parser.add_argument("-t", "--tests", help="Path to the YAML files of the force-style tests")
args = parser.parse_args()

TEST_DIR = os.path.dirname(os.path.realpath(__file__))
LAMMPS_DIR = os.path.realpath(os.path.join(TEST_DIR, '..', '..'))
src_dir = os.path.abspath(os.path.expanduser(args.src or os.path.join(LAMMPS_DIR, 'src')))
tests_dir = os.path.abspath(os.path.expanduser(args.tests or os.path.join(TEST_DIR, 'tests')))
for path in (src_dir, tests_dir):
    if not os.path.isdir(path):
        sys.exit(f"Folder {path} does not exist")

# file name prefix of the YAML files -> test program and the kind of prerequisite
# that needs the accelerator suffix (as registered in CMakeLists.txt next to this file)
PAIR = ('test_pair_style', 'pair')
GROUPS = {'mol-pair-': PAIR, 'atomic-pair-': PAIR, 'manybody-pair-': PAIR,
          'ellipsoid-pair-': PAIR, 'body-pair-': PAIR, 'spin-pair-': PAIR,
          'sph-pair-': PAIR, 'meso-pair-': PAIR, 'kspace-': PAIR,
          'bond-': ('test_bond_style', 'bond'), 'angle-': ('test_angle_style', 'angle'),
          'dihedral-': ('test_dihedral_style', 'dihedral'),
          'improper-': ('test_improper_style', 'improper'),
          'fix-timestep-': ('test_fix_timestep', 'fix'),
          'min-': ('test_min_style', 'minimize'),
          'compute-': ('test_output_style', None), 'fix-output-': ('test_output_style', None)}
# the tags this script manages, in the order it writes them
MANAGED = ('gpu', 'kokkos')

# registered styles: category -> set of style names (with suffixes)
style_pattern = re.compile(r'^\s*(\w+)Style\(\s*([\w/.-]+)\s*,', re.M)
styles = {}
for path, dirs, files in os.walk(src_dir):
    for name in files:
        if name.endswith('.h'):
            with open(os.path.join(path, name), errors='replace') as f:
                for kind, style in style_pattern.findall(f.read()):
                    styles.setdefault(kind.lower(), set()).add(style)

# the names of the test cases of each test program
test_pattern = re.compile(r'^TEST\(\s*\w+\s*,\s*(\w+)\s*\)', re.M)
cases = {}
for program, kind in set(GROUPS.values()):
    with open(os.path.join(TEST_DIR, program + '.cpp')) as f:
        cases[program] = set(test_pattern.findall(f.read()))

with open(os.path.join(TEST_DIR, 'CMakeLists.txt')) as f:
    cmake = f.read()
for prefix in GROUPS:
    if '${TEST_INPUT_FOLDER}/' + prefix + '*.yaml' not in cmake:
        print(f"Warning: CMakeLists.txt does not register the files {prefix}*.yaml")


def read_yaml(lines):
    """Return the tags, skip_tests, and prerequisites of a YAML test file"""
    tags, skip, prereq = [], set(), []
    block = None
    for line in lines:
        if line.startswith((' ', '\t')):
            if block == 'prerequisites' and line.split():
                prereq.append(tuple(line.split()[:2]))
            continue
        key, _, value = line.partition(':')
        block = key
        words = [w for w in value.split() if w not in ('!', '""', "''", '|')]
        if key == 'tags':
            tags = words
        elif key == 'skip_tests':
            skip = set(words)
    return tags, skip, prereq


def skipped(case, skip):
    """Whether a test case is skipped in all builds"""
    if case in skip:
        return True
    # the *_full variants of KOKKOS also honor the entry without _full
    if case.endswith('_full') and case[:-5] in skip:
        return True
    # vatom_only_kokkos runs only when the case of the KOKKOS back end of the
    # build (kokkos_serial, kokkos_omp, or kokkos_gpu) runs: that one counts
    return case == 'vatom_only_kokkos'


def has_variants(program, kind, prereq, suffix):
    """Whether all prerequisites that need the suffix have a variant with it"""
    if kind is None:
        return True
    for cat, style in prereq:
        if cat != kind:
            continue
        if kind == 'pair' and style.startswith('hybrid') and suffix != 'kk':
            continue
        if style + '/' + suffix not in styles.get(cat, set()):
            return False
    return True


def wanted(program, kind, skip, prereq):
    """The managed tags a YAML file should have"""
    tags = []
    if 'gpu' in cases[program] and not skipped('gpu', skip) \
       and has_variants(program, kind, prereq, 'gpu'):
        tags.append('gpu')
    kk_cases = [c for c in cases[program] if 'kokkos' in c and not skipped(c, skip)]
    if kk_cases and has_variants(program, kind, prereq, 'kk'):
        tags.append('kokkos')
    return tags


def new_tags(tags, want):
    """The tags line with the managed tags replaced, "generated" stays last"""
    keep = [t for t in tags if t not in MANAGED and t != 'generated']
    return keep + want + (['generated'] if 'generated' in tags else [])


changed = unknown = 0
counts = {t: 0 for t in MANAGED}
for name in sorted(os.listdir(tests_dir)):
    if not name.endswith('.yaml'):
        continue
    group = [p for p in GROUPS if name.startswith(p)]
    if not group:
        unknown += 1
        if args.verbose:
            print(f"{name}: not a file of a force-style test program, not checked")
        continue
    program, kind = GROUPS[group[0]]
    path = os.path.join(tests_dir, name)
    with open(path) as f:
        lines = f.read().splitlines(keepends=True)
    tags, skip, prereq = read_yaml(lines)
    want = wanted(program, kind, skip, prereq)
    for t in want:
        counts[t] += 1
    result = new_tags(tags, want)
    if args.verbose:
        print(f"{name}: {' '.join(want) or '-'}")
    if result == tags:
        continue
    changed += 1
    print(f"{name}: tags '{' '.join(tags)}' -> '{' '.join(result)}'")
    if not args.update:
        continue
    line = ('tags: ' + ' '.join(result)).rstrip() + '\n'
    pos = [i for i, text in enumerate(lines) if text.startswith('tags:')]
    if pos:
        lines[pos[0]] = line
    else:
        # the place where the test programs write it: after lammps_version
        after = [i for i, text in enumerate(lines) if text.startswith('lammps_version:')]
        lines.insert(after[0] + 1 if after else 1, line)
    with open(path, 'w') as f:
        f.write(''.join(lines))

print(f"{counts['gpu']} files with the tag gpu, {counts['kokkos']} with the tag kokkos, "
      f"{unknown} other files not checked")
if changed:
    print(f"{changed} files {'updated' if args.update else 'differ (update them with -u)'}")
if changed and not args.update:
    sys.exit(1)
