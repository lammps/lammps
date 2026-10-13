#!/usr/bin/env python3
"""Report which LAMMPS styles are covered by unit tests.

A style is counted in one of three groups:
  - it has a YAML reference test in the force-style test folder,
  - it is only used by some other unit test (other YAML harnesses like
    bpm, granular, or graphics, or the commands of C++, Python, and
    Fortran unit test sources and the LAMMPS inputs they read),
  - it is not used by any unit test.
Aliases (several style names for the same class) share their status;
deprecated and removed styles are not counted.
"""

import os, re, sys
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
parser.add_argument("-u", "--unittest",
                    help="Path to the LAMMPS unit test tree")

args = parser.parse_args()
verbose = args.verbose
src_dir = args.src
tests_dir = args.tests
unittest_dir = args.unittest

LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

if not src_dir:
    src_dir = os.path.join(LAMMPS_DIR , 'src')

if not tests_dir:
    tests_dir = os.path.join(LAMMPS_DIR, 'unittest', 'force-styles', 'tests')

if not unittest_dir:
    unittest_dir = os.path.join(LAMMPS_DIR, 'unittest')

try:
    src_dir = os.path.abspath(os.path.expanduser(src_dir))
    tests_dir = os.path.abspath(os.path.expanduser(tests_dir))
    unittest_dir = os.path.abspath(os.path.expanduser(unittest_dir))
except:                                 # lgtm [py/catch-base-exception]
    parser.print_help()
    sys.exit(1)

if not os.path.isdir(src_dir):
    sys.exit(f"LAMMPS source path {src_dir} does not exist")

if not os.path.isdir(tests_dir):
    sys.exit(f"LAMMPS test inputs path {tests_dir} does not exist")

if not os.path.isdir(unittest_dir):
    sys.exit(f"LAMMPS unit test path {unittest_dir} does not exist")

# style categories: name of the style macro prefix -> category
CATEGORIES = {'Pair': 'pair', 'Bond': 'bond', 'Angle': 'angle', 'Dihedral': 'dihedral',
              'Improper': 'improper', 'KSpace': 'kspace', 'Fix': 'fix',
              'Compute': 'compute', 'Minimize': 'min'}
# categories that have a hybrid style with sub-styles as arguments
HYBRID = ('pair', 'bond', 'angle', 'dihedral', 'improper')

style_pattern = re.compile(r'^\s*(\w+)Style\(\s*([\w/]+)\s*,\s*([^)]+?)\s*\)', re.M)
suffix_pattern = re.compile(r'(.+)/(gpu|intel|kk|omp|opt)(/host|/device)?$')
deprecated_pattern = re.compile(r'.*Deprecated$')

# per category: style name -> dict with class names, suffixes, deprecated flag
styles = {cat: {} for cat in CATEGORIES.values()}

def split_suffix(name):
    """Return base style name, accelerator suffix, and kokkos host/device flag"""
    m = suffix_pattern.match(name)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return name, None, None

print("Parsing style names from C++ tree in:    ", src_dir)

headers = []
for path, dirs, files in os.walk(src_dir):
    headers += [os.path.join(path, f) for f in files if f.endswith('.h')]

for header in sorted(headers):
    with open(header, errors='replace') as f:
        text = f.read()
    for kind, name, cls in style_pattern.findall(text):
        if kind not in CATEGORIES:
            continue
        # skip over internal styles w/o explicit documentation
        if name[0].isupper():
            continue
        base, suffix, hostdev = split_suffix(name)
        # /kk/host and /kk/device are the same class as /kk
        if hostdev:
            continue
        info = styles[CATEGORIES[kind]].setdefault(base, {'classes': set(), 'suffixes': set(),
                                                          'deprecated': False})
        if suffix:
            info['suffixes'].add(suffix)
        else:
            info['classes'].add(cls)
            if deprecated_pattern.match(cls):
                info['deprecated'] = True

# style names for the same class are aliases of each other
aliases = {}
for cat, registry in styles.items():
    by_class = {}
    for name, info in registry.items():
        if info['deprecated']:
            continue
        for cls in info['classes']:
            by_class.setdefault(cls, set()).add(name)
    aliases[cat] = {name: set() for name in registry}
    for names in by_class.values():
        for name in names:
            aliases[cat][name] |= names - {name}

# where a style was found: category -> style name -> set of files
reference = {cat: {} for cat in styles}
other = {cat: {} for cat in styles}

def record(found, cat, text, source):
    """Record the style (and hybrid sub-styles) at the start of text"""
    words = text.split()
    if not words:
        return
    if (cat in HYBRID) and words[0].startswith('hybrid'):
        candidates = words
    else:
        candidates = words[:1]
    for word in candidates:
        base, suffix, hostdev = split_suffix(word.strip('\'"'))
        if base in styles[cat]:
            found[cat].setdefault(base, set()).add(source)

# YAML reference tests: file pattern, category, and the line naming the style under test
reference_tests = [
    (re.compile(r'.+-pair-.+\.yaml$'), 'pair', re.compile(r'^pair_style:\s*(.*)$', re.M)),
    (re.compile(r'bond-.+\.yaml$'), 'bond', re.compile(r'^bond_style:\s*(.*)$', re.M)),
    (re.compile(r'angle-.+\.yaml$'), 'angle', re.compile(r'^angle_style:\s*(.*)$', re.M)),
    (re.compile(r'dihedral-.+\.yaml$'), 'dihedral',
     re.compile(r'^dihedral_style:\s*(.*)$', re.M)),
    (re.compile(r'improper-.+\.yaml$'), 'improper',
     re.compile(r'^improper_style:\s*(.*)$', re.M)),
    (re.compile(r'kspace-.+\.yaml$'), 'kspace', re.compile(r'^\s*kspace_style\s+(.*)$', re.M)),
    # the fix or compute style under test is the one created with the ID "test"; the
    # command may also be quoted, e.g. inside an "if" command
    (re.compile(r'fix-.+\.yaml$'), 'fix',
     re.compile(r'(?:^|["\'])\s*fix[ \t]+test[ \t]+\S+[ \t]+(.*)$', re.M)),
    (re.compile(r'compute-.+\.yaml$'), 'compute',
     re.compile(r'(?:^|["\'])\s*compute[ \t]+test[ \t]+\S+[ \t]+(.*)$', re.M)),
    (re.compile(r'min-.+\.yaml$'), 'min', re.compile(r'^\s*min_style\s+(.*)$', re.M)),
]

print("Parsing force-style YAML tests in:       ", tests_dir)

for yaml in sorted(os.listdir(tests_dir)):
    for pattern, cat, search in reference_tests:
        if pattern.match(yaml):
            with open(os.path.join(tests_dir, yaml), errors='replace') as f:
                text = f.read()
            for m in search.findall(text):
                record(reference, cat, m, yaml)

# commands in YAML files, LAMMPS inputs, and string constants in test sources.
# a command must be at the start of a line, of a string, or after an escaped newline.
cmd_start = r'(?:^|["\']|\\n)[ \t]*'
cmd_args = r'([^"\'\\\n]*)'
style_cmd = re.compile(cmd_start + r'(pair|bond|angle|dihedral|improper|kspace|min)_style'
                       + r'(?:[ \t]+|:[ \t]*)' + cmd_args, re.M)
fixcompute_cmd = re.compile(cmd_start + r'(fix|compute)[ \t]+[^\s"\']+[ \t]+[^\s"\']+[ \t]+'
                            + cmd_args, re.M)
# the Python module "cmd" wrapper: lmp.cmd.pair_style("lj/cut", 2.5)
# or lmp.cmd.fix("1", "all", "nve") and lmp.cmd.fix("1 all nve")
wrapper_style = re.compile(r'\.cmd\.(pair|bond|angle|dihedral|improper|kspace|min)_style'
                           + r'\(\s*["\']([^"\']*)')
wrapper_fixcompute = re.compile(r'\.cmd\.(fix|compute)\(\s*["\']([^"\']*)["\']'
                                + r'(?:\s*,\s*["\']([^"\']*)["\']\s*,\s*["\']([^"\']*))?')
# prerequisites block in YAML test files
prereq_block = re.compile(r'^prerequisites:[^\n]*\n((?:[ \t]+\S[^\n]*\n)*)', re.M)
prereq_line = re.compile(r'^[ \t]+(\w+)[ \t]+(\S+)', re.M)
prereq_category = {'pair': 'pair', 'bond': 'bond', 'angle': 'angle', 'dihedral': 'dihedral',
                   'improper': 'improper', 'kspace': 'kspace', 'fix': 'fix',
                   'compute': 'compute', 'minimize': 'min'}

source_ext = ('.yaml', '.cpp', '.h', '.c', '.py', '.f90')
skip_dirs = ('Testing', '__pycache__', '.git')

print("Parsing commands in unit tests in:       ", unittest_dir)

for path, dirs, files in os.walk(unittest_dir):
    dirs[:] = [d for d in dirs if d not in skip_dirs]
    for name in sorted(files):
        if not (name.endswith(source_ext) or name.startswith('in.')):
            continue
        fullname = os.path.join(path, name)
        if os.path.samefile(fullname, __file__):
            continue
        source = os.path.relpath(fullname, unittest_dir)
        with open(fullname, errors='replace') as f:
            text = f.read()
        for cat, cmd in style_cmd.findall(text):
            record(other, cat, cmd, source)
        for kind, cmd in fixcompute_cmd.findall(text):
            record(other, kind, cmd, source)
        for cat, style in wrapper_style.findall(text):
            record(other, cat, style, source)
        for kind, first, group, style in wrapper_fixcompute.findall(text):
            if style:
                record(other, kind, style, source)
            else:
                words = first.split()
                if len(words) > 2:
                    record(other, kind, ' '.join(words[2:]), source)
        if name.endswith('.yaml'):
            for block in prereq_block.findall(text):
                for kind, style in prereq_line.findall(block):
                    if kind in prereq_category:
                        record(other, prereq_category[kind], style, source)

def covered(cat, name, found):
    """Return the style name or an alias of it that was found, or None"""
    if name in found[cat]:
        return name
    for alias in sorted(aliases[cat][name]):
        if alias in found[cat]:
            return alias
    return None

display = {cat: kind for kind, cat in CATEGORIES.items()}
summary = []
for cat, registry in styles.items():
    num_deprecated = 0
    with_reference = []
    only_other = []
    missing = []
    for name in sorted(registry):
        if registry[name]['deprecated']:
            num_deprecated += 1
            continue
        match = covered(cat, name, reference)
        if match:
            with_reference.append(name)
            if verbose and (match != name):
                print(f"{cat} style {name} is tested as alias {match}")
            continue
        match = covered(cat, name, other)
        if match:
            only_other.append(name)
            if verbose:
                alias = f" as alias {match}" if (match != name) else ""
                print(f"{cat} style {name} is used{alias} in: "
                      + ", ".join(sorted(other[cat][match])))
            continue
        missing.append(name)

    total = len(registry) - num_deprecated
    summary.append((cat, total, len(with_reference), len(only_other), len(missing)))
    print(f"\n{display[cat]} styles: {total} total, {len(with_reference)} with YAML "
          f"reference test, {len(only_other)} used only in other unit tests, "
          f"{len(missing)} without any test")
    if num_deprecated:
        print(f"Not counted: {num_deprecated} deprecated or removed {cat} styles")
    print("Used only in other unit tests: ", only_other)
    print("No tests for: ", missing)

print(f"\n{'Summary':<10} {'total':>7} {'YAML':>7} {'other':>7} {'none':>7}")
for row in summary:
    print(f"{row[0]:<10} {row[1]:>7} {row[2]:>7} {row[3]:>7} {row[4]:>7}")
totals = [sum(row[i] for row in summary) for i in range(1, 5)]
print(f"{'all':<10} {totals[0]:>7} {totals[1]:>7} {totals[2]:>7} {totals[3]:>7}")
print(f"\nTotal styles without any test: {totals[3]} of {totals[0]}")
