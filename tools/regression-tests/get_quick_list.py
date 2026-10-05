#!/usr/bin/env python3
"""
Find all example input files containing commands changed in this branch versus develop.
Companion script to run_tests.py regression tester.
"""

import glob, os, re, sys, subprocess
from pathlib import Path

if sys.version_info < (3,5):
    raise BaseException("Must use at least Python 3.5")

# infer top level LAMMPS dir from filename
LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

# ----------------------------------------------------------------------

def changed_files_from_git(branch='develop'):
    """
    Return list of changed file from git.

    This function queries git to return the list of changed files on
    the current branch relative to a given branch (default is 'develop').

    param branch: branch to compare with
    type branch: string
    return: path names of files with changes relative to the repository root
    rtype: list of strings
    """

    # get list of changed files relative to the develop branch from git
    output = None
    try:
        output = subprocess.run('git diff --diff-filter=MA --name-status ' + branch,
                                shell=True, capture_output=True)
    except:
        pass

    # collect header files to check for styles
    # - skip files that don't end in '.h' or '.cpp'
    # - skip paths that don't start with 'src/'
    # - replace '.cpp' with '.h' w/o checking it exists
    headers = []
    # output will have a letter 'A' or 'M' for added or modified files followed by pathname
    # append iterms to list and return it
    if output:
        for changed in output.stdout.decode().split():
            if (changed == 'A') or (changed == 'M'): continue
            if not changed.startswith('src/'): continue
            if changed.endswith('.h'): headers.append(changed)
            if changed.endswith('.cpp'): headers.append(changed.replace('.cpp','.h'))
    return headers

# ----------------------------------------------------------------------

def get_command_from_header(headers, topdir="."):
    """
    Loop over list of header files and extract style names, if present.

    LAMMPS commands have macros XxxxStyle() that connects a string with a class.
    We search the header files for those macros, extract the string and append
    it to a list in a dictionary of different types of styles. We skip over known
    suffixes and deprecated commands.

    param headers: header files to check for commands
    type headers:
    return: dictionary with lists of style names
    rtype: dict
    """

    styles = {}
    styles['command'] = []
    styles['atom'] = []
    styles['compute'] = []
    styles['fix'] = []
    styles['pair'] = []
    styles['body'] = []
    styles['bond'] = []
    styles['angle'] = []
    styles['dihedral'] = []
    styles['improper'] = []
    styles['kspace'] = []
    styles['dump'] = []
    styles['region'] = []
    styles['integrate'] = []
    styles['minimize'] = []

    # some regex
    style_pattern = re.compile(r"(.+)Style\((.+),(.+)\)")
    upper = re.compile("[A-Z]+")
    gpu = re.compile("(.+)/gpu$")
    intel = re.compile("(.+)/intel$")
    kokkos = re.compile("(.+)/kk$")
    kokkos_skip = re.compile("(.+)/kk/(host|device)$")
    omp = re.compile("(.+)/omp$")
    opt = re.compile("(.+)/opt$")
    removed = re.compile("(.*)Deprecated$")

    for file in headers:
        # don't fail if file is not present
        try:
            with open(os.path.join(topdir,file)) as f:
                for line in f:
                    matches = style_pattern.findall(line)
                    for m in matches:
                        # skip over internal styles w/o explicit documentation
                        style = m[1]
                        if upper.match(style):
                            continue

                        # skip over suffix styles:
                        suffix = kokkos_skip.match(style)
                        if suffix:
                            continue
                        suffix = gpu.match(style)
                        if suffix:
                            continue
                        suffix = intel.match(style)
                        if suffix:
                            continue
                        suffix = kokkos.match(style)
                        if suffix:
                            continue
                        suffix = omp.match(style)
                        if suffix:
                            continue
                        suffix = opt.match(style)
                        if suffix:
                            continue
                        deprecated = removed.match(m[2])
                        if deprecated:
                            continue

                        # register style and suffix flags
                        if m[0] == 'Angle':
                            styles['angle'].append(style)
                        elif m[0] == 'Atom':
                            styles['atom'].append(style)
                        elif m[0] == 'Body':
                            styles['body'].append(style)
                        elif m[0] == 'Bond':
                            styles['bond'].append(style)
                        elif m[0] == 'Command':
                            styles['command'].append(style)
                        elif m[0] == 'Compute':
                            styles['compute'].append(style)
                        elif m[0] == 'Dihedral':
                            styles['dihedral'].append(style)
                        elif m[0] == 'Dump':
                            styles['dump'].append(style)
                        elif m[0] == 'Fix':
                            styles['fix'].append(style)
                        elif m[0] == 'Improper':
                            styles['improper'].append(style)
                        elif m[0] == 'Integrate':
                            styles['integrate'].append(style)
                        elif m[0] == 'KSpace':
                            styles['kspace'].append(style)
                        elif m[0] == 'Minimize':
                            styles['minimize'].append(style)
                        elif m[0] == 'Pair':
                            styles['pair'].append(style)
                        elif m[0] == 'Region':
                            styles['region'].append(style)
                        else:
                            pass
        # header file not found or not readable
        except:
            pass
    return styles

# ----------------------------------------------------------------------

# style macro prefixes and the corresponding style category, i.e. the kind of
# XxxStyle(name,Class) macros in the headers that register input script styles
STYLE_MACRO_CATEGORY = {
    'Angle': 'angle', 'Atom': 'atom', 'Body': 'body', 'Bond': 'bond',
    'Command': 'command', 'Compute': 'compute', 'Dihedral': 'dihedral',
    'Dump': 'dump', 'Fix': 'fix', 'Improper': 'improper', 'Integrate': 'integrate',
    'KSpace': 'kspace', 'Minimize': 'minimize', 'Pair': 'pair', 'Region': 'region',
}

def get_style_universe(topdir):
    """
    Collect all style names that exist anywhere in the LAMMPS distribution.

    Scans the headers in src/ and all package folders under src/ for the
    XxxStyle(name,Class) macros, including the packages that are not installed
    in any given build.  The same styles are hidden as in the -help output of
    the LAMMPS binary: internal styles starting with an upper-case letter and
    the KOKKOS host/device variants.  Comparing this universe with the styles
    reported by "lmp -help" tells which known styles a binary does not include.

    param topdir: top-level folder of the LAMMPS source tree
    type topdir: string
    return: nested dictionary: category -> style name -> package folder name
            (an empty string marks a style from the core in src/)
    rtype: dict
    """

    style_macro = re.compile(r'^\s*([A-Za-z]+)Style\(([^,\s]+)\s*,\s*[^)\s]+\)')
    kokkos_internal = re.compile(r'.*/kk/(host|device)$')
    package_dir = re.compile(r'^[A-Z][A-Z0-9-]*$')

    universe = {category: {} for category in STYLE_MACRO_CATEGORY.values()}
    src = os.path.join(topdir, 'src')
    folders = [('', src)]
    try:
        entries = sorted(os.listdir(src))
    except OSError:
        entries = []
    for entry in entries:
        if package_dir.match(entry) and entry not in ('STUBS', 'MAKE', 'DEPEND'):
            folders.append((entry, os.path.join(src, entry)))

    for package, folder in folders:
        for header in sorted(glob.glob(os.path.join(folder, '*.h'))):
            try:
                with open(header) as f:
                    for line in f:
                        m = style_macro.match(line)
                        if not m:
                            continue
                        category = STYLE_MACRO_CATEGORY.get(m.group(1))
                        style = m.group(2)
                        if not category:
                            continue
                        # skip the same styles that the -help output hides
                        if style[0].isupper() or kokkos_internal.match(style):
                            continue
                        universe[category].setdefault(style, package)
            except OSError:
                pass
    return universe

# ----------------------------------------------------------------------

def make_regex(styles):
    """Convert dictionary with styles into a regular expression to scan input files with

    This will construct a regular expression matching LAMMPS commands. Ignores continuation

    param styles: dictionary with style names
    type styles: dict
    return: combined regular expression string
    rtype: string
    """

    restring = "^\\s*("
    if len(styles['command']):
        restring += '(' + '|'.join(styles['command']) + ')|'
    if len(styles['atom']):
        restring += '(atom_style\\s+(' + '|'.join(styles['atom']) + '))|'
    if len(styles['compute']):
        restring += '(compute\\s+\\S+\\s+\\S+\\s+(' + '|'.join(styles['compute']) + '))|'
    if len(styles['fix']):
        restring += '(fix\\s+\\S+\\s+\\S+\\s+(' + '|'.join(styles['fix']) + '))|'
    if len(styles['pair']):
        restring += '(pair_style\\s+(' + '|'.join(styles['pair']) + '))|'
    if len(styles['body']):
        restring += '(atom_style\\s+body\\s+(' + '|'.join(styles['body']) + '))|'
    if len(styles['bond']):
        restring += '(bond_style\\s+(' + '|'.join(styles['bond']) + '))|'
    if len(styles['angle']):
        restring += '(angle_style\\s+(' + '|'.join(styles['angle']) + '))|'
    if len(styles['dihedral']):
        restring += '(dihedral_style\\s+(' + '|'.join(styles['dihedral']) + '))|'
    if len(styles['improper']):
        restring += '(improper_style\\s+(' + '|'.join(styles['improper']) + '))|'
    if len(styles['kspace']):
        restring += '(kspace_style\\s+(' + '|'.join(styles['kspace']) + '))|'
    if len(styles['dump']):
        restring += '(dump\\s+\\S+\\s+\\S+\\s+(' + '|'.join(styles['dump']) + '))|'
    if len(styles['region']):
        restring += '(region\\s+(' + '|'.join(styles['region']) + '))|'
    if len(styles['integrate']):
        restring += '(run_style\\s+(' + '|'.join(styles['integrate']) + '))|'
    if len(styles['minimize']):
        restring += '(min_style\\s+(' + '|'.join(styles['minimize']) + '))|'

    # replace last (pipe) character with closing parenthesis
    length = len(restring)
    restring = restring[:length-1] + ')'
    # return combined regex string
    if length > 5:
        return restring
    else:
        return None

# ----------------------------------------------------------------------

def get_examples_using_styles(regex, examples='examples'):
    """
    Loop through LAMMPS examples tree and find all files staring with 'in.'
    that have at least one line matching the regex.

    param regex: string pattern matching LAMMPS commands
    type regex: compiled regex
    param example: path where to start looking for examples recursively
    type example: string
    return: list of matching example inputs
    rtype: list of strings
    """

    commands = re.compile(regex)
    inputs = []
    for filename in Path(examples).rglob('in.*'):
        # leave out editor backup copies (e.g. in.melt~) of input scripts
        if str(filename).endswith('~'):
            continue
        with open(filename) as f:
            for line in f:
                if commands.match(line):
                    inputs.append(str(filename))
                    break
    return inputs

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__":

    headers = changed_files_from_git('origin/develop')
    styles = get_command_from_header(headers, LAMMPS_DIR)
    regex = make_regex(styles)
    if regex:
        inputs = get_examples_using_styles(regex, os.path.join(LAMMPS_DIR,'examples'))
    else:
        inputs = []
    print("Found changes to the following styles:")
    print("Commands: ", styles['command'])
    print("Atom styles: ", styles['atom'])
    print("Compute styles: ", styles['compute'])
    print("Fix styles: ", styles['fix'])
    print("Pair styles: ", styles['pair'])
    print("Body styles: ", styles['body'])
    print("Bond styles: ", styles['bond'])
    print("Angle styles: ", styles['angle'])
    print("Dihedral styles: ", styles['dihedral'])
    print("Improper styles: ", styles['improper'])
    print("Kspace styles: ", styles['kspace'])
    print("Dump styles: ", styles['dump'])
    print("Region styles: ", styles['region'])
    print("Integrate styles: ", styles['integrate'])
    print("Minimize styles: ", styles['minimize'])

    print("Example input files affected: ", len(inputs))
    print("inputs: ", inputs.sort())
