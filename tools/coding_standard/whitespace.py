#!/usr/bin/env python3
# Utility for detecting and fixing whitespace issues in LAMMPS
#
# Written by Richard Berger (Temple University)
from __future__ import print_function
import sys

if sys.version_info.major < 3:
    sys.exit('This script must be run with Python 3.5 or later')

if sys.version_info.minor < 5:
    sys.exit('This script must be run with Python 3.5 or later')

import os
import glob
import re
import yaml
import argparse
import shutil

DEFAULT_CONFIG = """
recursive: true
include:
    - cmake/**
    - doc
    - doc/src/**
    - fortran/**
    - python/**
    - src/**
    - lib/**
    - tools/coding_standard
    - tools/python
    - tools/lammps-gui
    - unittest/**
exclude:
    - lib/colvars/Install.py
    - lib/gpu/geryon/file_to_cstr.sh
    - lib/hdnnp
    - lib/kim
    - lib/kokkos
    - lib/latte
    - lib/machdyn
    - lib/mdi
    - lib/mscg
    - lib/pace
    - lib/plumed
    - lib/quip
    - lib/scafacos
    - lib/voronoi
    - src/Make.sh
patterns:
    - "*.c"
    - "*.cmake"
    - "*.cpp"
    - "*.h"
    - "*.md"
    - "*.py"
    - "*.rst"
    - "*.sh"
    - "*.f90"
    - ".gitignore"
    - "README"
    - "requirements.txt"
"""

def check_trailing_whitespace(f):
    pattern = re.compile(r'[^\n]*\s+\n$')
    last_line = "\n"
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        last_line = line
        lineno += 1

    return errors, last_line

def check_tabs(f):
    pattern = re.compile(r'[^\n]*\t+[^n]*\n$')
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        lineno += 1
    return errors

def check_file(path):
    encoding = 'UTF-8'
    last_line = "\n"
    whitespace_errors = set()
    tab_errors = set()
    try:
        with open(path, 'r') as f:
            whitespace_errors, last_line = check_trailing_whitespace(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                whitespace_errors, last_line = check_trailing_whitespace(f)
        except Exception:
            encoding = 'unknown'

    try:
        with open(path, 'r') as f:
            tab_errors = check_tabs(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                tab_errors = check_tabs(f)
        except Exception:
            encoding = 'unknown'

    return {
        'tab_errors': tab_errors,
        'whitespace_errors': whitespace_errors,
        'encoding': encoding,
        'eof_error': not last_line.endswith('\n')
    }

def fix_file(path, check_result):
    newfile = path + ".modified"
    tab_pat = re.compile(r'^([^\t]*)(\t+)(.*)$')
    with open(newfile, 'w', encoding='UTF-8') as out:
        with open(path, 'r', encoding=check_result['encoding']) as src:
            for line in src:
                match = tab_pat.match(line)
                if match:
                    # compute number of blanks assuming 8 character tab setting
                    num = 8*len(match.group(2))-len(match.group(1))%8
                    line = match.group(1) + " "*num + match.group(3)
                print(line.rstrip(), file=out)
    shutil.copymode(path, newfile)
    shutil.move(newfile, path)

def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = []

    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files += glob.glob(path, recursive=config['recursive'])
    for exclude in config['exclude']:
        files = [f for f in files if not f.startswith(os.path.join(directory,exclude))]

    for f in files:
        path = os.path.normpath(f)

        if verbose:
            print("Checking file:", path)

        result = check_file(path)

        has_resolvable_errors = False

        for lineno in result['whitespace_errors']:
            print("[Error] Trailing whitespace @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        for lineno in result['tab_errors']:
            print("[Error] Tab @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        if result['eof_error']:
            print("[Error] Missing newline at end of file @ {}".format(path))
            has_resolvable_errors = True

        if result['encoding'] == 'unknown':
            print("[Error] Unknown text encoding @ {}".format(path))
            has_resolvable_errors = False
            success = False
        elif result['encoding'] == 'ISO-8859-1':
            print("[Error] Found ISO-8859-1 encoding instead of UTF-8 @ {}".format(path))
            has_resolvable_errors = True

        if has_resolvable_errors:
            if fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                success = False

    return success

def main():
    parser = argparse.ArgumentParser(description='Utility for detecting and fixing whitespace issues in LAMMPS')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='automatically fix common issues')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('DIRECTORY', help='directory that should be checked')
    args = parser.parse_args()
    lammpsdir = os.path.abspath(os.path.expanduser(args.DIRECTORY))

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    if not check_folder(lammpsdir, config, args.fix, args.verbose):
        sys.exit(1)

if __name__ == "__main__":
    main()
