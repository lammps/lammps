#!/usr/bin/env python3
# Utility for detecting fmtlib related issues
#
# Currently it checks for the following issues
# use of fmt::print() instead of utils::print()
#
# Written by Axel Kohlmeyer (Temple University)
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
    - src/**
exclude:
    - "src/fmt/"
    - "src/fmtlib"
    - "src/utils"
patterns:
    - "*.h"
    - "*.cpp"
"""

def check_fmtprint(f):
    pattern = re.compile(r'[ \t\n\r]*fmt::print\(')
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        lineno += 1

    return errors

def check_file(path):
    if path.find('fmtlib.py') >= 0: return { 'fmtlib_errors' : '' }
    encoding = 'UTF-8'
    fmtprint_errors = set()
    try:
        with open(path, 'r') as f:
            fmtprint_errors = check_fmtprint(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                fmtprint_errors = check_fmtprint(f)
        except Exception:
            encoding = 'unknown'

    return {
        'fmtprint_errors': fmtprint_errors,
        'encoding': encoding
    }

def fix_file(path, check_result):
    if path.find('fmtlib.py') >= 0: return
    newfile = path + ".modified"
    pattern = re.compile(r'fmt::print\(', re.DOTALL)
    with open(newfile, 'w', encoding='UTF-8') as out:
        with open(path, 'r', encoding=check_result['encoding']) as src:
            filetxt = re.sub(pattern,'utils::print(', src.read());
            print(filetxt, end='', file=out)
    shutil.copymode(path, newfile)
    shutil.move(newfile, path)

def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = []

    # compile list of files to check
    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files += glob.glob(path, recursive=config['recursive'])

    # prune list of files to skip from list
    for pattern in config['exclude']:
        path = os.path.join(directory, pattern)
        remove = []
        for file in files:
            if path not in file: continue
            remove += [file]
        for rm in remove:
            files.remove(rm)

    for f in files:
        path = os.path.normpath(f)

        if verbose:
            print("Checking file:", path)

        result = check_file(path)

        has_resolvable_errors = False

        for lineno in result['fmtprint_errors']:
            print("[Error] Found LAMMPS fmt::print @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        if has_resolvable_errors:
            if fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                success = False

    return success

def main():
    parser = argparse.ArgumentParser(description='Utility for detecting and fixing fmtlib issues in LAMMPS')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='automatically fix URLs')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('DIRECTORY', help='directory (or file) that should be checked')
    args = parser.parse_args()
    lammpsdir = os.path.abspath(os.path.expanduser(args.DIRECTORY))

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    if os.path.isdir(lammpsdir):
        if not check_folder(lammpsdir, config, args.fix, args.verbose):
           sys.exit(1)
    else:
        success = True
        path = os.path.normpath(lammpsdir)

        if args.verbose:
            print("Checking file:", path)

        result = check_file(path)

        has_resolvable_errors = False

        for lineno in result['fmtprint_errors']:
            print("[Error] Found LAMMPS fmt::print @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        if has_resolvable_errors:
            if args.fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                success = False

        if not success:
            sys.exit(1)

if __name__ == "__main__":
    main()
