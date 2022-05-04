#!/usr/bin/env python3
# Utility for detecting leftover error/warning docs in header files
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
patterns:
    - "*.h"
"""

def check_errordocs(f):
    pattern = re.compile(r'^ */\* *ERROR/WARNING messag')
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        lineno += 1

    return errors

def check_file(path):
    if path.find('errordocs.py') >= 0: return { 'errordocs_errors' : '' }
    encoding = 'UTF-8'
    errordocs_errors = set()
    try:
        with open(path, 'r') as f:
            errordocs_errors = check_errordocs(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                errordocs_errors = check_errordocs(f)
        except Exception:
            encoding = 'unknown'

    return {
        'errordocs_errors': errordocs_errors,
        'encoding': encoding
    }

def fix_file(path, check_result):
    if path.find('errordocs.py') >= 0: return
    newfile = path + ".modified"
    pattern = re.compile(r'[\n\r ]* */\* *ERROR/WARNING *messag.*\*/', re.DOTALL)
    with open(newfile, 'w', encoding='UTF-8') as out:
        with open(path, 'r', encoding=check_result['encoding']) as src:
            filetxt = re.sub(pattern,'', src.read());
            print(filetxt, end='', file=out)
    shutil.copymode(path, newfile)
    shutil.move(newfile, path)

def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = []

    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files += glob.glob(path, recursive=config['recursive'])

    for f in files:
        path = os.path.normpath(f)

        if verbose:
            print("Checking file:", path)

        result = check_file(path)

        has_resolvable_errors = False

        for lineno in result['errordocs_errors']:
            print("[Error] Found LAMMPS errordocs @ {}:{}".format(path, lineno))
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

        for lineno in result['errordocs_errors']:
            print("[Error] Found LAMMPS errordocs @ {}:{}".format(path, lineno))
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
