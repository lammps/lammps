#!/usr/bin/env python3
# Utility for detecting and setting `.. versionadded::` and `.. versionchanged::``
#
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
    - doc/src/**
    - python
    - src/**
exclude:
    - src/Make.sh
patterns:
    - "*.c"
    - "*.cpp"
    - "*.h"
    - "*.md"
    - "*.py"
    - "*.rst"
    - "*.f90"
"""

def check_pending_tag(f):
    pattern = re.compile(r'^ *\.\. +(version(changed|added)|deprecated):: +TBD')
    last_line = "\n"
    lineno = 1
    errors = set()

    for line in f:
        if pattern.match(line):
            errors.add(lineno)
        last_line = line
        lineno += 1

    return errors, last_line

def check_file(path):
    encoding = 'UTF-8'
    pending_tags = set()
    try:
        with open(path, 'r') as f:
            pending_tags, last_line = check_pending_tag(f)
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                pending_tags, last_line = check_pending_tag(f)
        except Exception:
            encoding = 'unknown'

    return {
        'pending_tags': pending_tags,
        'encoding': encoding
    }

def check_folder(directory, config, verbose=False):
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

        for lineno in result['pending_tags']:
            print("[Error] Pending version tag @ {}:{}".format(path, lineno))
            success = False

        if result['encoding'] == 'unknown':
            print("[Error] Unknown text encoding @ {}".format(path))
            success = False
        elif result['encoding'] == 'ISO-8859-1':
            print("[Error] Found ISO-8859-1 encoding instead of UTF-8 @ {}".format(path))

    return success

def main():
    parser = argparse.ArgumentParser(description='Utility for detecting pending version markers in LAMMPS')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('DIRECTORY', help='directory that should be checked')
    args = parser.parse_args()
    lammpsdir = os.path.abspath(os.path.expanduser(args.DIRECTORY))

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    if not check_folder(lammpsdir, config, args.verbose):
        sys.exit(1)

if __name__ == "__main__":
    main()
