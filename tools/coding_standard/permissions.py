#!/usr/bin/env python3
# Utility for detecting and fixing file permission issues in LAMMPS
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
import yaml
import argparse
import stat

DEFAULT_CONFIG = """
recursive: true
include:
    - cmake/**
    - doc/src/**
    - python
    - src/**
    - examples/**
    - tools/coding_standard
    - tools/lammps-gui
patterns:
    - "*.c"
    - "*.cmake"
    - "*.cpp"
    - "*.h"
    - "*.jpg"
    - "*.md"
    - "*.pdf"
    - "*.png"
    - "*.rst"
    - "*.tex"
    - ".gitignore"
    - "README"
    - "in.*"
    - "requirements.txt"
"""

def has_executable_bit(path):
    st = os.stat(path)
    return bool(st.st_mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH))


def get_non_exec_mask(path):
    st = os.stat(path)
    return st.st_mode & ~(stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

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

        if has_executable_bit(path):
            print("[Error] Wrong file permissions @ {}".format(path))

            if fix:
                if os.access(path, os.W_OK):
                    print("Removing executable bit of file {}".format(path))
                    new_mask = get_non_exec_mask(path)
                    os.chmod(path, new_mask)
                else:
                    print("[Error] Can not write permissions of file {}".format(path))
                    success = False
            else:
                success = False
    return success


def main():
    parser = argparse.ArgumentParser(description='Utility for detecting and fixing file permission issues in LAMMPS')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='automatically fix permissions')
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
