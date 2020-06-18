#!/usr/bin/env python3
# Utility for detecting and fixing file permission issues in LAMMPS
#
# Written by Richard Berger (Temple University)
import os
import sys
import glob
import yaml
import argparse
import stat

DEFAULT_CONFIG = """
permission: "rw-r--r--"
recursive: true
include:
    - cmake/**
    - doc/src/**
    - python
    - src/**
    - examples/**
    - tools/coding_standard
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

def check_permission(path, mask):
    st = os.stat(path)
    return bool(stat.S_IMODE(st.st_mode) == mask)

def generate_permission_mask(line):
    assert(len(line) == 9)
    mask = 0

    # USER
    if line[0] == "r":
        mask |= stat.S_IRUSR
    if line[1] == "w":
        mask |= stat.S_IWUSR
    if line[2] == "x":
        mask |= stat.S_IXUSR

    # GROUP
    if line[3] == "r":
        mask |= stat.S_IRGRP
    if line[4] == "w":
        mask |= stat.S_IWGRP
    if line[5] == "x":
        mask |= stat.S_IXGRP

    # OTHER
    if line[6] == "r":
        mask |= stat.S_IROTH
    if line[7] == "w":
        mask |= stat.S_IWOTH
    if line[8] == "x":
        mask |= stat.S_IXOTH

    return mask

def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = []

    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files += glob.glob(path, recursive=config['recursive'])

    mask = generate_permission_mask(config['permission'])

    for f in files:
        path = os.path.normpath(f)

        if verbose:
            print("Checking file:", path)

        ok = check_permission(path, mask)

        if not ok:
            print("[Error] Wrong file permissions @ {}".format(path))

            if fix:
                if os.access(path, os.W_OK):
                    print("Changing permissions of file {} to '{}'".format(path, config['permission']))
                    os.chmod(path, mask)
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

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    if not check_folder(args.DIRECTORY, config, args.fix, args.verbose):
        sys.exit(1)

if __name__ == "__main__":
    main()
