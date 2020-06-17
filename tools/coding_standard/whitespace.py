#!/usr/bin/env python3
# Utility for detecting and fixing whitespace issues in LAMMPS
#
# Written by Richard Berger (Temple University)
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
    - python
    - src/**
    - tools/coding_standard
patterns:
    - "*.c"
    - "*.cmake"
    - "*.cpp"
    - "*.h"
    - "*.md"
    - "*.py"
    - "*.rst"
    - "*.sh"
    - ".gitignore"
    - "README"
    - "requirements.txt"
"""

def check_trailing_whitespace(f):
    pattern = re.compile(r'\s+\n$')
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
    last_line = "\n"
    whitespace_errors = set()
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

    return {
        'whitespace_errors': whitespace_errors,
        'encoding': encoding,
        'eof_error': not last_line.endswith('\n')
    }

def fix_file(path, check_result):
    newfile = path + ".modified"
    with open(newfile, 'w', encoding='UTF-8') as out:
        with open(path, 'r', encoding=check_result['encoding']) as src:
            for line in src:
                print(line.rstrip(), file=out)
    shutil.move(newfile, path)

def check_folder(directory, config, fix=False, verbose=False):
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

        for lineno in result['whitespace_errors']:
            print("[Error] Trailing whitespace @ {}:{}".format(path, lineno))
            has_resolvable_errors = True

        if result['eof_error']:
            print("[Error] Missing newline at end of file @ {}".format(path))
            has_resolvable_errors = True

        if result['encoding'] == 'unknown':
            print("[Error] Unknown text encoding @ {}".format(path))
            has_resolvable_errors = False
        elif result['encoding'] == 'ISO-8859-1':
            print("[Error] Found ISO-8859-1 encoding instead of UTF-8 @ {}".format(path))
            has_resolvable_errors = True

        if has_resolvable_errors and fix:
            print("Applying automatic fixes to file:", path)
            fix_file(path, result)


def main():
    parser = argparse.ArgumentParser(description='Utility for detecting and fixing whitespace issues in LAMMPS')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of a optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='automatically fix common issues')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('DIRECTORY', help='directory that should be checked')
    args = parser.parse_args()

    if args.config:
        with open(args.config, 'r') as cfile:
            config = yaml.load(cfile, Loader=yaml.FullLoader)
    else:
        config = yaml.load(DEFAULT_CONFIG, Loader=yaml.FullLoader)

    check_folder(args.DIRECTORY, config, args.fix, args.verbose)

if __name__ == "__main__":
    main()
