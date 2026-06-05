#!/usr/bin/env python3
# Utility for enforcing consistent ID nomenclature in the LAMMPS manual
#
# Standardizes the spelling of identifiers in doc/src/*.rst on the dashed
# form: atom-ID, molecule-ID, group-ID, region-ID, chunk-ID, fix-ID,
# compute-ID, dump-ID.  The command synopsis placeholder is qualified by the
# command name, e.g. "fix fix-ID group-ID nve" instead of "fix ID group-ID nve".
# See https://github.com/lammps/lammps/issues/4598
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
patterns:
    - "*.rst"
"""

# command synopsis lines (only applied inside code-block / parsed-literal blocks):
#   "   fix ID group-ID nve"   -> "   fix fix-ID group-ID nve"
#   "   compute ID group-ID t" -> "   compute compute-ID group-ID t"
#   "   region ID style args"  -> "   region region-ID style args"
SYNOPSIS_RE = re.compile(r'^(\s+)(fix|compute|dump|region|group|molecule) ID\b')

# placeholder definition bullets in the reference command docs:
#   "* ID = user-assigned name for the fix"            -> "* fix-ID = ..."
#   "* ID = user-assigned name for the computation"    -> "* compute-ID = ..."
#   "* ID = user-defined name of the group"            -> "* group-ID = ..."
#   "* ID = user-assigned name for the molecule template" -> "* molecule-ID = ..."
BULLET_DEF_RE = re.compile(
    r'^(\* )ID( = user-(?:assigned|defined) name (?:for|of) the '
    r'(fix|computation|dump|region|group|molecule template)\b)')
BULLET_DEF_MAP = {
    'fix': 'fix',
    'computation': 'compute',
    'dump': 'dump',
    'region': 'region',
    'group': 'group',
    'molecule template': 'molecule',
}

# per-style synopsis explanation bullets:
#   "* ID, group-ID are documented in :doc:`fix <fix>` command"
#       -> "* fix-ID, group-ID are documented in :doc:`fix <fix>` command"
BULLET_DOC_RE = re.compile(
    r'^(\* )ID(, group-ID are documented in (?:the )?:doc:`(fix|compute|dump)\b)')

# inline prose nomenclature, preserving capitalization and singular/plural.
# "ID"/"IDs" are matched case-sensitively (uppercase); the (?!\[) guard keeps
# code references such as "c_ID[I]" untouched (those also lack the leading word).
_WORDS = ['atom', 'molecule', 'group', 'region', 'chunk', 'grid', 'fix', 'compute', 'dump']
_WORD_ALT = '|'.join('[%s%s]%s' % (w[0].upper(), w[0], w[1:]) for w in _WORDS)
PROSE_RE = re.compile(r'\b(%s) (IDs?)\b(?!\[)' % _WORD_ALT)

# code-block / parsed-literal directive
DIRECTIVE_RE = re.compile(r'^(\s*)\.\. +(code-block|parsed-literal)::')


def transform_line(line, in_block):
    # 1. command synopsis lines (inside code blocks only, to avoid hitting
    #    prose list items that happen to start with a command word)
    if in_block:
        line = SYNOPSIS_RE.sub(r'\1\2 \2-ID', line)
    # 2. reference-doc placeholder definition bullets
    line = BULLET_DEF_RE.sub(
        lambda m: m.group(1) + BULLET_DEF_MAP[m.group(3)] + '-ID' + m.group(2), line)
    # 3. per-style "are documented in :doc:" bullets
    line = BULLET_DOC_RE.sub(
        lambda m: m.group(1) + m.group(3) + '-ID' + m.group(2), line)
    # 4. inline prose (entity terms + leftover command references)
    line = PROSE_RE.sub(lambda m: m.group(1) + '-' + m.group(2), line)
    return line


def process_lines(lines):
    """Return (transformed_lines, set_of_changed_linenos)."""
    out = []
    changed = set()
    in_block = False
    block_indent = 0
    lineno = 0
    for line in lines:
        lineno += 1
        stripped = line.rstrip('\n')
        indent = len(stripped) - len(stripped.lstrip(' '))
        blank = (stripped.strip() == '')

        # maintain code-block / parsed-literal state
        if in_block:
            if not blank and indent <= block_indent:
                in_block = False
        if not in_block:
            m = DIRECTIVE_RE.match(stripped)
            if m:
                in_block = True
                block_indent = len(m.group(1))
                out.append(line)
                continue

        newline = transform_line(line, in_block)
        if newline != line:
            changed.add(lineno)
        out.append(newline)
    return out, changed


def check_file(path):
    encoding = 'UTF-8'
    try:
        with open(path, 'r', encoding=encoding) as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        encoding = 'ISO-8859-1'
        try:
            with open(path, 'r', encoding=encoding) as f:
                lines = f.readlines()
        except Exception:
            return {'errors': set(), 'lines': [], 'encoding': 'unknown'}
    newlines, changed = process_lines(lines)
    return {'errors': changed, 'lines': newlines, 'encoding': encoding}


def fix_file(path, result):
    newfile = path + ".modified"
    with open(newfile, 'w', encoding='UTF-8') as out:
        out.writelines(result['lines'])
    shutil.copymode(path, newfile)
    shutil.move(newfile, path)


def check_folder(directory, config, fix=False, verbose=False):
    success = True
    files = []

    for base_path in config['include']:
        for pattern in config['patterns']:
            path = os.path.join(directory, base_path, pattern)
            files += glob.glob(path, recursive=config['recursive'])

    for f in sorted(set(files)):
        path = os.path.normpath(f)

        if verbose:
            print("Checking file:", path)

        result = check_file(path)

        if result['encoding'] == 'unknown':
            print("[Error] Unknown text encoding @ {}".format(path))
            success = False
            continue

        if result['errors']:
            for lineno in sorted(result['errors']):
                print("[Error] Inconsistent ID nomenclature @ {}:{}".format(path, lineno))
            if fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                success = False

    return success


def main():
    parser = argparse.ArgumentParser(
        description='Utility for enforcing consistent ID nomenclature in the LAMMPS manual')
    parser.add_argument('-c', '--config', metavar='CONFIG_FILE', help='location of an optional configuration file')
    parser.add_argument('-f', '--fix', action='store_true', help='automatically apply fixes')
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
        path = os.path.normpath(lammpsdir)
        if args.verbose:
            print("Checking file:", path)
        result = check_file(path)
        if result['errors']:
            for lineno in sorted(result['errors']):
                print("[Error] Inconsistent ID nomenclature @ {}:{}".format(path, lineno))
            if args.fix:
                print("Applying automatic fixes to file:", path)
                fix_file(path, result)
            else:
                sys.exit(1)


if __name__ == "__main__":
    main()
