#!/usr/bin/env python3
# LAMMPS Documentation Utilities
#
# Scan for duplicate and unused anchor labels in documentation files
#
# Copyright (C) 2017 Richard Berger
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import re
import sys
import argparse

# references to labels: :ref:`label`, :ref:`text <label>`, :numref:`...`,
# `text`_, `text <label_>`_, and label_
ref_role_pattern = re.compile(r':(?:ref|numref):`([^`]+)`')
phrase_ref_pattern = re.compile(r'`([^`]+)`_(?!_)')
name_ref_pattern = re.compile(r'(?<![\w`])([A-Za-z0-9][\w.+-]*)_(?![\w`])')
explicit_target_pattern = re.compile(r'<([^<>]+)>\s*$')

# labels that are link targets for URLs from outside the manual and thus
# not reported as unused: errNNNN are the destinations of the URLs printed
# with error messages (see utils::errorurl()), hintNN are the debugging hints
ignore_pattern = re.compile(r'^(err\d{4}|hint\d{2})$')

def normalize(label):
    """labels are case insensitive and whitespace is normalized like in Sphinx"""
    return ' '.join(label.split()).lower()

def find_references(text):
    """return the set of normalized label names referenced in text"""
    refs = set()
    for m in ref_role_pattern.finditer(text):
        target = explicit_target_pattern.search(m.group(1))
        refs.add(normalize(target.group(1) if target else m.group(1)))
    for m in phrase_ref_pattern.finditer(text):
        target = explicit_target_pattern.search(m.group(1))
        if not target:
            refs.add(normalize(m.group(1)))
        elif target.group(1).endswith('_'):
            refs.add(normalize(target.group(1)[:-1]))
    for m in name_ref_pattern.finditer(text):
        refs.add(normalize(m.group(1)))
    return refs

def main():
    parser = argparse.ArgumentParser(description='scan for duplicate and unused anchor labels in documentation files')
    parser.add_argument('files',  metavar='file', nargs='+', help='one or more files to scan')
    parser.add_argument('-r', '--refs', metavar='file', nargs='*', default=[],
                        help='additional files (e.g. C++ or Python sources) that may reference labels')
    parsed_args = parser.parse_args()

    anchor_pattern = re.compile(r'^\.\. _(.*):$')
    anchors = {}
    references = set()

    for filename in parsed_args.files:
        with open(filename, 'rt') as f:
            text = f.read()
        references |= find_references(text)
        for line_number, line in enumerate(text.splitlines()):
            m = anchor_pattern.match(line)
            if m:
                label = m.group(1)
                if label in anchors:
                    anchors[label].append((filename, line_number+1))
                else:
                    anchors[label] = [(filename, line_number+1)]

    for filename in parsed_args.refs:
        with open(filename, 'rt', encoding='utf-8', errors='replace') as f:
            references |= find_references(f.read())

    print("Found %d anchor labels" % len(anchors))

    unused = []
    for label in anchors:
        if ignore_pattern.match(label) or normalize(label.strip('`')) in references:
            continue
        for filename, line_number in anchors[label]:
            unused.append((filename, line_number, label))

    for filename, line_number, label in sorted(unused):
        print("Unused anchor label '%s' @ %s:%d" % (label, filename, line_number))
    if len(unused) > 0:
        print("Found %d unused anchor labels." % len(unused))
    else:
        print("No unused anchor labels.")

    count = 0

    for label in sorted(anchors.keys()):
        if len(anchors[label]) > 1:
            print(label)
            count += 1
            for filename, line_number in anchors[label]:
                print(" - %s:%d" % (filename, line_number))


    if count > 0:
        print("Found %d anchor label errors." % count)
        sys.exit(1)
    else:
        print("No anchor label errors.")

if __name__ == "__main__":
    main()
