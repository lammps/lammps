#! /usr/bin/env python3
# LAMMPS Documentation Utilities
#
# Scan for duplicate anchor labels in documentation files
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

def main():
    parser = argparse.ArgumentParser(description='scan for duplicate anchor labels in documentation files')
    parser.add_argument('files',  metavar='file', nargs='+', help='one or more files to scan')
    parsed_args = parser.parse_args()

    anchor_pattern = re.compile(r'^\.\. _(.*):$')
    anchors = {}

    for filename in parsed_args.files:
        #print("filename: %s" % filename)
        with open(filename, 'rt') as f:
            for line_number, line in enumerate(f):
                m = anchor_pattern.match(line)
                if m:
                    label = m.group(1)
                    #print("found label: %s" % label)
                    if label in anchors:
                        anchors[label].append((filename, line_number+1))
                    else:
                        anchors[label] = [(filename, line_number+1)]

    print("found %d anchor labels" % len(anchors))

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
