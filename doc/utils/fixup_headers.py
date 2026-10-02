#!/usr/bin/env python3
# Utility script to fix headers in doc pages and generate "Accelerator Styles" portion
import os
import shutil
import re
import argparse

index_pattern = re.compile(r"^.. index:: (compute|fix|pair_style|angle_style|bond_style|dihedral_style|improper_style|kspace_style)\s+([a-zA-Z0-9/_]+)$")
pattern = re.compile(r"^(compute|fix|pair_style|angle_style|bond_style|dihedral_style|improper_style|kspace_style)\s+([a-zA-Z0-9/_]+)\s+command$")

parser = argparse.ArgumentParser(description='Fixup headers in docs')
parser.add_argument('files', metavar='FILE', nargs='+', help='files to fix')
args = parser.parse_args()

for orig_file in args.files:
    command_type = "unknown"
    styles = []
    headings = {}
    new_file = f"{orig_file}.tmp"
    found_syntax = False

    with open(orig_file, 'r') as reader, open(new_file, 'w') as writer:
        for line in reader:
            if line.startswith("Syntax"):
                found_syntax = True
                break

            m = index_pattern.match(line)

            if not m:
                m = pattern.match(line)

            if m:
                command_type = m.group(1)
                style        = m.group(2)

                if style not in styles:
                    styles.append(style)

                base_name = '/'.join(style.split('/')[0:-1])
                ext = style.split('/')[-1]

                if ext not in ('omp', 'intel', 'kk', 'gpu', 'opt'):
                    if style not in headings:
                        headings[style] = []
                elif style not in headings[base_name]:
                    headings[base_name].append(style)

        if found_syntax:
            # write new header
            for s in styles:
                print(f".. index:: {command_type} {s}", file=writer)

            print(file=writer)

            for s, variants in headings.items():
                header = f"{command_type} {s} command"
                print(header, file=writer)
                print("="*len(header), file=writer)
                print(file=writer)

                if len(variants) > 0:
                    print("Accelerator Variants: ", end="", file=writer)
                    print(", ".join([f"*{v}*" for v in variants]), file=writer)
                    print(file=writer)

            # write rest of reader
            print(line, end="", file=writer)
            for line in reader:
                print(line, end="", file=writer)

    if found_syntax and len(styles) > 0:
        # override original file
        shutil.move(new_file, orig_file)
    else:
        os.remove(new_file)
