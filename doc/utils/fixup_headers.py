#!/usr/bin/env python3
# Utility script to fix headers in doc pages and generate "Accelerator Styles" portion
import shutil
import re
import argparse

pattern = re.compile(r"^(fix|pair_style)\s+([a-zA-Z0-9/_]+)\s+command$")

parser = argparse.ArgumentParser(description='Fixup headers in docs')
parser.add_argument('files', metavar='FILE', nargs='+', help='files to fix')
args = parser.parse_args()

for orig_file in args.files:
    command_type = "unknown"
    styles = []
    headings = {}
    new_file = f"{orig_file}.tmp"

    with open(orig_file, 'r') as reader, open(new_file, 'w') as writer:
        for line in reader:
            if line.startswith("Syntax"):
                break

            m = pattern.match(line)
            if m:
                command_type = m.group(1)
                style        = m.group(2)

                styles.append(style)

                base_name = '/'.join(style.split('/')[0:-1])
                ext = style.split('/')[-1]

                if ext not in ('omp', 'intel', 'kk', 'gpu', 'opt'):
                    headings[style] = []
                else:
                    headings[base_name].append(style)
                
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
                print("Accelerator Styles: ", end="", file=writer)
                print(", ".join([f"*{v}*" for v in variants]), file=writer)

            print(file=writer)

        # write rest of reader
        print(line, end="", file=writer)
        for line in reader:
            print(line, end="", file=writer)

    # override original file
    shutil.move(new_file, orig_file)
