#!/usr/bin/env python

"""
   Get rid of lines containing duplicate bonded nbody interactions in the 
   corresponding section of a LAMMPS data file (such as bonds, angles,
   dihedrals and impropers).  Duplicate lines which occur later are
   preserved and the earlier lines are erased.
   (This program reads from sys.stdin.  This program does not parse the entire 
    data file.  The text from the relevant section of the LAMMPS file should be 
   extracted in advance before it is sent to this program.)

"""

import sys
in_stream = sys.stdin

if len(sys.argv) == 2:
    n = int(sys.argv[1])
if (len(sys.argv) != 2) or (n < 1):
    sys.stderr.write('Error (remove_duplicates_nbody.py): expected a positive integer argument.\n')
    sys.exit(-1)

atom_ids_in_use = set([])

lines = in_stream.readlines()

# Start at the end of the file and read backwards.
# If duplicate lines exist, eliminate the ones that occur earlier in the file.
i = len(lines)
while i > 0:
    i -= 1
    line_orig = lines[i]
    line = line_orig.rstrip('\n')
    if '#' in line_orig:
        ic = line.find('#')
        line = line_orig[:ic]

    tokens = line.strip().split()
    if len(tokens) == 2+n:
        atom_ids = tuple(tokens[2:2+n])
        if atom_ids in atom_ids_in_use:
            del lines[i]
        else:
            atom_ids_in_use.add(atom_ids)
    elif len(tokens) == 0:
        del lines[i]

for line in lines:
    sys.stdout.write(line)
