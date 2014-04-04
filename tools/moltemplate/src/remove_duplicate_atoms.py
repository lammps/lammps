#!/usr/bin/env python

"""
   Get rid of lines containing duplicate copies of the same atom in the "Atoms"
   section of a LAMMPS data file.  Duplicate lines which occur later are
   preserved and the earlier lines are erased.
   The file is read from sys.stdin.  This program does not parse the entire 
   data file.  The text from the "Atoms" section of the LAMMPS file must 
   be extracted in advance before it is sent to this program.)

"""

import sys
in_stream = sys.stdin
f = None
fname = None
if len(sys.argv) == 2:
    fname = sys.argv[1]
    f = open(fname, 'r')
    in_stream = f

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
    if len(tokens) > 0:
        atom_id = tokens[0]
        if atom_id in atom_ids_in_use:
            del lines[i]
        else:
            atom_ids_in_use.add(atom_id)
    else:
        del lines[i]


for line in lines:
    sys.stdout.write(line)

if f != None:
    f.close()

