#!/usr/bin/env python

"""
   Reorder the atoms in the Angles section of a data file to make sure that 
   atoms have a "canonical order" (for example the first atom has a lower
   id than the last atom, for angle and dihedral interactions.
   (This helps us detect potential problems like dupicate Angle interactions.)

"""

import sys
from operator import itemgetter

g_program_name = __file__.split('/')[-1]

in_stream = sys.stdin


section_name = ''
if len(sys.argv) == 3:
    section_name = sys.argv[1]
    module_name = sys.argv[2].rstrip('.py')
else:
    sys.stderr.write('Usage Example:\n\n'
                     '      '+g_program_name+' Angles nbody_angles.py < angles.txt > new_angles.txt\n\n'
                     '      In this example \"angles.txt\" contains only the \"Angles\" section of\n'
                     '      a LAMMPS DATA file.  (Either a text-editor, or the \n'
                     '      \"extract_lammps_data.py\" script can be used to select a section from\n'
                     '      a LAMMPS DATA file\n\n'
                     'Error('+g_program_name+'): expected exactly one argument:\n'
                     '       \"Angles\",  \"Dihedrals\", or \"Impropers\"\n')
    exit(-1)

# Ordering rules are defined in a seperate module named
# nbody_angles.py, nbody_dihedrals.py, nbody_impropers.py
# Load that now.

g = __import__(module_name)  #defines g.bond_pattern, g.canonical_order

# This module defines the graph representing the bond pattern for this type
# of interaction.  (The number of vertices and edges for the graph corresponds
# to the number of atoms and bonds in this type of interaction.)
natoms = g.bond_pattern.GetNumVerts()
nbonds = g.bond_pattern.GetNumEdges()


for line_orig in in_stream:
    line = line_orig.rstrip('\n')
    comment = ''
    if '#' in line_orig:
        ic = line.find('#')
        line = line_orig[:ic]
        comment = ' '+line_orig[ic:].rstrip('\n')

    tokens = line.strip().split()
    swapped = False
    if len(tokens) == 2+natoms:
        all_integers = True
        abids_l = [[0 for i in range(0, natoms)],
                   [0 for i in range(0, nbonds)]]
        for i in range(0, natoms):
            if not tokens[2+i].isdigit():
                all_integers = False
        if all_integers:
            for i in range(0, natoms):
                abids_l[0][i] = int(tokens[2+i])
        else:
            for i in range(0, natoms):
                abids_l[0][i] = tokens[2+i]
        abids = g.canonical_order( (tuple(abids_l[0]), tuple(abids_l[1])) )
        for i in range(0, natoms):
            tokens[2+i] = str(abids[0][i])

    sys.stdout.write(' '.join(tokens)+comment+'\n')
