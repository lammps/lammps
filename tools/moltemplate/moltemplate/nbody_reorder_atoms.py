#!/usr/bin/env python

"""
   Reorder the atoms in the Angles section of a data file to make sure that
   atoms have a "canonical order" (for example the first atom has a lower
   id than the last atom, for angle and dihedral interactions.
   (This helps us detect potential problems like dupicate Angle interactions.)

"""

from operator import itemgetter
import importlib
import os
import sys
sys.path.append(os.getcwd())

g_program_name = __file__.split('/')[-1]


def main():
    in_stream = sys.stdin

    section_name = ''
    if len(sys.argv) == 3:
        section_name = sys.argv[1]
        bond_pattern_module_name = sys.argv[2]
        # If the file name ends in ".py", then strip off this suffix.
        # The next line does not work. Too lazy do care why.
        # bond_pattern_module_name=bond_pattern_module_name.rstrip('.py')
        # Do this instead
        pc = bond_pattern_module_name.rfind('.py')
        if pc != -1:
            bond_pattern_module_name = bond_pattern_module_name[0:pc]

    else:
        sys.stderr.write('Usage Example:\n\n'
                         '      ' + g_program_name + ' Angles nbody_angles.py < angles.txt > new_angles.txt\n\n'
                         '      In this example \"angles.txt\" contains only the \"Angles\" section of\n'
                         '      a LAMMPS DATA file.  (Either a text-editor, or the \n'
                         '      \"extract_lammps_data.py\" script can be used to select a section from\n'
                         '      a LAMMPS DATA file\n\n'
                         'Error(' + g_program_name +
                         '): expected exactly one argument:\n'
                         '       \"Angles\",  \"Dihedrals\", or \"Impropers\"\n')
        exit(-1)

    # Ordering rules are defined in a seperate module named
    # nbody_angles.py, nbody_dihedrals.py, nbody_impropers.py
    # Load that now.

    # search locations
    package_opts = [[bond_pattern_module_name, __package__],
                    ['nbody_alt_symmetry.'+bond_pattern_module_name,
                     __package__]]

    if __package__:
        for i, _ in enumerate(package_opts):
            package_opts[i][0] = '.' + package_opts[i][0]

    g = None
    for name, pkg in package_opts:
        try:
            # defines g.bond_pattern, g.canonical_order
            g = importlib.import_module(name, pkg)
            break
        except (ImportError, SystemError, ValueError):
            pass

    if g is None:
        sys.stderr.write('Error: Unable to locate file \"' +
                         bond_pattern_module_name + '\"\n'
                         '       (Did you mispell the file name?\n'
                         '        Check the \"nbody_alternate_symmetry/\" directory.)\n')
        sys.exit(-1)

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
            comment = ' ' + line_orig[ic:].rstrip('\n')

        tokens = line.strip().split()
        swapped = False
        if len(tokens) == 2 + natoms:
            all_integers = True
            abids_l = [[0 for i in range(0, natoms)],
                       [0 for i in range(0, nbonds)]]
            for i in range(0, natoms):
                if not tokens[2 + i].isdigit():
                    all_integers = False
            if all_integers:
                for i in range(0, natoms):
                    abids_l[0][i] = int(tokens[2 + i])
            else:
                for i in range(0, natoms):
                    abids_l[0][i] = tokens[2 + i]
            abids = g.canonical_order((tuple(abids_l[0]), tuple(abids_l[1])))
            for i in range(0, natoms):
                tokens[2 + i] = str(abids[0][i])

        sys.stdout.write(' '.join(tokens) + comment + '\n')

    return


if __name__ == '__main__':
    main()
