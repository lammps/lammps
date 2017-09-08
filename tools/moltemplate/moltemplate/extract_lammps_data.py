#!/usr/bin/env python
"""
extract_lammps_data.py is a simple script which extracts sections of text from
a LAMMPS data file.

Typical usage: 

extract_lammps_data.py SECTION_NAME < FILE.DATA > EXCERPT.TXT

This extracts a section from a LAMMPS data file and saves it in EXCERPT.TXT.

More general usage:

extract_lammps_data.py [-n] SECTION_LIST < FILE.DATA > EXCERPT.TXT

For more details, see "doc/utils/docs_extract_lammps_data.txt"
"""

import sys

lammps_data_sections = set(['Atoms',
                            'Masses',
                            'Bonds',
                            'Bond Coeffs',
                            'Angles',
                            'Angle Coeffs',
                            'Dihedrals',
                            'Dihedral Coeffs',
                            'Impropers',
                            'Improper Coeffs',
                            'BondBond Coeffs',          # class2 angles
                            'BondAngle Coeffs',         # class2 angles
                            'MiddleBondTorsion Coeffs',  # class2 dihedrals
                            'EndBondTorsion Coeffs',    # class2 dihedrals
                            'AngleTorsion Coeffs',      # class2 dihedrals
                            'AngleAngleTorsion Coeffs',  # class2 dihedrals
                            'BondBond13 Coeffs',        # class2 dihedrals
                            'AngleAngle Coeffs',        # class2 impropers
                            'Angles By Type',   # new. not standard LAMMPS
                            'Dihedrals By Type',  # new. not standard LAMMPS
                            'Angles By Type'])   # new. not standard LAMMPS


def DeleteComments(string,
                   escape='\\',
                   comment_char='#'):
    escaped_state = False
    for i in range(0, len(string)):
        if string[i] in escape:
            if escaped_state:
                escaped_state = False
            else:
                escaped_state = True
        elif string[i] == comment_char:
            if not escaped_state:
                return string[0:i]
    return string


def ExtractDataSection(f,
                       section_name,
                       comment_char='#',
                       include_section_name=False,
                       return_line_nums=False):

    inside_section = False
    if section_name in ('header', 'Header'):  # "Header" section includes beginning
        inside_section = True

    nonblank_encountered = False
    nonheader_encountered = False

    i = 0
    for line_orig in f:
        return_this_line = False
        line = DeleteComments(line_orig).strip()
        if line in lammps_data_sections:
            nonheader_encountered = True
        if section_name in ('header', 'Header'):
            # The "header" section includes all lines at the beginning of the
            # before any other section is encountered.
            if nonheader_encountered:
                return_this_line = False
            else:
                return_this_line = True
        elif line == section_name:
            inside_section = True
            nonblank_encountered = False
            if include_section_name:
                return_this_line = True
        # A block of blank lines (which dont immediately follow
        # the section_name) signal the end of a section:
        elif len(line) == 0:
            if inside_section and include_section_name:
                return_this_line = True
            if nonblank_encountered:
                inside_section = False
        elif line[0] != comment_char:
            if inside_section:
                nonblank_encountered = True
                return_this_line = True

        if return_this_line:
            if return_line_nums:
                yield i
            else:
                yield line_orig

        i += 1


def main():
    lines = sys.stdin.readlines()
    exclude_sections = False
    if sys.argv[1] == '-n':
        exclude_sections = True
        del sys.argv[1]

    if not exclude_sections:
        for section_name in sys.argv[1:]:
            for line in ExtractDataSection(lines, section_name):
                sys.stdout.write(line)
    else:
        line_nums_exclude = set([])
        for section_name in sys.argv[1:]:
            for line_num in ExtractDataSection(lines,
                                               section_name,
                                               include_section_name=True,
                                               return_line_nums=True):
                line_nums_exclude.add(line_num)
        for i in range(0, len(lines)):
            if i not in line_nums_exclude:
                sys.stdout.write(lines[i])

    return

if __name__ == "__main__":
    main()
