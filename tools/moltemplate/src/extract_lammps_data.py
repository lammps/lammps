#!/usr/bin/env python

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
                            'MiddleBondTorsion Coeffs', # class2 dihedrals
                            'EndBondTorsion Coeffs',    # class2 dihedrals
                            'AngleTorsion Coeffs',      # class2 dihedrals
                            'AngleAngleTorsion Coeffs', # class2 dihedrals
                            'BondBond13 Coeffs',        # class2 dihedrals
                            'AngleAngle Coeffs',        # class2 impropers
                            'Angles By Type',   # new. not standard LAMMPS
                            'Dihedrals By Type',# new. not standard LAMMPS
                            'Angles By Type'])   # new. not standard LAMMPS

def ExtractDataSection(f,
                       section_header, 
                       comment_char = '#',
                       include_header = False,
                       return_line_nums = False):
    inside_section = False
    nonblank_encountered = False
    i = 0
    for line_orig in f:
        return_this_line = False
        line = line_orig.strip()
        if line == section_header:
            inside_section = True
            nonblank_encountered = False
            if include_header:
                return_this_line = True
        # A block of blank lines (which dont immediately follow
        # the section header-name) signal the end of a section:
        elif len(line) == 0:
            if inside_section and include_header:
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


def FindDataSection(f,
                    section_header, 
                    comment_char = '#'):
    i_section_start = -1
    i_section_stop = -1
    inside_section = False
    nonblank_encountered = False
    i = 0
    for line_orig in f:
        line = line_orig.strip()
        if line == section_header:
            inside_section = True
            nonblank_encountered = False
        # A block of blank lines (which dont immediately follow
        # the section header-name) signal the end of a section:
        elif len(line) == 0:
            if nonblank_encountered:
                inside_section = False
                i_section_stop = i
                break
        elif line[0] != comment_char:
            if inside_section:
                if not nonblank_encountered:
                    i_section_start = i # <- first non-blank line
                nonblank_encountered = True
        i += 1

    if i_section_stop == -1:
        if i_section_start != -1:
            i_section_stop = i

    return (i_section_start, i_section_stop)


if __name__ == "__main__":

    import sys
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
                                               include_header=True,
                                               return_line_nums=True):
                line_nums_exclude.add(line_num)
        for i in range(0, len(lines)):
            if i not in line_nums_exclude:
                sys.stdout.write(lines[i])
