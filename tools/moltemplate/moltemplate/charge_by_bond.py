#!/usr/bin/env python
# Author: Andrew Jewett (jewett.aij at g mail)
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2017, California Institute of Technology
# All rights reserved.

g_program_name = __file__.split('/')[-1]   # = 'charge_by_bond.py'
g_date_str = '2017-10-03'
g_version_str = '0.14.0'


import sys
import re
from collections import defaultdict

try:
    from . import ttree_lex
    from .lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid
except (ImportError, SystemError, ValueError):
    # not installed as a package
    import ttree_lex
    from lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid



def LookupChargePairs(chargebyatomid,
                      # bond_ids,
                      # bond_types,
                      # bond_pairs,
                      lines_atoms,
                      lines_bonds,
                      lines_bond_list,
                      lines_chargebybond,
                      atom_style,
                      section_name,
                      prefix='',
                      suffix=''):
                      # bond_ids_offset=0):
                      # report_progress = False):
    """
    LookupChargePairs() looks up partial-charge pair contributions from the
    types of atoms participating in a bond.

    Output:
    ...It looks up the corresponding change in the partial charges for
       each pair of atoms and stores this in the "chargebyatomid" dictionary.

    Input (continued):
       This function requires:
    ...a list of bonded pairs of atoms
        stored in the lines_bonds variable (from the "Data Bond List"
                                       or "Data Bonds AtomId AtomId" sections)
    ...and a list of atom types
        stored in the lines_atoms variable (from the "Data Atoms" section)

    ...and a list of charge-pairs-as-a-function-of-atom-types
        stored in the lines_chargebybond (from the "Data Bonds By Type" section)

    """

    column_names = AtomStyle2ColNames(atom_style)
    i_atomid, i_atomtype, i_molid = ColNames2AidAtypeMolid(column_names)

    atomids = []
    atomtypes = []
    atomids2types = {}

    for iv in range(0, len(lines_atoms)):
        line = lines_atoms[iv].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) > 0:
            tokens = ttree_lex.SplitQuotedString(line)
            if ((len(tokens) <= i_atomid) or (len(tokens) <= i_atomtype)):
                sys.stderr.write("\"" + line + "\"\n")
                raise(ttree_lex.InputError(
                    'Error not enough columns on line ' + str(iv + 1) + ' of \"Atoms\" section.'))
            tokens = ttree_lex.SplitQuotedString(line)
            atomid = ttree_lex.EscCharStrToChar(tokens[i_atomid])
            atomids.append(atomid)
            atomtype = ttree_lex.EscCharStrToChar(tokens[i_atomtype])
            atomtypes.append(atomtype)
            atomids2types[atomid] = atomtype

    #assert(isinstance(bond_ids, list))
    #assert(isinstance(bond_types, list))
    #assert(isinstance(bond_pairs, list))
    #del bond_ids[:]
    #del bond_types[:]
    #del bond_pairs[:]
    bond_pairs = []

    for ie in range(0, len(lines_bond_list)):
        line = lines_bond_list[ie].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) == 0:
            continue
        tokens = ttree_lex.SplitQuotedString(line)
        if len(tokens) == 3:
            # bond_ids.append(ttree_lex.EscCharStrToChar(tokens[0]))
            bond_pairs.append((ttree_lex.EscCharStrToChar(tokens[1]),
                               ttree_lex.EscCharStrToChar(tokens[2])))
        else:
            raise(ttree_lex.InputError('Incorrect number of columns on line ' +
                                       str(ie + 1) + ' of \"' + section_name + '\" section.'))

    for ie in range(0, len(lines_bonds)):
        line = lines_bonds[ie].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) == 0:
            continue
        tokens = ttree_lex.SplitQuotedString(line)
        if len(tokens) == 4:
            # bond_ids.append(ttree_lex.EscCharStrToChar(tokens[0]))
            # bond_types.append(ttree_lex.EscCharStrToChar(tokens[1]))
            bond_pairs.append((ttree_lex.EscCharStrToChar(tokens[2]),
                               ttree_lex.EscCharStrToChar(tokens[3])))
        else:
            raise(ttree_lex.InputError('Incorrect number of columns on line ' +
                                       str(ie + 1) + ' of \"' + section_name + '\" section.'))

    # for ie in range(0, len(lines_bonds_atomid_atomid)):
    #    line = lines_bonds_atomid_atomid[ie].strip()
    #    if '#' in line:
    #        icomment = line.find('#')
    #        line = (line[:icomment]).strip()
    #    if len(line) == 0:
    #        continue
    #    tokens = ttree_lex.SplitQuotedString(line)
    #    if len(tokens) == 2:
    #        #bondid_n = bond_ids_offset + len(bond_ids) + 1
    #        #bond_ids.append(prefix+str(bondid_n)+suffix)
    #        bond_pairs.append( (ttree_lex.EscCharStrToChar(tokens[0]),
    #                            ttree_lex.EscCharStrToChar(tokens[1])) )
    #    else:
    #        raise(ttree_lex.InputError('Incorrect number of columns on line '+str(ie+1)+' of \"'+section_name+'\" section.'))

    #assert(len(bond_types) == 0)
    typepattern_to_chargepairs = []
    warning_unassigned_chargepairs = None

    for i in range(0, len(lines_chargebybond)):
        line = lines_chargebybond[i].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) > 0:
            tokens = ttree_lex.SplitQuotedString(line)

            if (len(tokens) != 4):
                raise(ttree_lex.InputError('Error: Wrong number of columns in the \"Charge Pairs By Type\" section of data file.\n'
                                           'Offending line:\n' +
                                           '\"' + line + '\"\n'
                                           'Expected 4 columns\n'))

            chargepair = (float(tokens[2]),
                          float(tokens[3]))

            typepattern = []

            for typestr in tokens[:2]:
                if ((len(typestr) >= 2) and
                        (typestr[0] == '/') and (typestr[-1] == '/')):
                    regex_str = typestr[1:-1]
                    typepattern.append(re.compile(regex_str))
                else:
                    typepattern.append(ttree_lex.EscCharStrToChar(typestr))

            typepattern_to_chargepairs.append([typepattern, chargepair])

    for atomid1, atomid2 in bond_pairs:

        if atomid1 not in atomids2types:
            raise ttree_lex.InputError('Error: atom \"' + atomid1 + '\" not defined in \"Data Atoms\".\n'
                                       '       This usually happens when the user mistypes one of the names of the\n'
                                       '       $atoms in either a \"Data Atoms\" or \"Data Bond List\" section.\n'
                                       '       To find out where the mistake occured, search the \n'
                                       '       \"ttree_assignments.txt\" file for:\n'
                                       '       \"' + atomid1 + '\"\n')

        if atomid2 not in atomids2types:
            raise ttree_lex.InputError('Error: atom \"' + atomid2 + '\" not defined in \"Data Atoms\".\n'
                                       '       This usually happens when the user mistypes one of the names of the\n'
                                       '       $atoms in either a \"Data Atoms\" or \"Data Bond List\" section.\n'
                                       '       To find out where the mistake occured, search the \n'
                                       '       \"ttree_assignments.txt\" file for:\n'
                                       '       \"' + atomid2 + '\"\n')

        atomtype1 = atomids2types[atomid1]
        atomtype2 = atomids2types[atomid2]

        found = False
        for typepattern, chargepair in typepattern_to_chargepairs:
            # use string comparisons to check if atom types match the pattern
            if ttree_lex.MatchesAll((atomtype1, atomtype2), typepattern):
                # ("MatchesAll()" defined in "ttree_lex.py")
                chargebyatomid[atomid1] += chargepair[0]
                chargebyatomid[atomid2] += chargepair[1]
                found = True
            elif ttree_lex.MatchesAll((atomtype2, atomtype1), typepattern):
                chargebyatomid[atomid1] += chargepair[1]
                chargebyatomid[atomid2] += chargepair[0]
                found = True
        if (not found) and (not warning_unassigned_chargepairs):
            warning_unassigned_chargepairs = (atomid1, atomid2)

    if warning_unassigned_chargepairs:
        sys.stderr.write('---------------------------------------------------------------------------\n'
                         'Warning: bonds found between atoms with no partial-charge rules.\n'
                         '         This means that somewhere you are using a force-field\n'
                         '         which assigns atomic charge according to the bonds these atoms\n'
                         '         participate in, AND at least one pair of bonded atoms does NOT have\n'
                         '         a rule defined to assign charges to that pair of atoms.\n'
                         '         This can happen if there is a problem with the force-field file\n'
                         '         OR if you are defining the charges for these atoms manually\n'
                         '         In the later case, it is not a problem.\n'
                         '         The first bond with this problem is between this pair of atoms:\n'
                         '            ' +
                         str(warning_unassigned_chargepairs[0]) + '\n'
                         '            ' +
                         str(warning_unassigned_chargepairs[1]) + '\n'
                         '---------------------------------------------------------------------------\n')

def main():
    """
    This is is a "main module" wrapper for invoking chargepairs_by_type.py
    as a stand alone program.  This program:
    This program reads a LAMMPS data file (or an excerpt of a LAMMPS)
    data file containing bonded many-body interactions by atom type
    (and bond type), and generates a list of atom charges in LAMMPS input
    script format consistent with those types (to the standard out).

    Typical Usage:

    chargepairs_by_type.py -atoms atoms.data \\
                           -bonds bonds.data \\
                           -chargebybond chargepairs_by_type.data \\
                           > list_of_atom_charges.in

    """

    #######  Main Code Below: #######
    sys.stderr.write(g_program_name + ' v' +
                     g_version_str + ' ' + g_date_str + ' ')
    if sys.version < '3':
        sys.stderr.write(' (python version < 3)\n')
    else:
        sys.stderr.write('\n')

    try:
        fname_atoms = None
        fname_bonds = None
        fname_bond_list = None
        fname_chargebybond = None
        section_name = 'Data Bond List'  # (This will be replaced later.)
        atom_style = 'full'
        prefix = ''
        suffix = ''
        bond_lack_types = False

        argv = [arg for arg in sys.argv]

        # Loop over the remaining arguments not processed yet.
        # These arguments are specific to the lttree.py program
        # and are not understood by ttree.py:
        i = 1
        while i < len(argv):
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
            if ((argv[i].lower() == '-?') or
                    (argv[i].lower() == '--?') or
                    (argv[i].lower() == '-help') or
                    (argv[i].lower() == '-help')):
                if i + 1 >= len(argv):
                    sys.stdout.write(man_page_text + '\n')
                    sys.exit(0)

            elif argv[i].lower() == '-atoms':
                if i + 1 >= len(argv):
                    raise ttree_lex.InputError('Error: ' + argv[i] + ' flag should be followed by a file name containing lines of\n'
                                               '       text which might appear in the "Atoms" section of a LAMMPS data file.\n')
                fname_atoms = argv[i + 1]
                del(argv[i:i + 2])

            elif argv[i].lower() == '-bonds':
                if i + 1 >= len(argv):
                    raise ttree_lex.InputError('Error: ' + argv[i] + ' flag should be followed by a file name containing lines of\n'
                                               '       text which might appear in the "Bonds" section of a LAMMPS data file.\n')
                fname_bonds = argv[i + 1]
                del(argv[i:i + 2])

            elif argv[i].lower() == '-bond-list':
                if i + 1 >= len(argv):
                    raise ttree_lex.InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')
                    # raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                    #                 '       text which might appear in the "Bonds No Types" section of a LAMMPS data file.\n')
                fname_bond_list = argv[i + 1]
                section_name = "Data Bond List"
                del(argv[i:i + 2])

            elif ((argv[i].lower() == '-chargebybond') or
                  (argv[i].lower() == '-chargesbybond') or
                  (argv[i].lower() == '-charge-by-bond') or
                  (argv[i].lower() == '-charges-by-bond') or
                  (argv[i].lower() == '-chargepairsbytype') or
                  (argv[i].lower() == '-chargepairs-by-type') or
                  (argv[i].lower() == '-charge-pairs-by-type')):
                if i + 1 >= len(argv):
                    raise ttree_lex.InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')

                    # raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name containing\n'
                    #                 '       text which might appear in the "'+section_name+' By Type" section\n'
                    #                 '       of a LAMMPS data file.\n')
                fname_chargebybond = argv[i + 1]
                del(argv[i:i + 2])

            elif ((argv[i].lower() == '-atom-style') or
                  (argv[i].lower() == '-atom_style')):
                if i + 1 >= len(argv):
                    raise ttree_lex.InputError('Error: ' + argv[i] + ' flag should be followed by a an atom_style name.\n'
                                               '       (Or single quoted string which includes a space-separated\n'
                                               '       list of column names.)\n')
                atom_style = argv[i + 1]
                del(argv[i:i + 2])

            elif argv[i][0] == '-':
                raise ttree_lex.InputError('Error(' + g_program_name + '):\n'
                                           'Unrecogized command line argument \"' + argv[i] + '\"\n')
            else:
                i += 1

        if len(argv) != 1:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"' + arg + '\"' for arg in argv[1:]]
            raise ttree_lex.InputError('Syntax Error(' + g_program_name + '):\n\n'
                                       '       Problem with argument list.\n'
                                       '       The remaining arguments are:\n\n'
                                       '         ' +
                                       (' '.join(problem_args)) + '\n\n'
                                       '       (The actual problem may be earlier in the argument list.)\n')

        #bond_types = []
        fatoms = open(fname_atoms, 'r')
        lines_bonds = []
        lines_bond_list = []
        fbonds = fbond_list = None
        try:
            if fname_bonds != None:
                fbonds = open(fname_bonds, 'r')
                lines_bonds = fbonds.readlines()
                fbonds.close()
        except IOError:
            pass
        try:
            if fname_bond_list != None:
                fbond_list = open(fname_bond_list, 'r')
                lines_bond_list = fbond_list.readlines()
                fbond_list.close()
        except IOError:
            pass
        if ((len(lines_bonds) == 0) and (len(lines_bond_list) == 0)):
            sys.stderr.write('Error(' + g_program_name + '): No bonds defined for this system\n'
                             '      (This error may be a bug in moltemplate.)\n')
        fchargebybond = open(fname_chargebybond, 'r')
        lines_atoms = fatoms.readlines()

        lines_chargebybond = fchargebybond.readlines()
        fatoms.close()
        fchargebybond.close()
        chargebyatomid = defaultdict(float)

        LookupChargePairs(chargebyatomid,
                          lines_atoms,
                          lines_bonds,
                          lines_bond_list,
                          lines_chargebybond,
                          atom_style,
                          section_name)

        for atomid, charge in chargebyatomid.items():
            sys.stdout.write('  set atom ' + str(atomid) +
                             ' charge ' + str(charge) + '\n')

    except (ValueError, ttree_lex.InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)



if __name__ == "__main__":
    main()

