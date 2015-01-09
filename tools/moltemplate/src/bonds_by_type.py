#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2011, Regents of the University of California
# All rights reserved.

"""
    bonds_by_type.py reads a LAMMPS data file (or an excerpt of a LAMMPS)
    data file containing bonded many-body interactions by atom type
    (and bond type), and generates a list of additional interactions
    in LAMMPS format consistent with those type (to the standard out).

    Typical Usage:

    bonds_by_type.py -atoms atoms.data \\
                     -bonds bonds.data \\
                     -bondsbytype bonds_by_type.data \\
                     > new_bonds.data

"""

#                     -bonds-ids-atom-pairs bonds_ids_atom_pairs.data \\

import sys
#from extract_lammps_data import *
#from nbody_by_type_lib import GenInteractions_str
import ttree_lex
#from ttree_lex import *
from lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid



def LookupBondTypes(bond_types,
                    bond_ids,
                    bond_pairs,
                    lines_atoms,
                    lines_bonds,
                    lines_bondsbytype,
                    atom_style,
                    section_name,
                    prefix='',
                    suffix='',
                    bond_ids_offset=0):
                    #report_progress = False):
    """
    LookupBondTypes() looks up bond types.

    Output:
    ...It looks up the corresponding type of each bond and store it in the 
       "bond_types" list.  (If the bond_ids were not specified by the user, 
       generate them and store them in the bond_ids list.)


    Input (continued):
       This function requires:
    ...a list of bonded pairs of atoms
        stored in the lines_bonds variable (from the "Data Bond List"
                                       or "Data Bonds AtomId AtomId" sections)
    ...and a list of atom types
        stored in the lines_atoms variable (from the "Data Atoms" section)
    ...and a list of bond-types-as-a-function-of-atom-types
        stored in the lines_bondsbytype (from the "Data Bonds By Type" section)

    Generated bond_ids (if applicable) are of the form 
      prefix + str(number) + suffix
        (where "number" begins at bond_ids_offset+1)

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
                sys.stderr.write("\""+line+"\"\n")
                raise(ttree_lex.InputError('Error not enough columns on line '+str(iv+1)+' of \"Atoms\" section.'))
            tokens = ttree_lex.SplitQuotedString(line)
            atomid = ttree_lex.EscCharStrToChar(tokens[i_atomid])
            atomids.append(atomid)
            atomtype = ttree_lex.EscCharStrToChar(tokens[i_atomtype])
            atomtypes.append(atomtype)
            atomids2types[atomid] = atomtype


    assert(isinstance(bond_ids, list))
    assert(isinstance(bond_types, list))
    assert(isinstance(bond_pairs, list))
    del bond_ids[:]
    del bond_types[:]
    del bond_pairs[:]

    for ie in range(0, len(lines_bonds)):

        line = lines_bonds[ie].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()

        if len(line) == 0:
            continue

        tokens = ttree_lex.SplitQuotedString(line)

        if section_name == "Data Bonds AtomId AtomId":
            if len(tokens) == 2:
                bondid_n = bond_ids_offset + len(bond_ids) + 1
                bond_ids.append(prefix+str(bondid_n)+suffix)
                bond_pairs.append( (ttree_lex.EscCharStrToChar(tokens[0]),
                                        ttree_lex.EscCharStrToChar(tokens[1])) )
            else:
                raise(ttree_lex.InputError('Incorrect number of columns on line '+str(ie+1)+' of \"'+section_name+'\" section.'))

        elif section_name == "Data Bond List":
            if len(tokens) == 3:
                bond_ids.append(ttree_lex.EscCharStrToChar(tokens[0]))
                bond_pairs.append( (ttree_lex.EscCharStrToChar(tokens[1]),
                                    ttree_lex.EscCharStrToChar(tokens[2])) )
            else:
                raise(ttree_lex.InputError('Incorrect number of columns on line '+str(ie+1)+' of \"'+section_name+'\" section.'))

        #elif section_name == "Data Bonds BondType AtomId AtomId":
        #    if len(tokens) == 3:
        #        bondid_n = bond_ids_offset + len(bond_ids) + 1
        #        bond_ids.append(prefix+str(bondid_n)+suffix)
        #        bond_types.append(ttree_lex.EscCharStrToChar(tokens[0]))
        #        bond_pairs.append( (ttree_lex.EscCharStrToChar(tokens[1]),
        #                            ttree_lex.EscCharStrToChar(tokens[2])) )
        #    else:
        #        raise(ttree_lex.InputError('Incorrect number of columns on line '+str(ie+1)+' of \"'+section_name+'\" section.'))

        else:
            raise(ttree_lex.InputError('Internal Error ('+g_program_name+'): Unknown section name: \"'+section_name+'\"'))


    assert(len(bond_types) == 0)
    typepattern_to_coefftypes = []

    for i in range(0, len(lines_bondsbytype)):
        line = lines_bondsbytype[i].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) > 0:
            tokens = ttree_lex.SplitQuotedString(line)

            if (len(tokens) != 3):
                raise(ttree_lex.InputError('Error: Wrong number of columns in the \"Bonds By Type\" section of data file.\n'
                                 'Offending line:\n'+
                                 '\"'+line+'\"\n'
                                 'Expected 3 columns\n'))

            coefftype = ttree_lex.EscCharStrToChar(tokens[0])
            typepattern = []

            for typestr in tokens[1:]:
                if ((len(typestr) >= 2) and 
                    (typestr[0] == '/') and (typestr[-1] == '/')):
                    regex_str = typestr[1:-1]
                    typepattern.append( re.compile(regex_str) )
                else:
                    typepattern.append(ttree_lex.EscCharStrToChar(typestr))

            typepattern_to_coefftypes.append([typepattern, coefftype])
            


    assert(len(bond_ids) == len(bond_pairs))

    for ie in range(0,len(bond_ids)):
        bond_types.append(None)

    for ie in range(0, len(bond_ids)):
        bondid = bond_ids[ie]
        (atomid1, atomid2) = bond_pairs[ie]

        #for n in range(0, len(typepattern_to_coefftypes)):
        for typepattern, coefftype in typepattern_to_coefftypes:

            if atomid1 not in atomids2types:
                raise ttree_lex.InputError('Error: atom \"'+atomid1+'\" not defined in \"Data Atoms\".\n'
                                           '       This usually happens when the user mistypes one of the names of the\n'
                                           '       $atoms in either a \"Data Atoms\" or \"Data Bond List\" section.\n'
                                           '       To find out where the mistake occured, search the \n'
                                           '       \"ttree_assignments.txt\" file for:\n'
                                           '       \"'+atomid1+'\"\n')

            if atomid2 not in atomids2types:
                raise ttree_lex.InputError('Error: atom \"'+atomid2+'\" not defined in \"Data Atoms\".\n'
                                           '       This usually happens when the user mistypes one of the names of the\n'
                                           '       $atoms in either a \"Data Atoms\" or \"Data Bond List\" section.\n'
                                           '       To find out where the mistake occured, search the \n'
                                           '       \"ttree_assignments.txt\" file for:\n'
                                           '       \"'+atomid2+'\"\n')

            atomtype1 = atomids2types[atomid1]
            atomtype2 = atomids2types[atomid2]

            # use string comparisons to check if atom types match the pattern
            if (ttree_lex.MatchesAll((atomtype1, atomtype2), typepattern) or
                ttree_lex.MatchesAll((atomtype2, atomtype1), typepattern)):
                # ("MatchesAll()" defined in "ttree_lex.py")

                bond_types[ie] = coefftype

    for ie in range(0, len(bond_ids)):
        if not bond_types[ie]:
            atomtype1 = atomids2types[atomid1]
            atomtype2 = atomids2types[atomid2]
            raise ttree_lex.InputError('Error: No bond types defined for the bond between\n'
                              '       atoms '+atomid1+' (type '+atomtype1+')\n'
                              '         and '+atomid2+' (type '+atomtype2+')\n')




if __name__ == "__main__":

    g_program_name = __file__.split('/')[-1]  # = 'nbody_by_type.py'
    g_date_str     = '2013-8-06'
    g_version_str  = '0.1'

    #######  Main Code Below: #######
    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')
    if sys.version < '3':
        sys.stderr.write(' (python version < 3)\n')
    else:
        sys.stderr.write('\n')

    try:
        fname_atoms = None
        fname_bonds = None
        fname_bondsbytype = None
        section_name = 'Data Bond List'  # (This will be replaced later.)
        atom_style = 'full'
        prefix=''
        suffix=''
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
                if i+1 >= len(argv):
                    sys.stdout.write(man_page_text+'\n')
                    sys.exit(0)

            elif argv[i].lower() == '-atoms':
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                                     '       text which might appear in the "Atoms" section of a LAMMPS data file.\n')
                fname_atoms = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-bonds':
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                                     '       text which might appear in the "Bonds" section of a LAMMPS data file.\n')
                fname_bonds = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-bond-list':
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name\n')
                    #raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                    #                 '       text which might appear in the "Bonds No Types" section of a LAMMPS data file.\n')
                fname_bonds = argv[i+1]
                section_name = "Data Bond List"
                del(argv[i:i+2])

            elif argv[i].lower() == '-bondsbytype':
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name\n')

                    #raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a file name containing\n'
                    #                 '       text which might appear in the "'+section_name+' By Type" section\n'
                    #                 '       of a LAMMPS data file.\n')
                fname_bondsbytype = argv[i+1]
                del(argv[i:i+2])

            elif ((argv[i].lower() == '-atom-style') or 
                (argv[i].lower() == '-atom_style')):
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a an atom_style name.\n'
                                     '       (Or single quoted string which includes a space-separated\n'
                                     '       list of column names.)\n')
                atom_style = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-prefix':
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a prefix string\n'
                                     '       (a string you want to appear to the left of the integer\n'
                                     '        which counts the bonded interactions you have generated.)\n')
                prefix = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-suffix':
                if i+1 >= len(argv):
                    raise ttree_lex.InputError('Error: '+argv[i]+' flag should be followed by a suffix string\n'
                                     '       (a string you want to appear to the right of the integer\n'
                                     '        which counts the bonded interactions you have generated.)\n')
                prefix = argv[i+1]
                del(argv[i:i+2])

            elif argv[i][0] == '-':
                raise ttree_lex.InputError('Error('+g_program_name+'):\n'
                                 'Unrecogized command line argument \"'+argv[i]+'\"\n')
            else:
                i += 1

        if len(argv) != 1:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"'+arg+'\"' for arg in argv[1:]]
            raise ttree_lex.InputError('Syntax Error('+g_program_name+'):\n\n'
                             '       Problem with argument list.\n'
                             '       The remaining arguments are:\n\n'
                             '         '+(' '.join(problem_args))+'\n\n'
                             '       (The actual problem may be earlier in the argument list.)\n')

        bond_types = []
        bond_ids = []
        bond_pairs = []

        fatoms = open(fname_atoms, 'r')
        fbonds = open(fname_bonds, 'r')
        fbondsbytype = open(fname_bondsbytype, 'r')
        lines_atoms = fatoms.readlines()
        lines_bonds = fbonds.readlines()
        lines_bondsbytype = fbondsbytype.readlines()
        fatoms.close()
        fbonds.close()
        fbondsbytype.close()

        LookupBondTypes(bond_types,
                        bond_ids,
                        bond_pairs,
                        lines_atoms,
                        lines_bonds,
                        lines_bondsbytype,
                        atom_style,
                        section_name,
                        prefix='',
                        suffix='')

        assert(len(bond_types) == len(bond_ids) == len(bond_pairs))

        ie=0
        N = len(bond_types)
        for ie in range(0, N):
            sys.stdout.write(bond_ids[ie] + ' ' +
                             bond_types[ie] + ' ' +
                             bond_pairs[ie][0] + ' ' +
                             bond_pairs[ie][1] + '\n')


    except (ValueError, ttree_lex.InputError) as err:
        sys.stderr.write('\n'+str(err)+'\n')
        sys.exit(-1)

