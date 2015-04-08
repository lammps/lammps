#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2012, Regents of the University of California
# All rights reserved.

"""
ltemplify.py

The "ltemplify.py" script can be used to convert existing LAMMPS
input script and data files into a single .lt file 
(which includes both topology and force-field information 
 for a single molecule in your system).

Example:

   ltemplify.py -name Mol file.in file.data > mol.lt

This creates a template for a new type of molecule (named "Mol"),
consisting of all the atoms in the lammps files you included,
and saves this data in a single ttree file ("mol.lt").
This file can be used with moltemplate (ttree) to
define large systems containing this molecule.

"""

import sys
from ttree_lex import *
from lttree_styles import *



def Intify(s):
    if s.isdigit():
        return int(s)
    elif s[0:2] == 'id':
        return int(s[2:])
    elif s[0:4] == 'type':
        return int(s[4:])
    else:
        return s


def StringToInterval(sel_str, slice_delim='*'):
    # Split a string into 1-3 tokens using the slice_delim and convert to int.
    # What a mess. I should rewrite this function

    i_slice = sel_str.find(slice_delim)

    if i_slice == -1:
        a = sel_str
        b = sel_str
        c = ''
    else:
        a  = sel_str[:i_slice]
        bc = sel_str[i_slice+len(slice_delim):]
        b = ''
        c = ''
        i_slice = bc.find(slice_delim)
        if i_slice == -1:
            b = bc
            c = ''
        else:
            b = bc[:i_slice]
            c = bc[i_slice+len(slice_delim):]

    if a == '':
        a = None
    elif a.isdigit():
        a = int(a)
    else:
        raise InputError('Error: invalid selection string \"'+
                         sel_str+'\"\n')

    if b == '':
        b = None
    elif b.isdigit():
        b = int(b)
    else:
        raise InputError('Error: invalid selection string \"'+
                         sel_str+'\"\n')

    if c == '':
        c = None
    elif c.isdigit():
        c = int(c)
    else:
        raise InputError('Error: invalid selection string \"'+
                         sel_str+'\"\n')

    if c == None:
        return (a,b)
    else:
        return (a,b,c)



# Selections are simply lists of 2-tuples (pairs) 

def LammpsSelectToIntervals(sel_str, slice_delim='*', or_delim=', '):

    """
    This function converts a string such as "1*4 6 9*12 50*70*10" into 
    a list of tuples, for example: [(1,4), (6,6), (9,12), (50,50), (60,60), (70,70)]
    In general, the of intervals has the form:
       [(a1,b1), (a2,b2), (a3,b3), ... ]
    
    An atom is considered to belong to this selection 
    if it happens to lie within the closed interval [a,b] 
    for any pair of a,b values in the list of intervals. 
    If for a given pair a,b, either a or b is "None", then that a or b 
    value is not used to disqualify membership in the interval. 
    (Similar to -infinity or +infinity.  In other words if a is set to None, 
     then to belong to the interval it is enough to be less than b.)

    """
    selection_list = []
    #tokens = sel_str.split(or_delim) <-- Not what we want when len(or_delim)>1
    tokens = LineLex.TextBlock2Lines(sel_str, or_delim, keep_delim=False)
    for token in tokens:
        token = token.strip()
        interval = StringToInterval(token, slice_delim)

        if len(interval)==2:
            # Normally, "interval" should be a tuple containing 2 entries
            selection_list.append(interval)
        else:
            assert(len(interval)==3)
            # Handle 1000:2000:10 notation
            # (corresponding to 1000, 1010, 1020, 1030, ..., 1990, 2000)
            a=interval[0]
            b=interval[1]
            incr=interval[2]
            i=a
            while i<=b:
                selection_list.append((i,i))
                i += incr

    return selection_list


def IntervalListToMinMax(interval_list):
    min_a = None
    max_b = None
    for (a,b) in interval_list:
        if ((not (type(a) is int)) or (not (type(b) is int))):
            return None,None #only integer min/max makes sense. otherwise skip

        if (min_a == None) or (a < min_a):
            min_a = a
        if (max_b == None) or (b > max_b):
            max_b = b
    return min_a, max_b


def MergeIntervals(interval_list):
    """
    A crude simple function that merges consecutive intervals in the list
    whenever they overlap.  (This function does not bother to compare 
    non-consecutive entries in the interval_list.)

    """
    i = 1
    while i < len(interval_list):
        if ((interval_list[i-1][1] == None) or
            (interval_list[i-1][1]+1 >= interval_list[i][0])):
            interval_list[i-1] = (interval_list[i-1][0], interval_list[i][1])
            del interval_list[i]
        else:
            i += 1


def BelongsToSel(i, sel):
    if (i == None) or (sel == None) or (len(sel) == 0):
        # If the user has not specified a selection for this category,
        # then by default all objects are accepted
        return True

    elif (type(i) is str):
        if i.isdigit():
            i = int(i)
        else:
            return True

    belongs = False
    for interval in sel:
        assert(len(interval) == 2)
        if interval[0]:
            if i >= interval[0]:
                if (interval[1] == None) or (i <= interval[1]):
                    belongs = True
                    break
        elif interval[1]:
            if i <= interval[1]:
                belongs = True
                break
        else:
            # In that case, the user entered something like "*"
            # which covers all possible numbers
            belongs = True
            break

    return belongs



try:

    g_program_name = __file__.split('/')[-1]  # = 'ltemplify.py'
    g_version_str  = '0.45'
    g_date_str     = '2014-12-17'
    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

    non_empty_output = False
    no_warnings = True
    indent = 2
    cindent = 0
    atomid_selection = []
    atomtype_selection = []
    molid_selection = []
    mol_name = ''

    min_sel_atomid = None
    min_sel_atomtype = None
    min_sel_bondid = None
    min_sel_bondtype = None
    min_sel_angleid = None
    min_sel_angletype = None
    min_sel_dihedralid = None
    min_sel_dihedraltype = None
    min_sel_improperid = None
    min_sel_impropertype = None

    max_sel_atomid = None
    max_sel_atomtype = None
    max_sel_bondid = None
    max_sel_bondtype = None
    max_sel_angleid = None
    max_sel_angletype = None
    max_sel_dihedralid = None
    max_sel_dihedraltype = None
    max_sel_improperid = None
    max_sel_impropertype = None

    needed_atomids = set([])
    needed_atomtypes = set([])
    needed_molids = set([])
    needed_bondids = set([])
    needed_bondtypes = set([])
    needed_angleids = set([])
    needed_angletypes = set([])
    needed_dihedralids = set([])
    needed_dihedraltypes = set([])
    needed_improperids = set([])
    needed_impropertypes = set([])

    min_needed_atomtype = None
    max_needed_atomtype = None
    min_needed_bondtype = None
    max_needed_bondtype = None
    min_needed_angletype = None
    max_needed_angletype = None
    min_needed_dihedraltype = None
    max_needed_dihedraltype = None
    min_needed_impropertype = None
    max_needed_impropertype = None

    min_needed_atomid = None
    max_needed_atomid = None
    min_needed_molid = None
    max_needed_molid = None
    min_needed_bondid = None
    max_needed_bondid = None
    min_needed_angleid = None
    max_needed_angleid = None
    min_needed_dihedralid = None
    max_needed_dihedralid = None
    min_needed_improperid = None
    max_needed_improperid = None


    # To process the selections, we need to know the atom style:
    atom_style_undefined = True

    i_atomid   = None
    i_atomtype = None
    i_molid    = None

    l_in_init     = []
    l_in_settings = []
    l_in_masses = []
    l_in_pair_coeffs = []
    l_in_bond_coeffs = []
    l_in_angle_coeffs = []
    l_in_dihedral_coeffs = []
    l_in_improper_coeffs = []
    l_in_group = []
    l_in_fix_shake = []
    l_in_fix_rigid = []
    l_in_fix_poems = []
    l_in_fix_qeq = []
    l_in_fix_qmmm = []
    l_data_masses = []
    l_data_bond_coeffs = []
    l_data_angle_coeffs = []
    l_data_dihedral_coeffs = []
    l_data_improper_coeffs = []
    l_data_pair_coeffs = []
    l_data_atoms = []
    l_data_velocities = []
    l_data_bonds = []
    l_data_angles = []
    l_data_dihedrals = []
    l_data_impropers = []

    # class2 force fields
    l_data_bondbond_coeffs = []
    l_data_bondangle_coeffs = []
    l_data_middlebondtorsion_coeffs = []
    l_data_endbondtorsion_coeffs = []
    l_data_angletorsion_coeffs = []
    l_data_angleangletorsion_coeffs = []
    l_data_bondbond13_coeffs = []
    l_data_angleangle_coeffs = []

    # non-point-like particles:
    l_data_ellipsoids = []
    l_data_lines = []
    l_data_triangles = []

    # automatic generation of bonded interactions by type:
    l_data_angles_by_type = []
    l_data_dihedrals_by_type = []
    l_data_impropers_by_type = []

    atoms_already_read = False
    some_pair_coeffs_read = False
    complained_atom_style_mismatch = False
    infer_types_from_comments = False


    argv = sys.argv

    i = 1

    while i < len(argv):

        #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')

        if argv[i] == '-columns':
            if i+1 >= len(argv):
                raise InputError('Error: the \"'+argv[i]+'\" argument should be followed by a quoted\n'
                                 '       string which contains a space-delimited list of the names of\n'
                                 '       of columns in the \"Atoms\" section of the LAMMPS data file.\n'
                                 '       If the list contains the symbols:\n'
                                 '    \"atom-ID\" or \"atomid\", they are interpreted\n'
                                 '       as unique atom ID numbers, and columns named\n'
                                 '    \"atom-type\" or \"atomtype\" are interpreted\n'
                                 '       as atom types.  Finally, columns named\n'
                                 '    \"molecule-ID\", \"molecule\", or \"mol-ID\", or \"mol\"\n'
                                 '       are interpreted as unique molecule id numbers.\n'
                                 'Example:\n'
                                 '    '+argv[i]+' \'atom-ID atom-type q polarizability molecule-ID x y z\'\n'
                                 '       defines a custom atom_style containing the properties\n'
                                 '            atom-ID atom-type q polarizability molecule-ID x y z\n'
                                 '    Make sure you enclose the entire list in quotes.\n');
            column_names = argv[i+1].strip('\"\'').strip().split()
            del argv[i:i+2]

        elif (argv[i] == '-ignore-comments'):
            infer_types_from_comments = False
            del argv[i:i+1]

        elif (argv[i] == '-infer-comments'):
            infer_types_from_comments = True
            del argv[i:i+1]

        elif ((argv[i] == '-name') or
              (argv[i] == '-molname') or
              (argv[i] == '-molecule-name') or
              (argv[i] == '-molecule_name')):
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a a molecule type name.\n')
            cindent = 2
            indent += cindent
            mol_name = argv[i+1]
            del argv[i:i+2]

        elif ((argv[i].lower() == '-atomstyle') or 
              (argv[i].lower() == '-atom_style') or 
              (argv[i].lower() == '-atom-style')):
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a an atom_style name.\n'
                                 '       (or single quoted string which includes a space-separated\n'
                                 '       list of column names).\n')
            atom_style_undefined = False
            column_names = AtomStyle2ColNames(argv[i+1])
            if (argv[i+1].strip().split()[0] in g_style_map):
                l_in_init.append((' '*indent) + 'atom_style ' + argv[i+1] + '\n')
            sys.stderr.write('\n    \"Atoms\" column format:\n')
            sys.stderr.write('    '+(' '.join(column_names))+'\n')
            i_atomid, i_atomtype, i_molid = ColNames2AidAtypeMolid(column_names)
            if i_molid:
                sys.stderr.write('      (i_atomid='+str(i_atomid+1)+', i_atomtype='+str(i_atomtype+1)+', i_molid='+str(i_molid+1)+')\n\n')
            else:
                sys.stderr.write('      (i_atomid='+str(i_atomid+1)+', i_atomtype='+str(i_atomtype+1)+')\n')
            del argv[i:i+2]

        elif ((argv[i].lower() == '-id') or 
              #(argv[i].lower() == '-a') or 
              #(argv[i].lower() == '-atoms') or 
              (argv[i].lower() == '-atomid') or 
              #(argv[i].lower() == '-atomids') or 
              (argv[i].lower() == '-atom-id')
              #(argv[i].lower() == '-atom-ids') or 
              #(argv[i].lower() == '-$atom') or
              #(argv[i].lower() == '-$atoms')
              ):
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a list of integers\n'
                                 '       (or strings).  These identify the group of atoms you want to\n'
                                 '       to include in the template you are creating.\n')
            atomid_selection += LammpsSelectToIntervals(argv[i+1])
            min_sel_atomid, max_sel_atomid = IntervalListToMinMax(atomid_selection)
            del argv[i:i+2]
        elif ((argv[i].lower() == '-type') or 
              #(argv[i].lower() == '-t') or 
              (argv[i].lower() == '-atomtype') or 
              (argv[i].lower() == '-atom-type')
              #(argv[i].lower() == '-atomtypes') or 
              #(argv[i].lower() == '-atom-types') or 
              #(argv[i].lower() == '-@atom') or 
              #(argv[i].lower() == '-@atoms') or 
              #(argv[i].lower() == '-@atomtype') or
              #(argv[i].lower() == '-@atomtypes')
              ):
            if i+1 >= len(argv):
                raise InputError('Error: '+argv[i]+' flag should be followed by a list of integers.\n'
                                 '       (or strings).  These identify the group of atom types you want to\n'
                                 '       to include in the template you are creating.\n')
            atomtype_selection += LammpsSelectToIntervals(argv[i+1])
            min_sel_atomtype, max_sel_atomtype = IntervalListToMinMax(atomtype_selection)
            del argv[i:i+2]
        elif ((argv[i].lower() == '-mol') or 
              #(argv[i].lower() == '-m') or 
              (argv[i].lower() == '-molid') or 
              #(argv[i].lower() == '-molids') or 
              (argv[i].lower() == '-mol-id') or 
              #(argv[i].lower() == '-mol-ids') or 
              #(argv[i].lower() == '-molecule') or 
              (argv[i].lower() == '-moleculeid') or 
              (argv[i].lower() == '-molecule-id')
              #(argv[i].lower() == '-molecules') or 
              #(argv[i].lower() == '-molecule-ids') or 
              #(argv[i].lower() == '-$mol') or
              #(argv[i].lower() == '-$molecule')
              ):
            if i+1 >= len(argv):
                sys.stderr.write('Error: '+argv[i]+' flag should be followed by a list of integers.\n'
                                 '       (or strings).  These identify the group of molecules you want to\n'
                                 '       include in the template you are creating.\n')
            molid_selection += LammpsSelectToIntervals(argv[i+1])
            del argv[i:i+2]
        else:
            i += 1


    # atom type names
    atomtypes_name2int = {}
    atomtypes_int2name = {}
    #atomids_name2int = {}  not needed
    atomids_int2name = {}
    atomids_by_type = {}


    if atom_style_undefined:
        # The default atom_style is "full"
        column_names = AtomStyle2ColNames('full')
        i_atomid, i_atomtype, i_molid = ColNames2AidAtypeMolid(column_names)

    #---------------------------------------------------------
    #-- The remaining arguments are files that the user wants
    #-- us to read and convert.  It is typical to have 
    #-- multiple input files, because LAMMPS users often
    #-- store their force field parameters in either the LAMMPS
    #-- data files and input script files, or both.
    #-- We want to search all of the LAMMPS input files in
    #-- order to make sure we extracted all the force field
    #-- parameters (coeff commands).
    #---------------------------------------------------------

    for i_arg in range(1,len(argv)):
        fname = argv[i_arg]
        try:
            lammps_file = open(fname, 'r')
        except IOError: 
            raise InputError('Error: unrecognized argument (\"'+fname+'\"),\n'
                             '       OR unable to open file:\n'
                             '\n'
                             '       \"'+fname+'\"\n'
                             '       for reading.\n'
                             '\n'
                             '       (If you were not trying to open a file with this name,\n'
                             '        then there is a problem in your argument list.)\n')

        sys.stderr.write('reading file \"'+fname+'\"\n')

        atomid2type = {}
        atomid2mol  = {}
        data_file_header_names = set(['LAMMPS Description',
                                      'Atoms', 'Masses', 'Velocities', 'Bonds',
                                      'Angles', 'Dihedrals', 'Impropers',
                                      'Pair Coeffs', 
                                      'Bond Coeffs', 'Angle Coeffs', 
                                      'Dihedral Coeffs', 'Improper Coeffs',
                                      #class2 force fields:
                                      'BondBond Coeffs', 'BondAngle Coeffs',
                                      'MiddleBondTorsion Coeffs', 'EndBondTorsion Coeffs',
                                      'AngleTorsion Coeffs', 'AngleAngleTorsion Coeffs',
                                      'BondBond13 Coeffs',
                                      'AngleAngle Coeffs',
                                      # non-point-like particles:
                                      'Ellipsoids', 'Triangles', 'Lines',
                                      #specifying bonded interactions by type:
                                      'Angles By Type', 'Dihedrals By Type', 'Impropers By Type'
                                      ])

        lex=LineLex(lammps_file, fname)
        lex.source_triggers = set(['include','import'])
        # set up lex to accept most characters in file names:
        lex.wordterminators = '(){}' + lex.whitespace
        # set up lex to understand the "include" statement:
        lex.source = 'include'
        lex.escape = '\\'

        while lex:
            infile = lex.infile
            lineno = lex.lineno
            line = lex.ReadLine()
            if (lex.infile != infile):
                infile = lex.infile
                lineno = lex.lineno

            #sys.stderr.write('  processing \"'+line.strip()+'\", (\"'+infile+'\":'+str(lineno)+')\n')

            if line == '':
                break

            tokens = line.strip().split()
            if (len(tokens) > 0):
                if ((tokens[0] == 'atom_style') and
                    atom_style_undefined):
                    
                    sys.stderr.write('  Atom Style found. Processing: \"'+line.strip()+'\"\n')
                    if atoms_already_read:
                        raise InputError('Error: The file containing the \"atom_style\" command must\n'
                                         '       come before the data file in the argument list.\n'
                                         '       (The templify program needs to know the atom style before reading\n'
                                         '       the data file.  Either change the order of arguments so that the\n'
                                         '       LAMMPS input script file is processed before the data file, or use\n'
                                         '       the \"-atom_style\" command line argument to specify the atom_style.)\n')

                    column_names = AtomStyle2ColNames(line.split()[1])
                    i_atomid, i_atomtype, i_molid = ColNames2AidAtypeMolid(column_names)

                    sys.stderr.write('\n    \"Atoms\" column format:\n')
                    sys.stderr.write('    '+(' '.join(column_names))+'\n')
                    if i_molid:
                        sys.stderr.write('        (i_atomid='+str(i_atomid+1)+', i_atomtype='+str(i_atomtype+1)+', i_molid='+str(i_molid+1)+')\n\n')
                    else:
                        sys.stderr.write('        (i_atomid='+str(i_atomid+1)+', i_atomtype='+str(i_atomtype+1)+')\n\n')
                    l_in_init.append((' '*indent)+line.lstrip())

                elif (tokens[0] in set(['units',
                                        'angle_style',
                                        'bond_style',
                                        'dihedral_style',
                                        'improper_style',
                                        'min_style',
                                        'pair_style',
                                        'pair_modify',
                                        'special_bonds',
                                        'kspace_style',
                                        'kspace_modify'])):
                    l_in_init.append((' '*indent)+line.lstrip())

                #if (line.strip() == 'LAMMPS Description'):
                #    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                #    # skip over this section
                #    while lex:
                #        line = lex.ReadLine()
                #        if line.strip() in data_file_header_names:
                #            lex.push_raw_text(line) # <- Save line for later
                #            break
                  
                elif (line.strip() == 'Atoms'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    atoms_already_read = True
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            if ((len(tokens) <= i_atomid) or
                                (len(tokens) <= i_atomtype) or
                                ((i_molid != None) and
                                 (len(tokens) <= i_molid))):
                                raise InputError('Error: The number of columns in the \"Atoms\" section does\n'
                                                 '       not match the atom_style (see column name list above).\n')
                            elif ((len(tokens) != len(column_names)) and 
                                (not complained_atom_style_mismatch)):
                                complained_atom_style_mismatch = True
                                sys.stderr.write('Warning: The number of columns in the \"Atoms\" section does\n'
                                                 '         not match the atom_style (see column name list above).\n')
                                # this is not a very serious warning.
                                #no_warnings = False <--no need. commenting out


                            atomid   = Intify(tokens[i_atomid])
                            atomtype = Intify(tokens[i_atomtype])

                            molid = None
                            if i_molid:
                                molid = Intify(tokens[i_molid])

                            atomid2type[atomid] = atomtype
                            if i_molid:
                                atomid2mol[atomid] = molid


                            if (BelongsToSel(atomid, atomid_selection) and
                                BelongsToSel(atomtype, atomtype_selection) and
                                BelongsToSel(molid, molid_selection)):

                                tokens[i_atomid] = '$atom:id'+tokens[i_atomid]
                                #tokens[i_atomid] = '$atom:'+atomids_int2name[atomid]
                                # fill atomtype_int2str[] with a default name (change later):
                                #tokens[i_atomtype] = '@atom:type'+tokens[i_atomtype]
                                atomtype_name = 'type'+tokens[i_atomtype]
                                atomtypes_int2name[atomtype] = atomtype_name
                                tokens[i_atomtype] = '@atom:'+atomtype_name

                                # I can't use atomids_int2names or atomtypes_int2names yet
                                # because they probably have not been defined yet.
                                # (Instead assign these names in a later pass.)

                                if i_molid:
                                    tokens[i_molid]    = '$mol:id'+tokens[i_molid]
                                l_data_atoms.append((' '*indent)+(' '.join(tokens)+'\n'))
                                needed_atomids.add(atomid)
                                needed_atomtypes.add(atomtype)
                                needed_molids.add(molid)

                    for atomtype in needed_atomtypes:
                        assert(type(atomtype) is int)
                        if ((min_needed_atomtype == None) or
                            (min_needed_atomtype > atomtype)):
                            min_needed_atomtype = atomtype
                        if ((max_needed_atomtype == None) or
                            (max_needed_atomtype < atomtype)):
                            max_needed_atomtype = atomtype

                    for atomid in needed_atomids:
                        assert(type(atomid) is int)
                        if ((min_needed_atomid == None) or
                            (min_needed_atomid > atomid)):
                            min_needed_atomid = atomid
                        if ((max_needed_atomid == None) or
                            (max_needed_atomid < atomid)):
                            max_needed_atomid = atomid

                    for molid in needed_molids:
                        assert(type(molid) is int)
                        if ((min_needed_molid == None) or
                            (min_needed_molid > molid)):
                            min_needed_molid = molid
                        if ((max_needed_molid == None) or
                            (max_needed_molid < molid)):
                            max_needed_molid = molid

                elif (line.strip() == 'Masses'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        # Read the next line of text but don't skip comments
                        comment_char_backup = lex.commenters
                        lex.commenters = ''
                        line_orig = lex.ReadLine()
                        lex.commenters = comment_char_backup

                        comment_text = ''
                        ic = line_orig.find('#')
                        line = line_orig[:ic]
                        if ic != -1:
                            comment_text = line_orig[ic+1:].strip()
                                
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break

                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            atomtype = Intify(tokens[0])
                            atomtype_name = str(atomtype)

                            if comment_text != '':
                                comment_tokens = comment_text.split()
                                # Assume the first word after the # is the atom type name
                                atomtype_name = comment_tokens[0]

                            if BelongsToSel(atomtype, atomtype_selection):
                                #tokens[0] = '@atom:type'+tokens[0]
                                l_data_masses.append((' '*indent)+(' '.join(tokens)+'\n'))
                                # infer atom type names from comment strings?
                                if infer_types_from_comments:
                                    if atomtype_name in atomtypes_name2int:
                                        raise InputError('Error: duplicate atom type names in mass section: \"'+atomtype_name+'\"\n'
                                                         '       (By default '+g_program_name+' attempts to infer atom type names from\n'
                                                         '       comments which appear in the \"Masses\" section of your data file.)\n'
                                                         '       You can avoid this error by adding the \"-ignore-comments\" argument.\n')
                                    atomtypes_name2int[atomtype_name] = atomtype
                                    atomtypes_int2name[atomtype] = atomtype_name
                                else:
                                    atomtypes_int2name[atomtype] = 'type'+str(atomtype)


                elif (line.strip() == 'Velocities'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            atomid  = Intify(tokens[0])
                            atomtype = None
                            if atomid in atomid2type:
                                atomtype = atomid2type[atomid]
                            moldid = None
                            if atomid in atomid2mol:
                                molid = atomid2mol[atomid]
                            if (BelongsToSel(atomid, atomid_selection) and
                                BelongsToSel(atomtype, atomtype_selection) and
                                BelongsToSel(molid, molid_selection)):
                                #tokens[0] = '$atom:id'+tokens[0]
                                tokens[0] = '$atom:'+atomids_int2name[atomid]
                                l_data_velocities.append((' '*indent)+(' '.join(tokens)+'\n'))

                # non-point-like-particles:
                elif (line.strip() == 'Ellipsoids'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            atomid  = Intify(tokens[0])
                            atomtype = None
                            if atomid in atomid2type:
                                atomtype = atomid2type[atomid]
                            moldid = None
                            if atomid in atomid2mol:
                                molid = atomid2mol[atomid]
                            if (BelongsToSel(atomid, atomid_selection) and
                                BelongsToSel(atomtype, atomtype_selection) and
                                BelongsToSel(molid, molid_selection)):
                                #tokens[0] = '$atom:id'+tokens[0]
                                tokens[0] = '$atom:'+atomids_int2name[atomid]
                                l_data_ellipsoids.append((' '*indent)+(' '.join(tokens)+'\n'))
                elif (line.strip() == 'Lines'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            atomid  = Intify(tokens[0])
                            atomtype = None
                            if atomid in atomid2type:
                                atomtype = atomid2type[atomid]
                            moldid = None
                            if atomid in atomid2mol:
                                molid = atomid2mol[atomid]
                            if (BelongsToSel(atomid, atomid_selection) and
                                BelongsToSel(atomtype, atomtype_selection) and
                                BelongsToSel(molid, molid_selection)):
                                #tokens[0] = '$atom:id'+tokens[0]
                                tokens[0] = '$atom:'+atomids_int2name[atomid]
                                l_data_lines.append((' '*indent)+(' '.join(tokens)+'\n'))
                elif (line.strip() == 'Triangles'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            atomid  = Intify(tokens[0])
                            atomtype = None
                            if atomid in atomid2type:
                                atomtype = atomid2type[atomid]
                            moldid = None
                            if atomid in atomid2mol:
                                molid = atomid2mol[atomid]
                            if (BelongsToSel(atomid, atomid_selection) and
                                BelongsToSel(atomtype, atomtype_selection) and
                                BelongsToSel(molid, molid_selection)):
                                #tokens[0] = '$atom:id'+tokens[0]
                                tokens[0] = '$atom:'+atomids_int2name[atomid]
                                l_data_triangles.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Bonds'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            if (len(tokens) < 4):
                                raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                                 '       Nonsensical line in Bonds section:\n'
                                                 '       \"'+line.strip()+'\"\n')
                            #tokens[0] = '$bond:id'+tokens[0]
                            #tokens[1] = '@bond:type'+tokens[1]
                            atomids = [None, None]
                            atomtypes = [None, None]
                            molids = [None, None]
                            in_selections = True
                            some_in_selection = False
                            for n in range(0,2):
                                atomids[n]  = Intify(tokens[2+n])
                                if atomids[n] in atomid2type:
                                    atomtypes[n] = atomid2type[atomids[n]]
                                if atomids[n] in atomid2mol:
                                    molids[n] = atomid2mol[atomids[n]]
                                if (BelongsToSel(atomids[n], atomid_selection) and
                                    BelongsToSel(atomtypes[n], atomtype_selection) and
                                    BelongsToSel(molids[n], molid_selection)):
                                    #tokens[2+n] = '$atom:id'+tokens[2+n]
                                    #tokens[2+n] = '$atom:'+atomids_int2name[atomids[n]]
                                    some_in_selection = True
                                else:
                                    in_selections = False
                            if in_selections:
                                l_data_bonds.append((' '*indent)+(' '.join(tokens)+'\n'))
                            elif some_in_selection:
                                sys.stderr.write('WARNING: SELECTION BREAKS BONDS\n')
                                sys.stderr.write('         (between atom ids: ')

                                for n in range(0,2):
                                    sys.stderr.write(str(atomids[n])+' ')
                                sys.stderr.write(')\n'
                                                 '         The atoms you selected are bonded\n'
                                                 '         to other atoms you didn\'t select.\n'
                                                 '         Are you sure you selected the correct atoms?\n')
                                no_warnings = False



                elif (line.strip() == 'Angles'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line == '':
                            break
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            if (len(tokens) < 5):
                                raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                                 '       Nonsensical line in Angles section:\n'
                                                 '       \"'+line.strip()+'\"\n')
                            #tokens[0] = '$angle:id'+tokens[0]
                            #tokens[1] = '@angle:type'+tokens[1]
                            atomids = [None, None, None]
                            atomtypes = [None, None, None]
                            molids = [None, None, None]
                            in_selections = True
                            some_in_selection = False
                            for n in range(0,3):
                                atomids[n]  = Intify(tokens[2+n])
                                if atomids[n] in atomid2type:
                                    atomtypes[n] = atomid2type[atomids[n]]
                                if atomids[n] in atomid2mol:
                                    molids[n] = atomid2mol[atomids[n]]
                                if (BelongsToSel(atomids[n], atomid_selection) and
                                    BelongsToSel(atomtypes[n], atomtype_selection) and
                                    BelongsToSel(molids[n], molid_selection)):
                                    #tokens[2+n] = '$atom:id'+tokens[2+n]
                                    #tokens[2+n] = '$atom:'+atomids_int2name[atomids[n]]
                                    some_in_selection = True
                                else:
                                    in_selections = False
                            if in_selections:
                                l_data_angles.append((' '*indent)+(' '.join(tokens)+'\n'))
                            elif some_in_selection:
                                sys.stderr.write('WARNING: SELECTION BREAKS ANGLES\n')
                                sys.stderr.write('         (between atom ids: ')
                                for n in range(0,3):
                                    sys.stderr.write(str(atomids[n])+' ')
                                sys.stderr.write(')\n'
                                                 '         The atoms you selected participate in 3-body \"Angle\"\n'
                                                 '         interactions with other atoms you didn\'t select.\n'
                                                 '         (They will be ignored.)\n'
                                                 '         Are you sure you selected the correct atoms?\n')
                                no_warnings = False


                elif (line.strip() == 'Dihedrals'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            if (len(tokens) < 6):
                                raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                                 '       Nonsensical line in Dihedrals section:\n'
                                                 '       \"'+line.strip()+'\"\n')
                            #tokens[0] = '$dihedral:id'+tokens[0]
                            #tokens[1] = '@dihedral:type'+tokens[1]
                            atomids = [None, None, None, None]
                            atomtypes = [None, None, None, None]
                            molids = [None, None, None, None]
                            in_selections = True
                            some_in_selection = False
                            for n in range(0,4):
                                atomids[n]  = Intify(tokens[2+n])
                                if atomids[n] in atomid2type:
                                    atomtypes[n] = atomid2type[atomids[n]]
                                if atomids[n] in atomid2mol:
                                    molids[n] = atomid2mol[atomids[n]]
                                if (BelongsToSel(atomids[n], atomid_selection) and
                                    BelongsToSel(atomtypes[n], atomtype_selection) and
                                    BelongsToSel(molids[n], molid_selection)):
                                    #tokens[2+n] = '$atom:id'+tokens[2+n]
                                    #tokens[2+n] = '$atom:'+atomids_int2name[atomids[n]]
                                    some_in_selection = True
                                else:
                                    in_selections = False
                            if in_selections:
                                l_data_dihedrals.append((' '*indent)+(' '.join(tokens)+'\n'))
                            elif some_in_selection:
                                sys.stderr.write('WARNING: SELECTION BREAKS DIHEDRALS\n')
                                sys.stderr.write('         (between atom ids: ')
                                for n in range(0,4):
                                    sys.stderr.write(str(atomids[n])+' ')
                                sys.stderr.write(')\n'
                                                 '         The atoms you selected participate in 4-body \"Dihedral\"\n'
                                                 '         interactions with other atoms you didn\'t select.\n'
                                                 '         (They will be ignored.)\n'
                                                 '         Are you sure you selected the correct atoms?\n')
                                no_warnings = False


                elif (line.strip() == 'Impropers'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            if (len(tokens) < 6):
                                raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                                 '       Nonsensical line in Impropers section:\n'
                                                 '       \"'+line.strip()+'\"\n')
                            #tokens[0] = '$improper:id'+tokens[0]
                            #tokens[1] = '@improper:type'+tokens[1]
                            atomids = [None, None, None, None]
                            atomtypes = [None, None, None, None]
                            molids = [None, None, None, None]
                            in_selections = True
                            some_in_selection = False
                            for n in range(0,4):
                                atomids[n]  = Intify(tokens[2+n])
                                if atomids[n] in atomid2type:
                                    atomtypes[n] = atomid2type[atomids[n]]
                                if atomids[n] in atomid2mol:
                                    molids[n] = atomid2mol[atomids[n]]
                                if (BelongsToSel(atomids[n], atomid_selection) and
                                    BelongsToSel(atomtypes[n], atomtype_selection) and
                                    BelongsToSel(molids[n], molid_selection)):
                                    #tokens[2+n] = '$atom:id'+tokens[2+n]
                                    #tokens[2+n] = '$atom:'+atomids_int2name[atomids[n]]
                                    some_in_selection = True
                                else:
                                    in_selections = False
                            if in_selections:
                                l_data_impropers.append((' '*indent)+(' '.join(tokens)+'\n'))
                            elif some_in_selection:
                                sys.stderr.write('WARNING: SELECTION BREAKS IMPROPERS\n')
                                sys.stderr.write('         (between atom ids: ')
                                for n in range(0,4):
                                    sys.stderr.write(str(atomids[n])+' ')
                                sys.stderr.write(')\n'
                                                 '         The atoms you selected participate in 4-body \"Improper\"\n'
                                                 '         interactions with other atoms you didn\'t select.\n'
                                                 '         (They will be ignored.)\n'
                                                 '         Are you sure you selected the correct atoms?\n')
                                no_warnings = False


                elif (line.strip() == 'Bond Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            #tokens[0] = '@bond:type'+tokens[0]
                            l_data_bond_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Angle Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            #tokens[0] = '@angle:type'+tokens[0]
                            l_data_angle_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Dihedral Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            #tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_dihedral_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Improper Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            #tokens[0] = '@improper:type'+tokens[0]
                            l_data_improper_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Pair Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    some_pair_coeffs_read = True
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            if (len(tokens) < 2):
                                raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                                 '       Nonsensical line in Pair Coeffs section:\n'
                                                 '       \"'+line.strip()+'\"\n')
                            atomtype_i_str = tokens[0]
                            if '*' in atomtype_i_str:
                                raise InputError('PROBLEM near or before '+ErrorLeader(infile, lineno)+'\n'
                                                 '         As of 2012-7, moltemplate forbids use of the "\*\" wildcard\n'
                                                 '         character in the \"Pair Coeffs\" section.\n')
                            else:
                                i = int(atomtype_i_str)
                                if ((not i) or 
                                    BelongsToSel(i, atomtype_selection)):
                                    i_str = '@atom:type'+str(i)
                                    tokens[0] = i_str
                                    l_data_pair_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (tokens[0] == 'pair_coeff'):
                    some_pair_coeffs_read = True
                    if (len(tokens) < 3):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical pair_coeff command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    l_in_pair_coeffs.append(' '*indent+line.strip())

                elif (tokens[0] == 'mass'):
                    some_pair_coeffs_read = True
                    if (len(tokens) < 3):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical \"mass\" command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    l_in_masses.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (tokens[0] == 'bond_coeff'):
                    if (len(tokens) < 2):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical bond_coeff command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    #tokens[1] = '@bond:type'+tokens[1]
                    l_in_bond_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (tokens[0] == 'angle_coeff'):
                    if (len(tokens) < 2):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical angle_coeff command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    #tokens[1] = '@angle:type'+tokens[1]
                    l_in_angle_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (tokens[0] == 'dihedral_coeff'):
                    if (len(tokens) < 2):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical dihedral_coeff command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    #tokens[1] = '@dihedral:type'+tokens[1]
                    l_in_dihedral_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))
                elif (tokens[0] == 'improper_coeff'):
                    if (len(tokens) < 2):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical improper_coeff command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    #tokens[1] = '@improper:type'+tokens[1]
                    l_in_improper_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))


                # -- class2 force fields --
                elif (line.strip() == 'BondBond Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@angle:type'+tokens[0]
                            l_data_bondbond_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'BondAngle Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@angle:type'+tokens[0]
                            l_data_bondangle_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'MiddleBondTorsion Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_middlebondtorsion_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'EndBondTorsion Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_endbondtorsion_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'AngleTorsion Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_angletorsion_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'AngleAngleTorsion Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_angleangletorsion_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'BondBond13 Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_bondbond13_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'AngleAngle Coeffs'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@improper:type'+tokens[0]
                            l_data_angleangle_coeffs.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Angles By Type'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@angle:type'+tokens[0]
                            l_data_angles_by_type.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Dihedrals By Type'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@dihedral:type'+tokens[0]
                            l_data_dihedrals_by_type.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (line.strip() == 'Impropers By Type'):
                    sys.stderr.write('  reading \"'+line.strip()+'\"\n')
                    while lex:
                        line = lex.ReadLine()
                        if line.strip() in data_file_header_names:
                            lex.push_raw_text(line) # <- Save line for later
                            break
                        tokens = line.strip().split()
                        if len(tokens) > 0:
                            tokens[0] = '@improper:type'+tokens[0]
                            l_data_impropers_by_type.append((' '*indent)+(' '.join(tokens)+'\n'))

                elif (tokens[0] == 'group'):
                    if (len(tokens) < 3):
                        raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                         '       Nonsensical group command:\n'
                                         '       \"'+line.strip()+'\"\n')
                    l_in_group.append((' '*indent)+(' '.join(tokens)+'\n'))
                
                elif ((tokens[0] == 'fix') and (len(tokens) >= 4)):
                    if (tokens[3].find('rigid') == 0):
                        if (len(tokens) < 6):
                            raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                             '       Nonsensical '+tokens[0]+' '+tokens[3]+' command:\n'
                                             '       \"'+line.strip()+'\"\n')
                        l_in_fix_rigid.append((' '*indent)+(' '.join(tokens)+'\n'))
                    elif (tokens[3].find('shake') == 0):
                        if (len(tokens) < 7):
                            raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                             '       Nonsensical '+tokens[0]+' '+tokens[3]+' command:\n'
                                             '       \"'+line.strip()+'\"\n')
                        l_in_fix_shake.append((' '*indent)+(' '.join(tokens)+'\n'))
                    elif (tokens[3].find('poems') == 0):
                        if (len(tokens) < 4):
                            raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                             '       Nonsensical '+tokens[0]+' '+tokens[3]+' command:\n'
                                             '       \"'+line.strip()+'\"\n')
                        l_in_fix_poems.append((' '*indent)+(' '.join(tokens)+'\n'))
                    elif (tokens[3].find('qeq') == 0):
                        if (len(tokens) < 8):
                            raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                             '       Nonsensical '+tokens[0]+' '+tokens[3]+' command:\n'
                                             '       \"'+line.strip()+'\"\n')
                        l_in_fix_qeq.append((' '*indent)+(' '.join(tokens)+'\n'))
                    elif (tokens[3].find('qmmm') == 0):
                        if (len(tokens) < 8):
                            raise InputError('Error: near or before '+ErrorLeader(infile, lineno)+'\n'
                                             '       Nonsensical '+tokens[0]+' '+tokens[3]+' command:\n'
                                             '       \"'+line.strip()+'\"\n')
                        l_in_fix_qmmm.append((' '*indent)+(' '.join(tokens)+'\n'))
                    elif (tokens[3].find('restrain') == 0):
                        sys.stderr('WARNING: fix \"'+tokens[3]+'\" commands are NOT understood by '+g_program_name+'.\n'
                                   '  If you need restraints, add them to your final .LT file (eg. \"system.lt\"),\n'
                                   '  (And be sure to use unique (full, long) moltemplate names for each $atom:.)\n'
                                   '  Ignoring line \"'+line.strip()+'\"\n')

                else:
                    sys.stderr.write('  Ignoring line \"'+line.strip()+'\"\n')

    sys.stderr.write('\n\n')

    sys.stderr.write('  processing \"Atoms\" section (')

    # post-processing:

    if len(l_data_masses) == 0:
        infer_types_from_comments = False

    # Pass 1 through l_data_atoms:
    # Now do a second-pass throught the "l_data_atoms" section, and 
    # finish dealing with "infer_types_from_comments".
    # During this pass, peplace the atomtype names and atomid names with 
    # atom type names which were inferred from comments read earlier.

    sys.stderr.write('pass1')
    for i in range(0, len(l_data_atoms)):
        tokens = l_data_atoms[i].split()
        atomid = tokens[i_atomid]
        if atomid.find('$atom:') == 0:
            atomid = atomid[6:]
            # convert to an integer
            atomid = Intify(atomid)

        if infer_types_from_comments:
            atomtype = tokens[i_atomtype]
            # remove the "@atom:" prefix (we will put it back later)
            if atomtype.find('@atom:') == 0:
                atomtype = atomtype[6:]
            # convert to an integer
            atomtype = Intify(atomtype)
            atomtype_name = atomtypes_int2name[atomtype]
            if atomtype in atomids_by_type:
                l_atomids = atomids_by_type[atomtype]
                prev_count = len(l_atomids)
                # lookup the most recently added atom of this type:
                #prev_atomid_name = l_atomids[-1]
                #ic = prev_atomid_name.rfind('_')
                #prev_count = int(prev_atomid_name[ic+1:])
                atomid_name = atomtype_name+'_'+str(prev_count+1)
                atomids_by_type[atomtype].append(atomid)
            else:
                atomids_by_type[atomtype] = [atomid]
                atomid_name = atomtype_name+'_1'
            atomids_int2name[atomid] = atomid_name
            #atomids_name2str[atomid_name] = atomid
        else:
            atomids_int2name[atomid] = 'id'+str(atomid)

    sys.stderr.write(', pass2')
    # Pass 2: If any atom types only appear once, simplify their atomid names.
    for i in range(0, len(l_data_atoms)):
        tokens = l_data_atoms[i].split()

        # remove the "@atom:" prefix (we will put it back later)
        atomtype = tokens[i_atomtype]
        if atomtype.find('@atom:') == 0:
            atomtype = atomtype[6:]
        atomtype = Intify(atomtype)
        if infer_types_from_comments:
            if len(atomids_by_type[atomtype]) == 1:
                atomid = tokens[i_atomid]
                if atomid.find('$atom:') == 0:
                    atomid = atomid[6:]
                atomid = Intify(atomid)
                atomtype_name = atomtypes_int2name[atomtype]
                atomids_int2name[atomid] = atomtype_name

    sys.stderr.write(', pass3')
    # Pass 3: substitute the atomid names and atom type names into l_data_atoms
    for i in range(0, len(l_data_atoms)):
        tokens = l_data_atoms[i].split()
        atomid = tokens[i_atomid]
        if atomid.find('$atom:') == 0:
            atomid = atomid[6:]
            # convert to an integer
            atomid = Intify(atomid)
        atomtype = tokens[i_atomtype]
        if atomtype.find('@atom:') == 0:
            atomtype = atomtype[6:]
        atomtype = Intify(atomtype)
        tokens = l_data_atoms[i].split()
        tokens[i_atomid] = '$atom:'+atomids_int2name[atomid]
        tokens[i_atomtype] = '@atom:'+atomtypes_int2name[atomtype]
        l_data_atoms[i] = (' '*indent)+(' '.join(tokens)+'\n')
    sys.stderr.write(')\n')


    if len(l_data_atoms) == 0:
        raise InputError('Error('+g_program_name+'): You have no atoms in you selection!\n'
                         '\n'
                         '       Either you have chosen a set of atoms, molecules, or atom types which\n'
                         '       does not exist, or there is a problem with (the format of) your\n'
                         '       arguments. Check the documentation and examples.\n')


    # --- Now delete items that were not selected from the other lists ---

    # --- MASSES ---

    # delete masses for atom types we don't care about anymore:
    i_line = 0
    while i_line < len(l_data_masses):
        line = l_data_masses[i_line]
        tokens = line.strip().split()
        atomtype = Intify(tokens[0])
        if ((not (atomtype in needed_atomtypes)) and
            (not ((len(atomtype_selection) > 0) and
                  BelongsToSel(atomtype, atomtype_selection)))):
            del l_data_masses[i_line]
        else:
            atomtype_name = atomtypes_int2name[atomtype]
            tokens[0] = '@atom:'+atomtype_name
            l_data_masses[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1



    # --- PAIR COEFFS ---

    # delete data_pair_coeffs for atom types we don't care about anymore:
    i_line = 0
    while i_line < len(l_data_pair_coeffs):
        line = l_data_pair_coeffs[i_line]
        tokens = line.strip().split()
        assert(len(tokens) > 0)
        split_colon = tokens[0].split(':')
        assert(len(split_colon) == 2)
        atomtype = Intify(split_colon[1])
        if ((not (atomtype in needed_atomtypes)) and
            (not ((len(atomtype_selection) > 0) and
                  BelongsToSel(atomtype, atomtype_selection)))):
            del l_data_pair_coeffs[i_line]
        else:
            i_line += 1

    # delete in_pair_coeffs for atom we don't care about anymore:
    i_line = 0
    while i_line < len(l_in_pair_coeffs):
        line = l_in_pair_coeffs[i_line]
        tokens = line.strip().split()
        atomtype_i_str = tokens[1]
        atomtype_j_str = tokens[2]
        #if (('*' in atomtype_i_str) or
        #    ('*' in atomtype_j_str)):
        #    sys.stderr.write('WARNING: near or before '+ErrorLeader(infile, lineno)+'\n'
        #                     '         pair_coeff command contains a \"*\" character.\n'
        #                     '         Keep in mind that using moltemplate.sh you can manually change the\n'
        #                     '         numbers assigned to each atom type (when using -a or -b).  Make sure\n'
        #                     '         nor to accidentally change the order of atom types in one of these\n'
        #                     '         pair_coeff commands.  For example, commands like\n'
        #                     '            pair_coeff 10*4 20*10 0.15 3.6\n'
        #                     '         can be generated by moltemplate.sh, however\n'
        #                     '         they may be rejected by LAMMPS (because LAMMPS prefers this\n'
        #                     '            pair_coeff 4*10 10*20 0.15 3.6)\n'
        #                     '         Later on, you may want to check to make sure moltemplate.sh\n'
        #                     '         is not doing this.  (Fortunately you never have to worry unless\n'
        #                     '         you are using the -a or -b arguments with moltemplate.sh)\n')

        if ('*' in atomtype_i_str):
            atomtype_i_tokens = atomtype_i_str.split('*')

            if atomtype_i_tokens[0] == '':
                if (min_sel_atomtype and 
                    (min_sel_atomtype < min_needed_atomtype)):
                    i_a = min_sel_atomtype
                else:
                    i_a = min_needed_atomtype
            else:
                i_a = Intify(atomtype_i_tokens[0])

            if atomtype_i_tokens[1] == '':
                if (max_sel_atomtype and 
                    (max_sel_atomtype > max_needed_atomtype)):
                    i_b = max_sel_atomtype
                else:
                    i_b = max_needed_atomtype
            else:
                i_b = Intify(atomtype_i_tokens[1])

        else:
            i_a = i_b = Intify(atomtype_i_str)

        i_a_final = None
        i_b_final = None
        for i in range(i_a, i_b+1):
            if ((i in needed_atomtypes) or (min_sel_atomtype <= i)):
                i_a_final = i
                break
        for i in reversed(range(i_a, i_b+1)):
            if ((i in needed_atomtypes) or (max_sel_atomtype >= i)):
                i_b_final = i
                break

        #if i_a_final and i_b_final:
        #    if i_a_final == i_b_final:
        #        i_str = '@atom:type'+str(i_a_final)
        #        tokens[1] = i_str
        #    else:
        #        i_str = '@{atom:type'+str(i_a_final)+'}*@{atom:type'+str(i_b_final)+'}'



        if ('*' in atomtype_j_str):
            atomtype_j_tokens = atomtype_j_str.split('*')

            if atomtype_j_tokens[0] == '':
                if (min_sel_atomtype and 
                    (min_sel_atomtype < min_needed_atomtype)):
                    j_a = min_sel_atomtype
                else:
                    j_a = min_needed_atomtype
            else:
                j_a = Intify(atomtype_j_tokens[0])

            if atomtype_j_tokens[1] == '':
                if (max_sel_atomtype and 
                    (max_sel_atomtype > max_needed_atomtype)):
                    j_b = max_sel_atomtype
                else:
                    j_b = max_needed_atomtype
            else:
                j_b = Intify(atomtype_j_tokens[1])

        else:
            j_a = j_b = Intify(atomtype_j_str)

        j_a_final = None
        j_b_final = None
        for j in range(j_a, j_b+1):
            if ((j in needed_atomtypes) or (min_sel_atomtype <= j)):
                j_a_final = j
                break
        for j in reversed(range(j_a, j_b+1)):
            if ((j in needed_atomtypes) or (max_sel_atomtype >= j)):
                j_b_final = j
                break

        #if j_a_final and j_b_final:
        #    if j_a_final == j_b_final:
        #        j_str = '@atom:type'+str(j_a_final)
        #        tokens[1] = j_str
        #    else:
        #        j_str = '@{atom:type'+str(j_a_final)+'}*@{atom:type'+str(j_b_final)+'}'



        if not (i_a_final and i_b_final and j_a_final and j_b_final):
            del l_in_pair_coeffs[i_line]
        elif (('*' in atomtype_i_str) or ('*' in atomtype_j_str)):
            del l_in_pair_coeffs[i_line]
            for i in range(i_a_final, i_b_final+1):
                for j in range(j_a_final, j_b_final+1):
                    if j >= i:
                        #tokens[1] = '@atom:type'+str(i)
                        #tokens[2] = '@atom:type'+str(j)
                        tokens[1] = '@atom:'+atomtypes_int2name[i]
                        tokens[2] = '@atom:'+atomtypes_int2name[j]
                        l_in_pair_coeffs.insert(i_line, 
                                                (' '*indent)+(' '.join(tokens)+'\n'))
                        i_line += 1
        else:
            #tokens[1] = '@atom:type'+tokens[1]
            #tokens[2] = '@atom:type'+tokens[2]
            tokens[1] = '@atom:'+atomtypes_int2name[int(tokens[1])]
            tokens[2] = '@atom:'+atomtypes_int2name[int(tokens[2])]
            l_in_pair_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1




    # delete mass commands for atom types we don't care about anymore:
    i_line = 0
    while i_line < len(l_in_masses):
        line = l_in_masses[i_line]
        tokens = line.strip().split()
        atomtype_i_str = tokens[1]
        #if (('*' in atomtype_i_str) or
        #    ('*' in atomtype_j_str)):
        #    sys.stderr.write('WARNING: near or before '+ErrorLeader(infile, lineno)+'\n'
        #                     '         pair_coeff command contains a \"*\" character.\n'
        #                     '         Keep in mind that using moltemplate.sh you can manually change the\n'
        #                     '         numbers assigned to each atom type (when using -a or -b).  Make sure\n'
        #                     '         nor to accidentally change the order of atom types in one of these\n'
        #                     '         pair_coeff commands.  For example, commands like\n'
        #                     '            pair_coeff 10*4 20*10 0.15 3.6\n'
        #                     '         can be generated by moltemplate.sh, however\n'
        #                     '         they may be rejected by LAMMPS (because LAMMPS prefers this\n'
        #                     '            pair_coeff 4*10 10*20 0.15 3.6)\n'
        #                     '         Later on, you may want to check to make sure moltemplate.sh\n'
        #                     '         is not doing this.  (Fortunately you never have to worry unless\n'
        #                     '         you are using the -a or -b arguments with moltemplate.sh)\n')

        if ('*' in atomtype_i_str):
            atomtype_i_tokens = atomtype_i_str.split('*')

            if atomtype_i_tokens[0] == '':
                if (min_sel_atomtype and 
                    (min_sel_atomtype < min_needed_atomtype)):
                    i_a = min_sel_atomtype
                else:
                    i_a = min_needed_atomtype
            else:
                i_a = Intify(atomtype_i_tokens[0])

            if atomtype_i_tokens[1] == '':
                if (max_sel_atomtype and 
                    (max_sel_atomtype > max_needed_atomtype)):
                    i_b = max_sel_atomtype
                else:
                    i_b = max_needed_atomtype
            else:
                i_b = Intify(atomtype_i_tokens[1])

        else:
            i_a = i_b = Intify(atomtype_i_str)

        i_a_final = None
        i_b_final = None
        for i in range(i_a, i_b+1):
            if ((i in needed_atomtypes) or (min_sel_atomtype <= i)):
                i_a_final = i
                break
        for i in reversed(range(i_a, i_b+1)):
            if ((i in needed_atomtypes) or (max_sel_atomtype >= i)):
                i_b_final = i
                break
        #if i_a_final and i_b_final:
        #    if i_a_final == i_b_final:
        #        i_str = '@atom:type'+str(i_a_final)
        #        tokens[1] = i_str
        #    else:
        #        i_str = '@{atom:type'+str(i_a_final)+'}*@{atom:type'+str(i_b_final)+'}'

        if not (i_a_final and i_b_final and j_a_final and j_b_final):
            del l_in_masses[i_line]
        elif ('*' in atomtype_i_str):
            del l_in_masses[i_line]
            for i in range(i_a_final, i_b_final+1):
                #tokens[1] = '@atom:type'+str(i)
                tokens[1] = '@atom:'+atomtypes_int2name[i]
                # CONTINUEHERE: CHECK THAT THIS IS WORKING
                l_in_masses.insert(i_line, (' '*indent)+(' '.join(tokens)+'\n'))
                i_line += 1
        else:
            assert(i_a == i_b)
            #tokens[1] = '@atom:type'+str(i_a)
            tokens[1] = '@atom:'+atomtypes_int2name[i_a]
            # CONTINUEHERE: CHECK THAT THIS IS WORKING
            l_in_masses[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1




    # --- BONDS AND BOND COEFFS ---

    # delete lines from data_bonds if they involve atoms we don't care about
    i_line = 0
    while i_line < len(l_data_bonds):
        line = l_data_bonds[i_line]
        tokens = line.strip().split()
        assert(len(tokens) == 4)

        bondid   = Intify(tokens[0])
        bondtype = Intify(tokens[1])
        atomid1  = Intify(tokens[2])
        atomid2  = Intify(tokens[3])
        #if ((atomid1 in needed_atomids) and
        #    (atomid2 in needed_atomids)):
        tokens[0] = '$bond:id'+str(bondid)
        tokens[1] = '@bond:type'+str(bondtype)
        #tokens[2] = '$atom:id'+str(atomid1)
        #tokens[3] = '$atom:id'+str(atomid2)
        tokens[2] = '$atom:'+atomids_int2name[atomid1]
        tokens[3] = '$atom:'+atomids_int2name[atomid2]
        needed_bondids.add(bondid)
        needed_bondtypes.add(bondtype)
        l_data_bonds[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
        i_line += 1
        #else:
        #    del l_data_bonds[i_line]

    # delete data_bond_coeffs for bondtypes we don't care about anymore:
    i_line = 0
    while i_line < len(l_data_bond_coeffs):
        line = l_data_bond_coeffs[i_line]
        tokens = line.strip().split()
        bondtype = Intify(tokens[0])
        if (not (bondtype in needed_bondtypes)):
            del l_data_bond_coeffs[i_line]
        else:
            tokens[0] = '@bond:type'+str(bondtype)
            l_data_bond_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1

    # delete in_bond_coeffs for bondtypes we don't care about anymore:
    for bondtype in needed_bondtypes:
        assert(type(bondtype) is int)
        if ((min_needed_bondtype == None) or
            (min_needed_bondtype > bondtype)):
            min_needed_bondtype = bondtype
        if ((max_needed_bondtype == None) or
            (max_needed_bondtype < bondtype)):
            max_needed_bondtype = bondtype
    for bondid in needed_bondids:
        assert(type(bondid) is int)
        if ((min_needed_bondid == None) or
            (min_needed_bondid > bondid)):
            min_needed_bondid = bondid
        if ((max_needed_bondid == None) or
            (max_needed_bondid < bondid)):
            max_needed_bondid = bondid


    i_line = 0
    while i_line < len(l_in_bond_coeffs):
        line = l_in_bond_coeffs[i_line]
        tokens = line.strip().split()
        bondtype_str = tokens[1]

        if ('*' in bondtype_str):
            bondtype_tokens = bondtype_str.split('*')

            if bondtype_tokens[0] == '':
                i_a = min_needed_bondtype
            else:
                i_a = Intify(bondtype_tokens[0])

            if bondtype_tokens[1] == '':
                i_b = max_needed_bondtype
            else:
                i_b = Intify(bondtype_tokens[1])

        else:
            i_a = Intify(bondtype_str)
            i_b = i_a

        if i_a < min_needed_bondtype:
            i_a = min_needed_bondtype
        if i_b > max_needed_bondtype:
            i_b = max_needed_bondtype

        #if i_a == i_b:
        #    i_str = '@bond:type'+str(i_a)
        #    tokens[1] = i_str
        #else:
        #    i_str = '@{bond:type'+str(j_a)+'}*@{bond:type'+str(j_b)+'}'

        if ('*' in bondtype_str):
            del l_in_bond_coeffs[i_line]
            for i in range(i_a, i_b+1):
                if (i in needed_bondtypes):
                    tokens[1] = '@bond:type'+str(i)
                    l_in_bond_coeffs.insert(i_line, 
                                            (' '*indent)+(' '.join(tokens)+'\n'))
                    i_line += 1
        else:
            if i_a < i_b:
                raise InputError('Error: number of bond types in data file is not consistent with the\n'
                                 '       number of bond types you have define bond_coeffs for.\n')
            if (i_a == i_b) and (i_a in needed_bondtypes):
                tokens[1] = '@bond:type'+str(i_a)
                l_in_bond_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
                i_line += 1
            else:
                del l_in_bond_coeffs[i_line]





    # --- ANGLES AND ANGLE COEFFS ---

    # delete lines from data_angles if they involve atoms we don't care about
    i_line = 0
    while i_line < len(l_data_angles):
        line = l_data_angles[i_line]
        tokens = line.strip().split()
        assert(len(tokens) == 5)
        
        angleid   = Intify(tokens[0])
        angletype = Intify(tokens[1])
        atomid1  = Intify(tokens[2])
        atomid2  = Intify(tokens[3])
        atomid3  = Intify(tokens[4])
        #if ((atomid1 in needed_atomids) and
        #    (atomid2 in needed_atomids)):
        tokens[0] = '$angle:id'+str(angleid)
        tokens[1] = '@angle:type'+str(angletype)
        #tokens[2] = '$atom:id'+str(atomid1)
        #tokens[3] = '$atom:id'+str(atomid2)
        #tokens[4] = '$atom:id'+str(atomid3)
        tokens[2] = '$atom:'+atomids_int2name[atomid1]
        tokens[3] = '$atom:'+atomids_int2name[atomid2]
        tokens[4] = '$atom:'+atomids_int2name[atomid3]
        needed_angleids.add(angleid)
        needed_angletypes.add(angletype)
        l_data_angles[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
        i_line += 1
        #else:
        #    del l_data_angles[i_line]

    # delete data_angle_coeffs for angletypes we don't care about anymore:
    i_line = 0
    while i_line < len(l_data_angle_coeffs):
        line = l_data_angle_coeffs[i_line]
        tokens = line.strip().split()
        angletype = Intify(tokens[0])
        if (not (angletype in needed_angletypes)):
            del l_data_angle_coeffs[i_line]
        else:
            tokens[0] = '@angle:type'+str(angletype)
            l_data_angle_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1

    # delete in_angle_coeffs for angletypes we don't care about anymore:
    for angletype in needed_angletypes:
        assert(type(angletype) is int)
        if ((min_needed_angletype == None) or
            (min_needed_angletype > angletype)):
            min_needed_angletype = angletype
        if ((max_needed_angletype == None) or
            (max_needed_angletype < angletype)):
            max_needed_angletype = angletype
    for angleid in needed_angleids:
        assert(type(angleid) is int)
        if ((min_needed_angleid == None) or
            (min_needed_angleid > angleid)):
            min_needed_angleid = angleid
        if ((max_needed_angleid == None) or
            (max_needed_angleid < angleid)):
            max_needed_angleid = angleid

    i_line = 0
    while i_line < len(l_in_angle_coeffs):
        line = l_in_angle_coeffs[i_line]
        tokens = line.strip().split()
        angletype_str = tokens[1]

        if ('*' in angletype_str):
            angletype_tokens = angletype_str.split('*')

            if angletype_tokens[0] == '':
                i_a = min_needed_angletype
            else:
                i_a = Intify(angletype_tokens[0])

            if angletype_tokens[1] == '':
                i_b = max_needed_angletype
            else:
                i_b = Intify(angletype_tokens[1])

        else:
            i_a = i_b = Intify(angletype_str)

        if i_a < min_needed_angletype:
            i_a = min_needed_angletype
        if i_b > max_needed_angletype:
            i_b = max_needed_angletype

        #if i_a == i_b:
        #    i_str = '@angle:type'+str(i_a)
        #    tokens[1] = i_str
        #else:
        #    i_str = '@{angle:type'+str(j_a)+'}*@{angle:type'+str(j_b)+'}'

        if ('*' in angletype_str):
            del l_in_angle_coeffs[i_line]
            for i in range(i_a, i_b+1):
                if (i in needed_angletypes):
                    tokens[1] = '@angle:type'+str(i)
                    l_in_angle_coeffs.insert(i_line, 
                                             (' '*indent)+(' '.join(tokens)+'\n'))
                    i_line += 1
        else:
            if i_a < i_b:
                raise InputError('Error: number of angle types in data file is not consistent with the\n'
                                 '       number of angle types you have define angle_coeffs for.\n')
            if (i_a == i_b) and (i_a in needed_angletypes):
                tokens[1] = '@angle:type'+str(i_a)
                l_in_angle_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
                i_line += 1
            else:
                del l_in_angle_coeffs[i_line]



    # --- DIHEDRALS AND DIHEDRAL COEFFS ---

    # delete lines from data_dihedrals if they involve atoms we don't care about
    i_line = 0
    while i_line < len(l_data_dihedrals):
        line = l_data_dihedrals[i_line]
        tokens = line.strip().split()
        assert(len(tokens) == 6)
        
        dihedralid   = Intify(tokens[0])
        dihedraltype = Intify(tokens[1])
        atomid1  = Intify(tokens[2])
        atomid2  = Intify(tokens[3])
        atomid3  = Intify(tokens[4])
        atomid4  = Intify(tokens[5])
        #if ((atomid1 in needed_atomids) and
        #    (atomid2 in needed_atomids)):
        tokens[0] = '$dihedral:id'+str(dihedralid)
        tokens[1] = '@dihedral:type'+str(dihedraltype)
        #tokens[2] = '$atom:id'+str(atomid1)
        #tokens[3] = '$atom:id'+str(atomid2)
        #tokens[4] = '$atom:id'+str(atomid3)
        #tokens[5] = '$atom:id'+str(atomid4)
        tokens[2] = '$atom:'+atomids_int2name[atomid1]
        tokens[3] = '$atom:'+atomids_int2name[atomid2]
        tokens[4] = '$atom:'+atomids_int2name[atomid3]
        tokens[5] = '$atom:'+atomids_int2name[atomid4]

        needed_dihedralids.add(dihedralid)
        needed_dihedraltypes.add(dihedraltype)
        l_data_dihedrals[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
        i_line += 1
        #else:
        #    del l_data_dihedrals[i_line]

    # delete data_dihedral_coeffs for dihedraltypes we don't care about anymore:
    i_line = 0
    while i_line < len(l_data_dihedral_coeffs):
        line = l_data_dihedral_coeffs[i_line]
        tokens = line.strip().split()
        dihedraltype = Intify(tokens[0])
        if (not (dihedraltype in needed_dihedraltypes)):
            del l_data_dihedral_coeffs[i_line]
        else:
            tokens[0] = '@dihedral:type'+str(dihedraltype)
            l_data_dihedral_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1

    # delete in_dihedral_coeffs for dihedraltypes we don't care about anymore:
    for dihedraltype in needed_dihedraltypes:
        assert(type(dihedraltype) is int)
        if ((min_needed_dihedraltype == None) or
            (min_needed_dihedraltype > dihedraltype)):
            min_needed_dihedraltype = dihedraltype
        if ((max_needed_dihedraltype == None) or
            (max_needed_dihedraltype < dihedraltype)):
            max_needed_dihedraltype = dihedraltype
    for dihedralid in needed_dihedralids:
        assert(type(dihedralid) is int)
        if ((min_needed_dihedralid == None) or
            (min_needed_dihedralid > dihedralid)):
            min_needed_dihedralid = dihedralid
        if ((max_needed_dihedralid == None) or
            (max_needed_dihedralid < dihedralid)):
            max_needed_dihedralid = dihedralid

    i_line = 0
    while i_line < len(l_in_dihedral_coeffs):
        line = l_in_dihedral_coeffs[i_line]
        tokens = line.strip().split()
        dihedraltype_str = tokens[1]

        if ('*' in dihedraltype_str):
            dihedraltype_tokens = dihedraltype_str.split('*')

            if dihedraltype_tokens[0] == '':
                i_a = min_needed_dihedraltype
            else:
                i_a = Intify(dihedraltype_tokens[0])

            if dihedraltype_tokens[1] == '':
                i_b = max_needed_dihedraltype
            else:
                i_b = Intify(dihedraltype_tokens[1])

        else:
            i_a = i_b = Intify(dihedraltype_str)

        if i_a < min_needed_dihedraltype:
            i_a = min_needed_dihedraltype
        if i_b > max_needed_dihedraltype:
            i_b = max_needed_dihedraltype

        #if i_a == i_b:
        #    i_str = '@dihedral:type'+str(i_a)
        #    tokens[1] = i_str
        #else:
        #    i_str = '@{dihedral:type'+str(j_a)+'}*@{dihedral:type'+str(j_b)+'}'

        if ('*' in dihedraltype_str):
            del l_in_dihedral_coeffs[i_line]
            for i in range(i_a, i_b+1):
                if (i in needed_dihedraltypes):
                    tokens[1] = '@dihedral:type'+str(i)
                    l_in_dihedral_coeffs.insert(i_line, 
                                                (' '*indent)+(' '.join(tokens)+'\n'))
                    i_line += 1
        else:
            if i_a < i_b:
                raise InputError('Error: number of dihedral types in data file is not consistent with the\n'
                                 '       number of dihedral types you have define dihedral_coeffs for.\n')
            if (i_a == i_b) and (i_a in needed_dihedraltypes):
                tokens[1] = '@dihedral:type'+str(i_a)
                l_in_dihedral_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
                i_line += 1
            else:
                del l_in_dihedral_coeffs[i_line]



    # --- IMPROPERS AND IMPROPER COEFFS ---

    # delete lines from data_impropers if they involve atoms we don't care about
    i_line = 0
    while i_line < len(l_data_impropers):
        line = l_data_impropers[i_line]
        tokens = line.strip().split()
        assert(len(tokens) == 6)
        
        improperid   = Intify(tokens[0])
        impropertype = Intify(tokens[1])
        atomid1  = Intify(tokens[2])
        atomid2  = Intify(tokens[3])
        atomid3  = Intify(tokens[4])
        atomid4  = Intify(tokens[5])
        #if ((atomid1 in needed_atomids) and
        #    (atomid2 in needed_atomids)):
        tokens[0] = '$improper:id'+str(improperid)
        tokens[1] = '@improper:type'+str(impropertype)
        #tokens[2] = '$atom:id'+str(atomid1)
        #tokens[3] = '$atom:id'+str(atomid2)
        #tokens[4] = '$atom:id'+str(atomid3)
        #tokens[5] = '$atom:id'+str(atomid4)
        tokens[2] = '$atom:'+atomids_int2name[atomid1]
        tokens[3] = '$atom:'+atomids_int2name[atomid2]
        tokens[4] = '$atom:'+atomids_int2name[atomid3]
        tokens[5] = '$atom:'+atomids_int2name[atomid4]

        needed_improperids.add(improperid)
        needed_impropertypes.add(impropertype)
        l_data_impropers[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
        i_line += 1
        #else:
        #    del l_data_impropers[i_line]

    # delete data_improper_coeffs for impropertypes we don't care about anymore:
    i_line = 0
    while i_line < len(l_data_improper_coeffs):
        line = l_data_improper_coeffs[i_line]
        tokens = line.strip().split()
        impropertype = Intify(tokens[0])
        if (not (impropertype in needed_impropertypes)):
            del l_data_improper_coeffs[i_line]
        else:
            tokens[0] = '@improper:type'+str(impropertype)
            l_data_improper_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
            i_line += 1

    # delete in_improper_coeffs for impropertypes we don't care about anymore:
    for impropertype in needed_impropertypes:
        assert(type(impropertype) is int)
        if ((min_needed_impropertype == None) or
            (min_needed_impropertype > impropertype)):
            min_needed_impropertype = impropertype
        if ((max_needed_impropertype == None) or
            (max_needed_impropertype < impropertype)):
            max_needed_impropertype = impropertype
    for improperid in needed_improperids:
        assert(type(improperid) is int)
        if ((min_needed_improperid == None) or
            (min_needed_improperid > improperid)):
            min_needed_improperid = improperid
        if ((max_needed_improperid == None) or
            (max_needed_improperid < improperid)):
            max_needed_improperid = improperid

    i_line = 0
    while i_line < len(l_in_improper_coeffs):
        line = l_in_improper_coeffs[i_line]
        tokens = line.strip().split()
        impropertype_str = tokens[1]

        if ('*' in impropertype_str):
            impropertype_tokens = impropertype_str.split('*')

            if impropertype_tokens[0] == '':
                i_a = min_needed_impropertype
            else:
                i_a = Intify(impropertype_tokens[0])

            if impropertype_tokens[1] == '':
                i_b = max_needed_impropertype
            else:
                i_b = Intify(impropertype_tokens[1])

        else:
            i_a = i_b = Intify(impropertype_str)

        if i_a < min_needed_impropertype:
            i_a = min_needed_impropertype
        if i_b > max_needed_impropertype:
            i_b = max_needed_impropertype

        #if i_a == i_b:
        #    i_str = '@improper:type'+str(i_a)
        #    tokens[1] = i_str
        #else:
        #    i_str = '@{improper:type'+str(j_a)+'}*@{improper:type'+str(j_b)+'}'

        if ('*' in impropertype_str):
            del l_in_improper_coeffs[i_line]
            for i in range(i_a, i_b+1):
                if (i in needed_impropertypes):
                    tokens[1] = '@improper:type'+str(i)
                    l_in_improper_coeffs.insert(i_line, 
                                                (' '*indent)+(' '.join(tokens)+'\n'))
                    i_line += 1
        else:
            if i_a < i_b:
                raise InputError('Error: number of improper types in data file is not consistent with the\n'
                                 '       number of improper types you have define improper_coeffs for.\n')
            if (i_a == i_b) and (i_a in needed_impropertypes):
                tokens[1] = '@improper:type'+str(i_a)
                l_in_improper_coeffs[i_line] = (' '*indent)+(' '.join(tokens)+'\n')
                i_line += 1
            else:
                del l_in_improper_coeffs[i_line]

    # --- GROUPS ---

    # Now parse through all of the "group" commands and try and figure
    # out if any of these groups contain any of the atoms we are keeping.
    # If so, then save the group and write it out.
    # (I hate trying to parse this kind of text.)

    if len(l_in_group) > 0:
        sys.stderr.write('\n'
                         ' --groups--  Attempting to parse \"group\" commands.\n'
                         '         This may cause '+g_program_name+' to crash.\n'
                         '         If so, comment out all group commands in your input script(s), and\n'
                         '         try again.  (And please report the error. -Andrew 2014-10-30)\n')

    i_line = 0
    groups_needed = set(['all'])
    while i_line < len(l_in_group):
        line = l_in_group[i_line]
        tokens = line.strip().split()
        delete_this_command = False
        explicit_definition = False
        if len(tokens) < 3:
            delete_this_command = True
        group_name = tokens[1]
        specifier_style = tokens[2]
        str_logical = ''
        str_selection = ''
        if specifier_style[0:4] == 'type':
            str_logical+=specifier_style[4:]
            explicit_definition = True
            specifier_style = 'type'
        elif specifier_style == 'id':
            str_logical+=specifier_style[2:]
            explicit_definition = True
            specifier_style = 'id'
        elif specifier_style == 'molecule':
            str_logical+=specifier_style[8:]
            specifier_style = 'molecule'
            explicit_definition = True

        if explicit_definition:
            i_token_sel_min = 3
            if len(tokens) <= i_token_sel_min:
                sys.stderr.write('WARNING: possible syntax error on this line:\n'
                                 +'        '+l_in_group[i_line]+'\n')
                delete_this_command = True
            if str_logical == '':
                str_logical = tokens[i_token_sel_min]
                if not str_logical[0].isdigit():
                    i_token_sel_min += 1
                    if len(tokens) <= i_token_sel_min:
                        tokens.append('')
            else:
                tokens.insert(i_token_sel_min, str_logical)

            i_token_sel_max = len(tokens)-1

            for i in range(i_token_sel_min, len(tokens)):
                if tokens[i].isdigit():
                    break
                else:
                    i_token_sel_max = i

            assert(len(tokens) > i_token_sel_min)

            if str_logical[0:2] in ('<=','>=','==','!=','<>'):
                tokens[i_token_sel_min] = str_logical[2:] + tokens[i_token_sel_min]
                str_logical = str_logical[0:2]
                if str_logical == '<=':
                    l_group_selection = [ (None,int(tokens[i_token_sel_min])) ]
                elif str_logical == '>=':
                    l_group_selection = [ (int(tokens[i_token_sel_min]),None) ]
                elif str_logical == '==':
                    l_group_selection = [ (int(tokens[i_token_sel_min]),
                                           int(tokens[i_token_sel_min])) ]
                elif str_logical == '!=':
                    l_group_selection = [ (None,int(tokens[i_token_sel_min])-1),
                                          (int(tokens[i_token_sel_min])+1,None)]
                elif str_logical == '<>':
                    l_group_selection = [ (int(tokens[i_token_sel_min]),
                                           int(tokens[i_token_sel_max])) ]

            elif str_logical[0:1] in ('<','>'):
                tokens[i_token_sel_min] = str_logical[1:] + tokens[i_token_sel_min]
                str_logical = str_logical[0:1]
                if str_logical == '<':
                    l_group_selection = [(None,int(tokens[i_token_sel_min])-1)]
                elif str_logical == '>':
                    l_group_selection = [(int(tokens[i_token_sel_min])+1,None)]
            else:
                str_selection = ' '.join(tokens[i_token_sel_min:i_token_sel_max+1])
                l_group_selection = LammpsSelectToIntervals(str_selection,
                                                            slice_delim=':',
                                                            or_delim=' ')

            mn, mx = IntervalListToMinMax(l_group_selection)
            if mn == None:
                mn = 1
            filtered_selection=[]
            if specifier_style == 'type':
                if mx == None:
                    mx = max_needed_atomtype
                for i in range(mn, mx+1):
                    if (BelongsToSel(i, l_group_selection)
                        and (i in needed_atomtypes)):
                        filtered_selection.append((i,i))
            elif specifier_style == 'id':
                if mx == None:
                    mx = max_needed_atomid
                for i in range(mn, mx+1):
                    if (BelongsToSel(i, l_group_selection)
                        and (i in needed_atomids)):
                        filtered_selection.append((i,i))
            elif specifier_style == 'molecule':
                if mx == None:
                    mx = max_needed_molid
                for i in range(mn, mx+1):
                    if (BelongsToSel(i, l_group_selection)
                        and (i in needed_molids)):
                        filtered_selection.append((i,i))


            MergeIntervals(filtered_selection)


            if len(filtered_selection) > 0:
                
                tokens = ['group', group_name, specifier_style]
                for interval in filtered_selection:
                    a = interval[0]
                    b = interval[1]

                    if specifier_style == 'type':
                        if a == b:
                            tokens.append('@atom:type'+str(a))
                        else:
                            tokens.append('@{atom:type'+str(a)+
                                                    '}:@{atom:type'+str(b)+'}')

                    if specifier_style == 'id':
                        if a == b:
                            tokens.append('$atom:id'+str(a))
                        else:
                            tokens.append('${atom:id'+str(a)
                                                    +'}:${atom:id'+str(b)+'}')

                    if specifier_style == 'molecule':
                        if a == b:
                            tokens.append('$mol:id'+str(a))
                        else:
                            tokens.append('${mol:id'+str(a)+
                                                    '}:${mol:id'+str(b)+'}')


                # Commenting out next two lines.  (This is handled later.)
                #l_in_group[i_line] = ' '.join(tokens)
                #groups_needed.add(group_name)

            else:
                delete_this_command = True


        else:
            if len(tokens) > 3:
                if tokens[2] == 'union':
                    i_token = 3
                    while i_token < len(tokens):
                        if not (tokens[i_token] in groups_needed):
                            del tokens[i_token]
                        else:
                            i_token += 1
                    # if none of the groups contain atoms we need,
                    # then delete the entire command
                    if len(tokens) <= 3:
                        delete_this_command = True
                elif tokens[2] == 'intersect':
                    i_token = 3
                    while i_token < len(tokens):
                        if not (tokens[i_token] in groups_needed):
                            # if any of the groups we need are empty
                            # then delete the command
                            delete_this_command = True
                            break
                        i_token += 1
                elif (tokens[2] == 'subtract') and (len(tokens) >= 5):
                    if not (tokens[3] in groups_needed):
                        delete_this_command = True
                    i_token = 4
                    while i_token < len(tokens):
                        if not (tokens[i_token] in groups_needed):
                            del tokens[i_token]
                        else:
                            i_token += 1
                else:
                    # Otherwise I don't recongize the syntax of this
                    # group command.  In that case, I just delete it.
                    delete_this_command = True

            elif tokens[2] == 'clear':
                pass
            elif tokens[2] == 'delete':
                pass
            else:
                delete_this_command = True
        if delete_this_command:
            sys.stderr.write('WARNING: Ignoring line \n\"'+l_in_group[i_line].rstrip()+'\"\n')
            del l_in_group[i_line]
        else:
            groups_needed.add(group_name)
            l_in_group[i_line] = (' '*indent) + ' '.join(tokens) + '\n'
            i_line += 1



    # --- fix rigid ---

    i_line = 0
    while i_line < len(l_in_fix_rigid):
        line = l_in_fix_rigid[i_line]
        tokens = line.strip().split()
        if len(tokens) < 4:
            break
        fixid = tokens[1]
        group_name = tokens[2]
        delete_this_command = True
        assert(tokens[3].find('rigid') == 0)
        if group_name in groups_needed:
            delete_this_command = False

        if delete_this_command:
            sys.stderr.write('WARNING: Ignoring line \n\"'+l_in_fix_rigid[i_line].rstrip()+'\"\n')
            del l_in_fix_rigid[i_line]
        else:
            l_in_fix_rigid[i_line] = (' '*indent) + ' '.join(tokens) + '\n'
            i_line += 1



    # --- fix shake ---

    i_line = 0
    while i_line < len(l_in_fix_shake):
        line = l_in_fix_shake[i_line]
        tokens = line.strip().split()
        if len(tokens) < 4:
            break
        fixid = tokens[1]
        group_name = tokens[2]
        delete_this_command = True
        assert(tokens[3].find('shake') == 0)

        #  parse the list of angle types
        #i_token = tokens.index('a')
        for i_token in range(0, len(tokens)):
            if tokens[i_token] == 'a':
                break
        if i_token != len(tokens):
            i_token += 1
            while (i_token < len(tokens)) and tokens[i_token].isdigit():
                # delete angle types from the list which 
                # do not belong to the selection
                btype=int(tokens[i_token])
                if int(tokens[i_token]) in needed_angletypes:
                    tokens[i_token] = '@angle:type'+tokens[i_token]
                    i_token += 1
                    delete_this_command = False
                else:
                    del tokens[i_token]

        #  parse the list of bond types
        #i_token = tokens.index('b')
        for i_token in range(0, len(tokens)):
            if tokens[i_token] == 'b':
                break
        if i_token != len(tokens):
            i_token += 1
            while (i_token < len(tokens)) and tokens[i_token].isdigit():
                # delete bond types from the list which 
                # do not belong to the selection
                btype=int(tokens[i_token])
                if int(tokens[i_token]) in needed_bondtypes:
                    tokens[i_token] = '@bond:type'+tokens[i_token]
                    i_token += 1
                    delete_this_command = False
                else:
                    del tokens[i_token]

        #  parse the list of atom types
        # i_token = tokens.index('t')
        for i_token in range(0, len(tokens)):
            if tokens[i_token] == 't':
                break
        if i_token != len(tokens):
            i_token += 1
            while (i_token < len(tokens)) and tokens[i_token].isdigit():
                # delete atom types from the list which 
                # do not belong to the selection
                btype=int(tokens[i_token])
                if int(tokens[i_token]) in needed_atomtypes:
                    tokens[i_token] = '@atom:type'+tokens[i_token]
                    i_token += 1
                    delete_this_command = False
                else:
                    del tokens[i_token]


        #  Selecting atoms by mass feature should still work, so we 
        #  don't need to delete or ignore these kinds of commands.
        #for i_token in range(0, len(tokens)):
        #    if tokens[i_token] == 'm':
        #        break
        #if i_token != len(tokens):
        #    delete_this_command = True

        if 'mol' in tokens:
            delete_this_command = True

        if not (group_name in groups_needed):
            delete_this_command = True

        if delete_this_command:
            sys.stderr.write('WARNING: Ignoring line \n\"'+l_in_fix_shake[i_line].rstrip()+'\"\n')
            del l_in_fix_shake[i_line]
        else:
            l_in_fix_shake[i_line] = (' '*indent) + ' '.join(tokens) + '\n'
            i_line += 1


    # --- fix poems ---

    i_line = 0
    while i_line < len(l_in_fix_poems):
        line = l_in_fix_poems[i_line]
        tokens = line.strip().split()
        if len(tokens) < 4:
            break
        fixid = tokens[1]
        group_name = tokens[2]
        delete_this_command = True
        assert(tokens[3].find('poems') == 0)
        if group_name in groups_needed:
            delete_this_command = False
        if tokens[4] != 'molecule':
            delete_this_command = True
            sys.stderr.write('WARNING: '+g_program_name+' ONLY supports \"fix poems\" commands\n'
                             '         which use the \"molecule\" keyword.\n')
        if tokens[4] == 'file':
            sys.stderr.write('         If you want use external files with fix poems, then you will have to\n'
                             '         generate the file yourself.  You ask use moltemplate to generate\n'
                             '         this file for you, by manually adding a section at the end of your\n'
                             '         final .LT file (eg. \"system.lt\") which resembles the following:\n\n'
                             'write(\"poems_file.txt\") {\n'
                             '  1 1 $atom:idname1a $atom:idname2a $atom:idname3a ...\n'
                             '  2 1 $atom:idname1b $atom:idname2b $atom:idname3b ...\n'
                             '  3 1 $atom:idname1c $atom:idname2c $atom:idname3c ...\n'
                             '  : :   etc...\n'
                             '}\n\n'
                             '      ...where $atom:idname1a, $atom:idname2a, ... are moltemplate-compatible\n'
                             '         unique (full,long) id-names for the atoms in each rigid body.\n'
                             '         This will insure the atom-id numbers in this file are correct.\n'

                             '         See the documentation for fix poems for details.\n')
                             

        if delete_this_command:
            sys.stderr.write('WARNING: Ignoring line \n\"'+l_in_fix_poems[i_line].rstrip()+'\"\n')
            del l_in_fix_poems[i_line]
        else:
            l_in_fix_poems[i_line] = (' '*indent) + ' '.join(tokens) + '\n'
            i_line += 1



    # --- fix qeq ---

    i_line = 0
    while i_line < len(l_in_fix_qeq):
        line = l_in_fix_qeq[i_line]
        tokens = line.strip().split()
        if len(tokens) < 4:
            break
        fixid = tokens[1]
        group_name = tokens[2]
        delete_this_command = True
        assert(tokens[3].find('qeq') == 0)
        if group_name in groups_needed:
            delete_this_command = False

        if delete_this_command:
            sys.stderr.write('WARNING: Ignoring line \n\"'+l_in_fix_qeq[i_line].rstrip()+'\"\n')
            del l_in_fix_qeq[i_line]
        else:
            l_in_fix_qeq[i_line] = (' '*indent) + ' '.join(tokens) + '\n'
            i_line += 1



    # --- fix qmmm ---

    i_line = 0
    while i_line < len(l_in_fix_qmmm):
        line = l_in_fix_qmmm[i_line]
        tokens = line.strip().split()
        if len(tokens) < 4:
            break
        fixid = tokens[1]
        group_name = tokens[2]
        delete_this_command = True
        assert(tokens[3].find('qmmm') == 0)
        if group_name in groups_needed:
            delete_this_command = False

        if delete_this_command:
            sys.stderr.write('WARNING: Ignoring line \n\"'+l_in_fix_qmmm[i_line].rstrip()+'\"\n')
            del l_in_fix_qmmm[i_line]
        else:
            l_in_fix_qmmm[i_line] = (' '*indent) + ' '.join(tokens) + '\n'
            i_line += 1






    ########################################
    ###  Now begin writing the template. ###
    ########################################

    if not some_pair_coeffs_read:
        sys.stderr.write('Warning: No \"pair coeffs\" set.\n'
                         '         (No interactions between non-bonded atoms defined.)\n')
        no_warnings = False

    #sys.stderr.write('Writing ttree data to standard out.\n'
    #                 '       You can redirect this to a file using:\n'+
    #                 '   '+' '.join(sys.argv)+' > filename.ttree\n'
    #                 '        ----------------------\n')

    if mol_name != '':
        sys.stdout.write(mol_name + ' {\n')

    if len(l_in_init) > 0:
        sys.stdout.write('\n### LAMMPS commands for initialization\n'
                         '### (These can be overridden later.)\n\n')
        l_in_init.insert(0, (' '*cindent)+'write_once(\"'+in_init+'\") {\n')
        l_in_init.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_init))
    if len(l_in_settings) > 0:
        sys.stdout.write('\n### LAMMPS commands for settings\n'
                         '### (These can be overridden later.)\n\n')
        l_in_settings.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_settings.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_settings))
        non_empty_output = True
    if len(l_in_masses) > 0:
        l_in_masses.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_masses.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_masses))
        non_empty_output = True
    if len(l_in_pair_coeffs) > 0:
        l_in_pair_coeffs.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_pair_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_pair_coeffs))
        non_empty_output = True
    if len(l_in_bond_coeffs) > 0:
        l_in_bond_coeffs.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_bond_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_bond_coeffs))
        non_empty_output = True
    if len(l_in_angle_coeffs) > 0:
        l_in_angle_coeffs.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_angle_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_angle_coeffs))
        non_empty_output = True
    if len(l_in_dihedral_coeffs) > 0:
        l_in_dihedral_coeffs.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_dihedral_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_dihedral_coeffs))
        non_empty_output = True
    if len(l_in_improper_coeffs) > 0:
        l_in_improper_coeffs.insert(0, (' '*cindent)+'write_once(\"'+in_settings+'\") {\n')
        l_in_improper_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_improper_coeffs))
        non_empty_output = True

    if non_empty_output:
        sys.stdout.write('\n\n### DATA sections\n\n')

    if len(l_data_masses) > 0:
        l_data_masses.insert(0, (' '*cindent)+'write_once(\"'+data_masses+'\") {\n')
        l_data_masses.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_masses))
        non_empty_output = True
    if len(l_data_bond_coeffs) > 0:
        l_data_bond_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_bond_coeffs+'\") {\n')
        l_data_bond_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_bond_coeffs))
        non_empty_output = True
    if len(l_data_angle_coeffs) > 0:
        l_data_angle_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_angle_coeffs+'\") {\n')
        l_data_angle_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_angle_coeffs))
        non_empty_output = True
    if len(l_data_dihedral_coeffs) > 0:
        l_data_dihedral_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_dihedral_coeffs+'\") {\n')
        l_data_dihedral_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_dihedral_coeffs))
        non_empty_output = True
    if len(l_data_improper_coeffs) > 0:
        l_data_improper_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_improper_coeffs+'\") {\n')
        l_data_improper_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_improper_coeffs))
        non_empty_output = True
    if len(l_data_pair_coeffs) > 0:
        l_data_pair_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_pair_coeffs+'\") {\n')
        l_data_pair_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_pair_coeffs))
        non_empty_output = True

    # class2 force fields:
    if len(l_data_bondbond_coeffs) > 0:
        l_data_bondbond_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_bondbond_coeffs+'\") {\n')
        l_data_bondbond_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_bondbond_coeffs))
        non_empty_output = True
    if len(l_data_bondangle_coeffs) > 0:
        l_data_bondangle_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_bondangle_coeffs+'\") {\n')
        l_data_bondangle_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_bondangle_coeffs))
        non_empty_output = True
    if len(l_data_middlebondtorsion_coeffs) > 0:
        l_data_middlebondtorsion_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_middlebondtorsion_coeffs+'\") {\n')
        l_data_middlebondtorsion_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_middlebondtorsion_coeffs))
        non_empty_output = True
    if len(l_data_endbondtorsion_coeffs) > 0:
        l_data_endbondtorsion_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_endbondtorsion_coeffs+'\") {\n')
        l_data_endbondtorsion_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_endbondtorsion_coeffs))
        non_empty_output = True
    if len(l_data_angletorsion_coeffs) > 0:
        l_data_angletorsion_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_angletorsion_coeffs+'\") {\n')
        l_data_angletorsion_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_angletorsion_coeffs))
        non_empty_output = True
    if len(l_data_angleangletorsion_coeffs) > 0:
        l_data_angleangletorsion_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_angleangletorsion_coeffs+'\") {\n')
        l_data_angleangletorsion_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_angleangletorsion_coeffs))
        non_empty_output = True
    if len(l_data_bondbond13_coeffs) > 0:
        l_data_bondbond13_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_bondbond13_coeffs+'\") {\n')
        l_data_bondbond13_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_bondbond13_coeffs))
        non_empty_output = True
    if len(l_data_angleangle_coeffs) > 0:
        l_data_angleangle_coeffs.insert(0, (' '*cindent)+'write_once(\"'+data_angleangle_coeffs+'\") {\n')
        l_data_angleangle_coeffs.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_angleangle_coeffs))
        non_empty_output = True

    # automatic generation of bonded interactions by type:
    if len(l_data_angles_by_type) > 0:
        l_data_angles_by_type.insert(0, (' '*cindent)+'write_once(\"'+data_angles_by_type+'\") {\n')
        l_data_angles_by_type.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_angles_by_type))
        non_empty_output = True
    if len(l_data_dihedrals_by_type) > 0:
        l_data_dihedrals_by_type.insert(0, (' '*cindent)+'write_once(\"'+data_dihedrals_by_type+'\") {\n')
        l_data_dihedrals_by_type.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_dihedrals_by_type))
        non_empty_output = True
    if len(l_data_impropers_by_type) > 0:
        l_data_impropers_by_type.insert(0, (' '*cindent)+'write_once(\"'+data_impropers_by_type+'\") {\n')
        l_data_impropers_by_type.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_impropers_by_type))
        non_empty_output = True

    if len(l_data_atoms) > 0:
        l_data_atoms.insert(0, (' '*cindent)+'write(\"'+data_atoms+'\") {\n')
        l_data_atoms.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_atoms))
        non_empty_output = True
    else:
        sys.stderr.write('Warning: missing \"Atoms\" section.\n'
                         '         (Did you include a LAMMPS data file in your argument list?)\n')
        no_warnings = False

    # non-point-like particles
    if len(l_data_ellipsoids) > 0:
        l_data_ellipsoids.insert(0, (' '*cindent)+'write(\"'+data_ellipsoids+'\") {\n')
        l_data_ellipsoids.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_ellipsoids))
    if len(l_data_lines) > 0:
        l_data_lines.insert(0, (' '*cindent)+'write(\"'+data_lines+'\") {\n')
        l_data_lines.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_lines))
    if len(l_data_triangles) > 0:
        l_data_triangles.insert(0, (' '*cindent)+'write(\"'+data_triangles+'\") {\n')
        l_data_triangles.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_triangles))

    if len(l_data_velocities) > 0:
        l_data_velocities.insert(0, (' '*cindent)+'write(\"'+data_velocities+'\") {\n')
        l_data_velocities.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_velocities))
    if len(l_data_bonds) > 0:
        l_data_bonds.insert(0, (' '*cindent)+'write(\"'+data_bonds+'\") {\n')
        l_data_bonds.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_bonds))
        non_empty_output = True
    if len(l_data_angles) > 0:
        l_data_angles.insert(0, (' '*cindent)+'write(\"'+data_angles+'\") {\n')
        l_data_angles.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_angles))
        non_empty_output = True
    if len(l_data_dihedrals) > 0:
        l_data_dihedrals.insert(0, (' '*cindent)+'write(\"'+data_dihedrals+'\") {\n')
        l_data_dihedrals.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_dihedrals))
        non_empty_output = True
    if len(l_data_impropers) > 0:
        l_data_impropers.insert(0, (' '*cindent)+'write(\"'+data_impropers+'\") {\n')
        l_data_impropers.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_data_impropers))
        non_empty_output = True 

    if len(l_in_group) > 0:
        no_warnings = False
        l_in_group.insert(0, (' '*cindent)+'write(\"'+in_settings+'\") {\n')
        l_in_group.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_group))
        sys.stderr.write('######################################################\n'
                         'WARNING: One or more \"group\" commands appear to refer to relevant atoms.\n'
                         '         Please check to make sure that the group(s) generated by\n'
                         '         '+g_program_name+' contain the correct atoms.  (-Andrew 2014-10-30)\n'
                         '######################################################\n')
        assert(non_empty_output)

    if len(l_in_fix_rigid) > 0:
        no_warnings = False
        l_in_fix_rigid.insert(0, (' '*cindent)+'write(\"'+in_settings+'\") {\n')
        l_in_fix_rigid.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_fix_rigid))
        sys.stderr.write('WARNING: \"fix rigid\" style command(s) applied to selected atoms.\n'
                         '         Please make sure that the fix group(s) are defined correctly.\n'
                         '######################################################\n')
        assert(non_empty_output)

    if len(l_in_fix_shake) > 0:
        no_warnings = False
        l_in_fix_shake.insert(0, (' '*cindent)+'write(\"'+in_settings+'\") {\n')
        l_in_fix_shake.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_fix_shake))
        sys.stderr.write('WARNING: \"fix shake\" style command(s) applied to selected atoms.\n'
                         '         Please check to make sure that the fix group(s) are defined correctly,\n'

                         '         and also check that the atom, bond, and angle types are correct.\n'
                         '######################################################\n')
        assert(non_empty_output)

    if len(l_in_fix_poems) > 0:
        no_warnings = False
        l_in_fix_poems.insert(0, (' '*cindent)+'write(\"'+in_settings+'\") {\n')
        l_in_fix_poems.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_fix_poems))
        sys.stderr.write('WARNING: \"fix poems\" style command(s) applied to selected atoms.\n'
                         '         Please make sure that the fix group(s) are defined correctly.\n'
                         '######################################################\n')
        assert(non_empty_output)

    if len(l_in_fix_qeq) > 0:
        no_warnings = False
        l_in_fix_qeq.insert(0, (' '*cindent)+'write(\"'+in_settings+'\") {\n')
        l_in_fix_qeq.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_fix_qeq))
        sys.stderr.write('WARNING: \"fix qeq\" style command(s) applied to selected atoms.\n'
                         '         Please make sure that the fix group(s) are defined correctly.\n'
                         '######################################################\n')
        assert(non_empty_output)

    if len(l_in_fix_qmmm) > 0:
        no_warnings = False
        l_in_fix_qmmm.insert(0, (' '*cindent)+'write(\"'+in_settings+'\") {\n')
        l_in_fix_qmmm.append((' '*cindent)+'}\n')
        sys.stdout.write('\n')
        sys.stdout.write(''.join(l_in_fix_qmmm))
        sys.stderr.write('WARNING: \"fix qmmm\" style command(s) applied to selected atoms.\n'
                         '         Please make sure that the fix group(s) are defined correctly.\n'
                         '######################################################\n')
        assert(non_empty_output)



    if mol_name != '':
        sys.stdout.write('\n} # end of \"'+mol_name+'\" type definition\n')



    #if non_empty_output and no_warnings:
    if non_empty_output:
        sys.stderr.write('WARNING: The '+g_program_name+' script has not been rigorously tested.\n'
                         '         Exotic (many-body) pair-styles and pair-styles with\n'
                         '         unusual syntax (such as hbond/dreiding) are not understood\n'
                         '         by '+g_program_name+' (...although they are supported by moltemplate).\n'
                         '         Please look over the resulting LT file and check for errors.\n'
                         '         Convert any remaining atom, bond, angle, dihedral, or improper id\n'
                         '         or type numbers to the corresponding $ or @-style counter variables.\n'
                         '         Feel free to report any bugs you find. (-Andrew Jewett 2014-10-30)\n')



except (ValueError, InputError) as err:
    sys.stderr.write('\n'+str(err)+'\n')
    sys.exit(-1)
