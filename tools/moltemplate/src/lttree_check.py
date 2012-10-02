#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2011, Regents of the University of California
# All rights reserved.

"""
ltree_check.py

The original lttree-file format supports any variables or file names.
However if you plan to process lttree-files using lttree.py, then 
lttree.py requires that variables and file names obey certain naming 
conventions.
This program makes an attempt to check that the variables and file names
which appear in an "lttree" file are not mispelled (or miscapitlised).
(This is not the pretiest code I've ever written.)

"""


import sys
#from ttree import *
from lttree_styles import *
from lttree import *
from ttree_lex import InputError

if sys.version < '2.6':
    raise InputError('Error: Alas, you must upgrade to a newever version of python.')


def CheckVarNames(prefix, descr_str, suffix, srcloc):
    """ Check the name of variables in a lttree-file to confirm 
        that they follow the conventions used by lttree.  
        Almost any variable/category name is permitted, except for 
        names which closely match those reserved by lttree. 

    """

    cat_name, cat_ptkns, leaf_ptkns = \
        DescrToCatLeafPtkns(descr_str,
                            srcloc)

    if (cat_name.lower()=='mol'):
        if (cat_name != 'mol'):
              raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                               'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                               'closely with a reserved lttree variable category.\n'
                               'Perhaps you meant \"mol\"?')

    elif (cat_name.lower()=='group'):
        if (cat_name != 'group'):
              raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                               'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                               'closely with a reserved lttree variable category.\n'
                               'Perhaps you meant \"group\"?')
    elif (cat_name.lower()=='fix'):
        if (cat_name != 'fix'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                             'closely with a reserved lttree variable category.\n'
                             'Use \"fix\" instead.')
    elif (cat_name.lower()=='atom'):
        if (cat_name != 'atom'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Illegal lttree variable category: \"'+cat_name+'\"\n'+
                             'Use \"atom\" instead.')
    elif (cat_name.lower()=='bond'):
        if (cat_name != 'bond'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                             'closely with a reserved lttree variable category.\n'
                             'Use \"bond\" instead.')
    elif (cat_name.lower()=='angle'):
        if (cat_name != 'angle'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                             'closely with a reserved lttree variable category.\n'
                             'Use \"angle\" instead.')
    elif (cat_name.lower()=='dihedral'):
        if (cat_name != 'dihedral'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                             'closely with a reserved lttree variable category.\n'
                             'Use \"dihedral\" instead.')
    elif (cat_name.lower()=='improper'):
        if (cat_name != 'improper'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Variable category: \"'+cat_name+'\" does not match, yet overlaps\n'+
                             'closely with a reserved lttree variable category.\n'
                             'Use \"improper\" instead.')
    else:
        sys.stderr.write('-----------------------------------------------------\n'+
                         'WARNING: in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                         '  Unrecognised template variable category: \"'+cat_name+'\"\n'+
                         '-----------------------------------------------------\n')






def CheckDataFileNames(file_name,
                       srcloc,
                       write_command,
                       fnames_found):
    N_data_prefix = len(data_prefix)
    #data_prefix_no_space = data_prefix.rstrip()
    N_data_prefix_no_space = len(data_prefix)

    section_name = file_name[N_data_prefix:]

    if ((section_name.lower() == 'atom') or
        (section_name.lower() == 'atoms')):
        if (file_name != data_atoms):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_atoms+'\"?')

        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')


    elif ((section_name.lower() == 'velocities') or
          (section_name.lower() == 'velocity')):
        if (file_name != data_velocities):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_velocities+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')


    elif ((section_name.lower() == 'mass') or 
          (section_name.lower() == 'masses')):
        if (file_name != data_masses):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_masses+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'ellipsoids') or 
          (section_name.lower() == 'ellipsoid') or 
          (section_name.lower() == 'elipsoids') or 
          (section_name.lower() == 'elipsoid')):
        if (file_name != data_ellipsoids):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_ellipsoids+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'triangle') or 
          (section_name.lower() == 'triangles')):
        if (file_name != data_triangles):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_triangles+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'line') or 
          (section_name.lower() == 'lines')):
        if (file_name != data_lines):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_lines+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('pair coef') == 0) or
          (section_name.lower().find('pair_coef') == 0) or
          (section_name.lower().find('paircoef') == 0) or
          (section_name.lower().find('pair by type') == 0) or
          (section_name.lower().find('pair bytype') == 0) or
          (section_name.lower().find('pair_by_type') == 0) or
          (section_name.lower().find('pair_bytype') == 0) or
          (section_name.lower().find('pairbytype') == 0)):
        if (file_name != data_pair_coeffs):
            err_msg = 'Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+\
                      'Output file name (\"'+file_name+'\") does not match,\n'+\
                      'yet overlaps closely with reserved lttree-file name.\n'+\
                      'Perhaps you meant \"'+data_pair_coeffs+'\"?'
            if ((section_name.lower().find('by type') != -1) or
                (section_name.lower().find('by_type') != -1) or
                (section_name.lower().find('bytype') != -1)):
                err_msg += '\n    (Note: "pair" parameters are always assigned by type.\n'+\
                           '     There\'s no need to specify \"by type\")'
            raise InputError(err_msg)
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('bond coef') == 0) or
          (section_name.lower().find('bond_coef') == 0) or
          (section_name.lower().find('bondcoef') == 0)):
        if (file_name != data_bond_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_bond_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('angle coef') == 0) or
          (section_name.lower().find('angle_coef') == 0) or
          (section_name.lower().find('anglecoef') == 0)):
        if (file_name != data_angle_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_angle_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('dihedral coef') == 0) or
          (section_name.lower().find('dihedral_coef') == 0) or
          (section_name.lower().find('dihedralcoef') == 0)):
        if (file_name != data_dihedral_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_dihedral_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('improper coef') == 0) or
          (section_name.lower().find('improper_coef') == 0) or
          (section_name.lower().find('impropercoef') == 0)):
        if (file_name != data_improper_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_improper_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')


    # -- class2 data sections --
    elif ((section_name.lower().find('bondbond coef') == 0) or
          (section_name.lower().find('bondbond_coef') == 0) or
          (section_name.lower().find('bondbondcoef') == 0)):
        if (file_name != data_bondbond_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_bondbond_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('bondangle coef') == 0) or
          (section_name.lower().find('bondangle_coef') == 0) or
          (section_name.lower().find('bondanglecoef') == 0)):
        if (file_name != data_bondangle_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_bondangle_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('middlebondtorsion coef') == 0) or
          (section_name.lower().find('middlebondtorsion_coef') == 0) or
          (section_name.lower().find('middlebondtorsioncoef') == 0) or
          (section_name.lower().find('middlebondtorision coef') == 0) or
          (section_name.lower().find('middlebondtorision_coef') == 0) or
          (section_name.lower().find('middlebondtorisioncoef') == 0)):
        if (file_name != data_middlebondtorsion_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_middlebondtorsion_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('endbondtorsion coef') == 0) or
          (section_name.lower().find('endbondtorsion_coef') == 0) or
          (section_name.lower().find('endbondtorsioncoef') == 0) or
          (section_name.lower().find('endbondtorision coef') == 0) or
          (section_name.lower().find('endbondtorision_coef') == 0) or
          (section_name.lower().find('endbondtorisioncoef') == 0)):
        if (file_name != data_endbondtorsion_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_endbondtorsion_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('angletorsion coef') == 0) or
          (section_name.lower().find('angletorsion_coef') == 0) or
          (section_name.lower().find('angletorsioncoef') == 0) or
          (section_name.lower().find('angletorision coef') == 0) or
          (section_name.lower().find('angletorision_coef') == 0) or
          (section_name.lower().find('angletorisioncoef') == 0)):
        if (file_name != data_angletorsion_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_angletorsion_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('angleangletorsion coef') == 0) or
          (section_name.lower().find('angleangletorsion_coef') == 0) or
          (section_name.lower().find('angleangletorsioncoef') == 0) or
          (section_name.lower().find('angleangletorision coef') == 0) or
          (section_name.lower().find('angleangletorision_coef') == 0) or
          (section_name.lower().find('angleangletorisioncoef') == 0)):
        if (file_name != data_angleangletorsion_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_angleangletorsion_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('bondbond13 coef') == 0) or
          (section_name.lower().find('bondbond13_coef') == 0) or
          (section_name.lower().find('bondbond13coef') == 0)):
        if (file_name != data_bondbond13_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_bondbond13_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower().find('angleangle coef') == 0) or
          (section_name.lower().find('angleangle_coef') == 0) or
          (section_name.lower().find('angleanglecoef') == 0)):
        if (file_name != data_angleangle_coeffs):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_angleangle_coeffs+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')






    elif ((section_name.lower() == 'angles by type') or
          (section_name.lower() == 'angles bytype') or          
          (section_name.lower() == 'angles_by_type') or
          (section_name.lower() == 'angles_bytype') or
          (section_name.lower() == 'anglesbytype') or
          (section_name.lower() == 'angle by type') or
          (section_name.lower() == 'angle bytype') or          
          (section_name.lower() == 'angle_by_type') or
          (section_name.lower() == 'angle_bytype') or
          (section_name.lower() == 'anglebytype')):
        if (file_name != data_angles_by_type):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_angles_by_type+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'dihedrals by type') or
          (section_name.lower() == 'dihedrals bytype') or          
          (section_name.lower() == 'dihedrals_by_type') or
          (section_name.lower() == 'dihedrals_bytype') or
          (section_name.lower() == 'dihedralsbytype') or
          (section_name.lower() == 'dihedral by type') or
          (section_name.lower() == 'dihedral bytype') or          
          (section_name.lower() == 'dihedral_by_type') or
          (section_name.lower() == 'dihedral_bytype') or
          (section_name.lower() == 'dihedralbytype')):
        if (file_name != data_dihedrals_by_type):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_dihedrals_by_type+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'impropers by type') or
          (section_name.lower() == 'impropers bytype') or          
          (section_name.lower() == 'impropers_by_type') or
          (section_name.lower() == 'impropers_bytype') or
          (section_name.lower() == 'impropersbytype') or
          (section_name.lower() == 'improper by type') or
          (section_name.lower() == 'improper bytype') or          
          (section_name.lower() == 'improper_by_type') or
          (section_name.lower() == 'improper_bytype') or
          (section_name.lower() == 'improperbytype')):
        if (file_name != data_impropers_by_type):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_impropers_by_type+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')




    elif ((section_name.lower() == 'bonds') or
          (section_name.lower() == 'bond')):
        if (file_name != data_bonds):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_bonds+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'angles') or
          (section_name.lower() == 'angle')):
        if (file_name != data_angles):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_angles+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'dihedrals') or
          (section_name.lower() == 'dihedral')):
        if (file_name != data_dihedrals):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_dihedrals+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'impropers') or
          (section_name.lower() == 'improper')):
        if (file_name != data_impropers):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_impropers+'\"?')
        elif (write_command == 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write(\"'+file_name+'\") instead.\n')

    elif ((section_name.lower() == 'box boundaries') or
          (section_name.lower() == 'box boundary') or
          (section_name.lower() == 'boundaries') or
          (section_name.lower() == 'boundary') or
          (section_name.lower() == 'boundary conditions') or
          (section_name.lower() == 'periodic boundaries') or
          (section_name.lower() == 'periodic boundary conditions') or
          (section_name.lower() == 'periodic_boundaries') or
          (section_name.lower() == 'periodic_boundary_conditions') or
          (section_name.lower() == 'pbc')):
        if ((file_name != data_boundary) and
            (file_name != data_pbc)):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+data_boundary+'\"?\n'
                             '(Specify periodic boundary conditions this way.)')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')
        elif (file_name == data_pbc):
            sys.stderr.write('WARNING: write_once(\"'+data_pbc+'\") is depreciated.\n'
                             '     Use write_once(\"'+data_boundary+'\") instead.\n')




def CheckFileNames(file_name,
                   srcloc,
                   write_command,
                   file_names_found):
    """ 
    Check the write() or write_once() statements in a 
    lttree-file to make sure that the files being written 
    follow the conventions used by lttree.  
    Almost any file name is permitted, except for file names 
    which closely match those reserved by lttree. 

    """

    file_names_found.add(file_name)

    N_data_prefix = len(data_prefix)
    #data_prefix_no_space = data_prefix.rstrip()
    N_data_prefix_no_space = len(data_prefix_no_space)

    if ((file_name[:N_data_prefix].lower() == data_prefix.lower()) and
        (file_name[:N_data_prefix] != data_prefix)):
        raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                         'The beginning of output file (\"'+file_name+'\")\n'
                         'does not match yet overlaps closely with a reserved lttree-file name prefix.\n'
                         '(\"'+data_prefix+'\").  Perhaps you meant \"'+data_prefix+filename[N_data_prefix:]+'\"?')

    # check did they forget the space?
    if (file_name[:N_data_prefix_no_space] == data_prefix_no_space):
        if (file_name[:N_data_prefix] == data_prefix):

            CheckDataFileNames(file_name,
                               srcloc,
                               write_command,
                               file_names_found)

        else:
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'The beginning of output file (\"'+file_name+'\")\n'
                             'does not match yet overlaps closely with a reserved lttree-file name prefix.\n'
                             '(\"'+data_prefix+'\").  Perhaps you meant \"'+data_prefix+filename[N_data_prefix_no_space:]+'\"?')


    elif ((file_name.lower() == 'box boundaries') or
          (file_name.lower() == 'box boundary') or
          (file_name.lower() == 'boundaries') or
          (file_name.lower() == 'boundary') or
          (file_name.lower() == 'boundary conditions') or
          (file_name.lower() == 'periodic boundaries') or
          (file_name.lower() == 'periodic boundary conditions') or
          (file_name.lower() == 'periodic_boundaries') or
          (file_name.lower() == 'periodic_boundary_conditions') or
          (file_name.lower() == 'pbc')):
        # In that case (for one thing) they forgot the data_prefix
        raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                         'Output file name (\"'+file_name+'\") does not match,\n'
                         'yet overlaps closely with reserved lttree-file name.\n'
                         'Perhaps you meant \"'+data_boundary+'\"?\n'
                         '(Specify periodic boundary conditions this way.)')



    elif ((file_name.lower() == 'init') or
          (file_name.lower() == 'in init') or
          (file_name.lower() == 'ininit') or
          (file_name.lower() == 'initialize') or
          (file_name.lower() == 'in initialize') or
          (file_name.lower() == 'ininitialize')):
        if (file_name != in_init):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+in_init+'\"?')
        elif (write_command != 'write_once'):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'When using moltemplate.sh to build LAMMPS input files, you probably do not\n'
                             'want to use the '+write_command+'() command with \"'+file_name+'\".\n'
                             'You should probably use write_once(\"'+file_name+'\") instead.\n')

    elif ((file_name.lower() == 'settings') or
          (file_name.lower() == 'in settings') or
          (file_name.lower() == 'insettings')):

        if (file_name != in_settings):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+in_settings+'\"?')

    elif ((file_name.lower() == 'set_coords') or
          (file_name.lower() == 'set coords') or
          (file_name.lower() == 'setcoords') or
          (file_name.lower() == 'in set_coords') or
          (file_name.lower() == 'in set coords') or
          (file_name.lower() == 'in setcoords')):
        if (file_name != in_set_coords):
            raise InputError('Probable typo in '+ErrorLeader(srcloc.infile,srcloc.lineno)+'\n\n'+
                             'Output file name (\"'+file_name+'\") does not match,\n'
                             'yet overlaps closely with reserved lttree-file name.\n'
                             'Perhaps you meant \"'+in_set_coords+'\"?')







def CheckSyntax(lex):

    """ Parse() builds a static tree of StaticObjs by parsing text file.
    -The "lex" argument is afile or input stream which has been converted 
     to a "TemplateLexer" object (similar to the python's built-in shlex lexer).
    """

    fnames_found = set([])
    prematurely_read_token = None

    while True:

        if prematurely_read_token == None:
            command = lex.get_token()
        else:
            command = prematurely_read_token
        prematurely_read_token = None

        #print('Parse(): token = \"'+command+'\", '+lex.error_leader())

        if command == lex.eof:
            #print('Parse(): EOF encountered\n')
            break

        if ((command == 'write') or (command == 'write_once')):
            open_paren  = lex.get_token()

            #print('Parse():     open_paren=\"'+open_paren+'\"')

            if open_paren=='{': 
                # ..then the user neglected to specify the "file_name" file-name
                # argument.  In that case, supply the default, ''.
                # (which is shorthand for the standard out in this case)
                open_curly = open_paren[0]
                open_paren  = ''
                close_paren = ''
                file_name   = ''
                srcloc  = lex.GetSrcLoc()
            else:
                file_name   = lex.get_token()
                if file_name == ')':
                    file_name == ''
                    close_paren = ')'
                else:
                    close_paren = lex.get_token()
                open_curly  = lex.get_token()
                srcloc      = lex.GetSrcLoc()

            if ((open_curly != '{') or 
                ((open_paren == '')  and (close_paren != '')) or
                ((open_paren == '(') and (close_paren != ')'))):
                raise InputError('Error: in '+lex.error_leader()+'\n\n'
                                 'Syntax error at beginning of '+command+' command.')

            file_name = RemoveOuterQuotes(file_name, lex.quotes)
            # The previous line is similar to:
            #file_name = file_name.strip(lex.quotes)


            CheckFileNames(file_name, lex.GetSrcLoc(), command, fnames_found)


            templ_contents = lex.ReadTemplate()


            for entry in templ_contents:
                if (type(entry) is VarRef):

                    CheckVarNames(entry.prefix,
                                  entry.descr_str, 
                                  entry.suffix,
                                  entry.srcloc)


    #if (data_velocities not in fnames_found):
    #    sys.stderr.write('-------------------------------------------------\n'
    #                     'WARNING: \"'+data_velocities+'\" file not found\n'
    #                     '-------------------------------------------------\n')
    #if (data_pair_coeffs not in fnames_found):
    #    sys.stderr.write('-------------------------------------------------\n'
    #                     'WARNING: \"'+data_pair_coeffs+'\" file not found\n'
    #                     '-------------------------------------------------\n')
    if (data_atoms not in fnames_found):
        sys.stderr.write('WARNING: \"'+data_atoms+'\" file not found\n')
    if (data_masses not in fnames_found):
        sys.stderr.write('WARNING: \"'+data_masses+'\" file not found\n')
    #if (data_bonds not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_bonds+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_angles not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_angles+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_dihedrals not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_dihedrals+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_impropers not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_impropers+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_bond_coeffs not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_bond_coeffs+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_angle_coeffs not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_angle_coeffs+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_dihedral_coeffs not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_dihedral_coeffs+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    #if (data_improper_coeffs not in fnames_found):
    #    sys.stderr.write('--------------------------------------------------\n'
    #                     'WARNING: \"'+data_imrpoper_coeffs+'\" file not found\n'
    #                     '--------------------------------------------------\n')
    if (in_init not in fnames_found):
        sys.stderr.write('WARNING: \"'+in_init+'\" file not found\n')
    
    if (in_settings not in fnames_found):
        sys.stderr.write('WARNING: \"'+in_settings+'\" file not found\n')




def LttreeCheckParseArgs(argv, settings):
                 
    LttreeParseArgs(sys.argv, settings)

    if __name__ == "__main__":
        # Instantiate the lexer we will be using.
        #  (The lexer's __init__() function requires an openned file.
        #   Assuming __name__ == "__main__", then the name of that file should
        #   be the last remaining (unprocessed) argument in the argument list.)
        if len(argv) == 1:
            raise InputError('Error: This program requires at least one argument\n'
                             '       the name of a file containing ttree template commands\n')
        elif len(argv) == 2:
            try:
                settings.lex = TemplateLexer(open(argv[1], 'r'), argv[1]) # Parse text from file
            except IOError: 
                sys.stderr.write('Error: unable to open file\n'
                                 '       \"'+argv[1]+'\"\n'
                                 '       for reading.\n')
                sys.exit(1)
            del(argv[1:2])

        else:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"'+arg+'\"' for arg in argv[1:]]
            raise InputError('Syntax Error('+__file__+'):\n\n'
                             '       Unrecognized argument.\n'
                             '         (That or there is some other problem with the argument list.)\n'
                             '       The problem begins with these arguments:\n'
                             '         '+(' '.join(problem_args))+'\n\n'
                             '       (The actual problem may be earlier in the argument list.\n'
                             '       If these arguments are source files, then keep in mind\n'
                             '       that this program can not parse multiple source files.)\n'
                             '       Check the syntax of the entire argument list.\n')





#######  control flow begins here: #######


if __name__ == "__main__":

    g_program_name = 'lttree_check.py'
    g_version_str  = '0.23'
    g_date_str     = '2012-7-31'
    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

    try:
        # Parse the argument list and instantiate the lexer we will be using:

        #settings = BasicUISettings()
        #BasicUIParseArgs(sys.argv, settings)
        settings = LttreeSettings()
        LttreeCheckParseArgs(sys.argv, settings)

        #Invoke syntax checker:
        CheckSyntax(settings.lex);

        # If no errors so far, then print an encouraging message.
        #sys.stdout.write(g_program_name+': No obvious lttree syntax errors\nwere detected in file \"'+lex.infile+'\".\n')
        exit(0)
        # Alternatley, you can use lex.error_leader() in place of lex.infile().

    except (ValueError, InputError) as err:
        # (note to self: I need to learn more about exceptions.  
        #  some examples at: http://www.doughellmann.com/PyMOTW/shlex/)
        sys.stderr.write('\n'+str(err)+'\n')
        sys.exit(1)

