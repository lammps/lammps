#!/usr/bin/env python

"""
    lttree_postprocess.py
    This is a stand-alone python script which checks the files created by
    lttree.py to insure that the standard instance-variables ($variables)
    have all been defined.  This script performs a task which is very similar
    to the task performed by lttree_check.py.  This script attempts to detect
    mistakes in the names of $atom, $bond, $angle, $dihedral, $improper, & $mol
    variables.

"""

import sys
from lttree_styles import *
from ttree_lex import ExtractCatName

g_program_name = __file__.split('/')[-1]  # = 'lttree_postprocess.py'
g_version_str  = '0.4'
g_date_str     = '2012-12-12'
atom_style = 'full'
ttree_assignments_fname = 'ttree_assignments.txt'
defined_mols = set([])
defined_atoms = set([])
defined_bonds = set([])
defined_angles = set([])
defined_dihedrals = set([])
defined_impropers = set([])

g_no_check_msg = \
  '(To override this error, run moltemplate using the \"-nocheck\" argument.)\n'

if len(sys.argv) > 1:
    for i in range(0,len(sys.argv)):
        if ((sys.argv[i].lower() == '-atomstyle') or 
            (sys.argv[i].lower() == '-atom-style') or 
            (sys.argv[i].lower() == '-atom_style')):
            if i+1 >= len(sys.argv):
                raise InputError('Error('+g_program_name+'): The '+sys.argv[i]+' flag should be followed by a LAMMPS\n'
                                 '       atom_style name (or single quoted string containing a space-separated\n'
                                 '       list of column names such as: atom-ID atom-type q x y z molecule-ID.)\n')

            atom_style = sys.argv[i+1]
        elif ((sys.argv[i].lower() == '-ttreeassignments') or
              (sys.argv[i].lower() == '-ttree-assignments') or
              (sys.argv[i].lower() == '-ttree_assignments')):
            if i+1 >= len(sys.argv):
                raise InputError('Error('+g_program_name+'): The '+sys.argv[i]+' flag should be followed by \n'
                                 '       a file containing the variable bindings created by ttree/moltemplate.\n')
            ttree_assignments_fname = sys.argv[i+1]
        else:
            pass # ignore other arguments (they are intended for lttree.py)


atom_column_names = AtomStyle2ColNames(atom_style)
i_atomid = 0
i_molid = -1
for i in range(0,len(atom_column_names)):
    if atom_column_names[i].lower() == 'atom-id':
        i_atomid = i
    elif atom_column_names[i].lower() == 'molecule-id':
        i_molid = i
i_max_column = max(i_atomid, i_molid)


# The following variables are defined in "lttree_styles.py"
#data_atoms="Data Atoms"
#data_masses="Data Masses"
#data_velocities="Data Velocities"
#data_bonds="Data Bonds"
#data_angles="Data Angles"
#data_dihedrals="Data Dihedrals"
#data_impropers="Data Impropers"


sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+'\n')

try:
    # ------------ defined_atoms ------------
    try:
        f = open(data_atoms+'.template', 'r')
    except:
        raise InputError('Error('+g_program_name+'): Unable to open file\n'+
                         '\"'+data_atoms+'.template\"\n'
                         '       for reading.  (Do your files lack a \"'+data_atoms+'\" section?)\n'
                         +g_no_check_msg+'\n')

    for line_orig in f:
        ic = line_orig.find('#')
        if ic != -1:
            line = line_orig[:ic]
        else:
            line = line_orig.rstrip('\n')

        tokens = line.strip().split()
        if len(tokens) == 0:
            pass
        elif len(tokens) <= i_max_column:
            raise InputError('Error('+g_program_name+'): The following line from\n'
                             '     "\"'+data_atoms+'.template\" has bad format:\n\n'
                             +line_orig+'\n'
                             '     This my probably an internal error. (Feel free to contact the developer.)\n'
                             +g_no_check_msg+'\n')
        else:
            defined_atoms.add(tokens[i_atomid])
            if i_molid != -1:
                defined_mols.add(tokens[i_molid])

    f.close()


    # ------------ defined_bonds ------------
    try:
        f = open(data_bonds+'.template', 'r')

        for line_orig in f:
            ic = line_orig.find('#')
            if ic != -1:
                line = line_orig[:ic]
            else:
                line = line_orig.rstrip('\n')

            tokens = line.strip().split()
            if len(tokens) == 0:
                pass
            elif len(tokens) < 4:
                raise InputError('Error('+g_program_name+'): The following line from\n'
                                 '     "\"'+data_bonds+'.template\" has bad format:\n\n'
                                 +line_orig+'\n'
                                 '     This my probably an internal error. (Feel free to contact the developer.)\n'
                                 +g_no_check_msg+'\n')
            else:
                defined_bonds.add(tokens[0])
        f.close()
    except:
        pass # Defining bonds (stored in the data_bonds file) is optional


    # ------------ defined_angles ------------
    try:
        f = open(data_angles+'.template', 'r')
        for line_orig in f:
            ic = line_orig.find('#')
            if ic != -1:
                line = line_orig[:ic]
            else:
                line = line_orig.rstrip('\n')

            tokens = line.strip().split()
            if len(tokens) == 0:
                pass
            elif len(tokens) < 5:
                raise InputError('Error('+g_program_name+'): The following line from\n'
                                 '     "\"'+data_angles+'.template\" has bad format:\n\n'
                                 +line_orig+'\n'
                                 '     This my probably an internal error. (Feel free to contact the developer.)\n'
                                 +g_no_check_msg+'\n')
            else:
                defined_angles.add(tokens[0])
        f.close()
    except:
        pass # Defining angles (stored in the data_angles file) is optional


    # ------------ defined_dihedrals ------------
    try:
        f = open(data_dihedrals+'.template', 'r')
        for line_orig in f:
            ic = line_orig.find('#')
            if ic != -1:
                line = line_orig[:ic]
            else:
                line = line_orig.rstrip('\n')

            tokens = line.strip().split()
            if len(tokens) == 0:
                pass
            elif len(tokens) < 6:
                raise InputError('Error('+g_program_name+'): The following line from\n'
                                 '     "\"'+data_dihedrals+'.template\" has bad format:\n\n'
                                 +line_orig+'\n'
                                 '     This my probably an internal error. (Feel free to contact the developer.)\n'
                                 +g_no_check_msg+'\n')
            else:
                defined_dihedrals.add(tokens[0])
        f.close()
    except:
        pass #Defining dihedrals (stored in the data_dihedrals file) is optional


    # ------------ defined_impropers ------------
    try:
        f = open(data_impropers+'.template', 'r')

        for line_orig in f:
            ic = line_orig.find('#')
            if ic != -1:
                line = line_orig[:ic]
            else:
                line = line_orig.rstrip('\n')

            tokens = line.strip().split()
            if len(tokens) == 0:
                pass
            elif len(tokens) < 6:
                raise InputError('Error('+g_program_name+'): The following line from\n'
                                 '     "\"'+data_impropers+'.template\" has bad format:\n\n'
                                 +line_orig+'\n'
                                 '     This my probably an internal error. (Feel free to contact the developer.)\n'
                                 +g_no_check_msg+'\n')
            else:
                defined_impropers.add(tokens[0])
        f.close()
    except:
        pass #Defining impropers (stored in the data_impropers file) is optional




    # ---- Check ttree_assignments to make sure variables are defined ----

    try:
        f = open(ttree_assignments_fname, 'r')
    except:
        raise InputError('Error('+g_program_name+'): Unable to open file\n'+
                         '\"'+ttree_assignments_fname+'\"\n'
                         '       for reading.  (Do your files lack a \"'+data_atoms+'\" section?)\n'
                         +g_no_check_msg+'\n')

    for line_orig in f:

        ic = line_orig.find('#')
        if ic != -1:
            line = line_orig[:ic]
            usage_location_str = 'near ' + line_orig[ic+1:]
        else:
            line = line_orig.rstrip('\n')
            usage_location_str = ''
        
        tokens = line.strip().split()
        if len(tokens) == 0:
            pass
        if len(tokens) > 0:
            # This file contains a list of variables of the form:
            #
            # @/atom:MoleculeType1:C    1
            # @/atom:MoleculeType1:H    2
            # @/atom:MoleculeType2:N    3
            # $/atom:molecule1:N1    1
            # $/atom:molecule1:C1    2
            #   :
            # $/atom:molecule1141:CH    13578
            # $/atom:molecule1142:N3    13579
            #   :
            # We only care about instance variables (which use the '$' prefix)
            # Lines corresponding to static variables (which use the '@' prefix)
            # are ignored during this pass.
            i_prefix = tokens[0].find('$')
            if i_prefix != -1:
                descr_str = tokens[0][i_prefix+1:]
                cat_name = ExtractCatName(descr_str)

                if ((cat_name == 'atom') and
                    (tokens[0] not in defined_atoms)):
                    raise InputError('Error('+g_program_name+'): '+usage_location_str+'\n'+
                                     '      Reference to undefined $atom:\n\n'
                                     '            '+tokens[0]+'     (<--full name)\n\n'+
                                     '      (If that atom belongs to a molecule (or other subunit), make sure that\n'+
                                     '       you specified the correct path which leads to it (using / and ..))\n\n'+
                                     g_no_check_msg)

                elif ((cat_name == 'bond') and
                    (tokens[0] not in defined_bonds)):
                    raise InputError('Error('+g_program_name+'): '+usage_location_str+'\n'+
                                     '      Reference to undefined $bond:\n\n'
                                     '            '+tokens[0]+'     (<--full name)\n\n'+
                                     '      (If that bond belongs to a molecule (or other subunit), make sure that\n'+
                                     '       you specified the correct path which leads to it (using / and ..))\n\n'+
                                     g_no_check_msg)

                elif ((cat_name == 'angle') and
                    (tokens[0] not in defined_angles)):
                    raise InputError('Error('+g_program_name+'): '+usage_location_str+'\n'+
                                     '      Reference to undefined $angle:\n\n'+
                                     '            '+tokens[0]+'     (<--full name)\n\n'+
                                     '      (If that angle belongs to a molecule (or other subunit), make sure that\n'+
                                     '       you specified the correct path which leads to it (using / and ..))\n\n'+
                                     g_no_check_msg)

                elif ((cat_name == 'dihedral') and
                    (tokens[0] not in defined_dihedrals)):
                    raise InputError('Error('+g_program_name+'): '+usage_location_str+'\n\n'+
                                     '      Reference to undefined $dihedral:\n\n'
                                     '            '+tokens[0]+'     (<--full name)\n\n'+
                                     '    (If that dihedral belongs to a molecule (or other subunit), make sure that\n'+
                                     '     you specified the correct path which leads to it (using / and ..))\n\n'+
                                     g_no_check_msg)

                elif ((cat_name == 'improper') and
                    (tokens[0] not in defined_impropers)):
                    raise InputError('Error('+g_program_name+'): '+usage_location_str+'\n'+
                                     '      Reference to undefined $improper:\n\n'
                                     '            '+tokens[0]+'     (<--full name)\n\n'+
                                     '    (If that improper belongs to a molecule (or other subunit), make sure that\n'+
                                     '     you specified the correct path which leads to it (using / and ..))\n\n'+
                                     g_no_check_msg)

                # I used to generate an error when a users defines a $mol 
                # variable but does not associate any atoms with it (or if the
                # user systematically deletes all the atoms in that molecule),
                # but I stopped this practice.
                # I don't think there is any real need to complain if some
                # molecule id numbers are undefined.  LAMMPS does not care.
                #
                #elif ((cat_name == 'mol') and
                #    (tokens[0] not in defined_mols)):
                #    raise InputError('Error('+g_program_name+'): '+usage_location_str+'\n'+
                #                     '      Reference to undefined $mol (molecule-ID) variable:\n\n'
                #                     '            '+tokens[0]+'     (<--full name)\n\n'+
                #                     '    (If that molecule is part of a larger molecule, then make sure that\n'+
                #                     '     you specified the correct path which leads to it (using / and ..))\n\n'+
                #                     g_no_check_msg)


    f.close()

    sys.stderr.write(g_program_name+': -- No errors detected. --\n')
    exit(0)



except (ValueError, InputError) as err:
    sys.stderr.write('\n'+str(err)+'\n')
    sys.exit(1)

