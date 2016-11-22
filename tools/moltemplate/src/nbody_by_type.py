#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2011, Regents of the University of California
# All rights reserved.

man_page_text = """

    nbody_by_type.py reads a LAMMPS data file (or an excerpt of a LAMMPS)
    data file containing bonded many-body interactions by atom type
    (and bond type), and generates a list of additional interactions
    in LAMMPS format consistent with those type (to the standard out).

    Typical Usage:

    nbody_by_type.py X < old.data > new.data

    --or--

    nbody_by_type.py X \\
                     -atoms atoms.data \\
                     -bonds bonds.data \\
                     -nbody X.data \\
                     -nbodybytype X_by_type.data \\
                     > new_X.data

    In both cases "X" denotes the interaction type, which 
    is either "Angles", "Dihedrals", or "Impropers".
    (Support for other interaction types can be added by the user. See below.)

    -------- Example 1 -------

    nbody_by_type.py X < old.data > new.data

    In this example, nbody_by_type.py reads a LAMMPS data file 
    "orig.data", and extracts the relevant section ("Angles", 
    "Dihedrals", or "Impropers").  It also looks a section named "X By Type",
       (eg. "Angles By type", "Impropers By type", "Impropers By type")
    which contains a list of criteria for automatically defining additional 
    interactions of that type.  For example, this file might contain:

    Angle By Type

    7 1 2 1 * *
    8 2 2 * * *
    9 3 4 3 * *

    The first column is an interaction type ID.
    The next 3 columns are atom type identifiers.
    The final 2 columns are bond type identifiers.
    The * is a wildcard symbol indicating there is no preference for bond types
    in this example.  (Optionally, regular expressions can also be used to
    define a type match, by enclosing the atom or bond type in / slashes.)

        The first line tells us to that there should be a 3-body "Angle" 
    interaction of type "7" whenever an atom of type 1 is bonded to an atom
    of type "2", which is bonded to another atom of type "1" again.
    The second line tells us that an angle is defined whenever three atoms 
    are bonded together and the first two are of type "2".
    (Redundant angle interactions are filtered.)

        New interactions are created for every group of bonded 
    atoms which match these criteria if they are bonded together 
    in the relevant way for that interaction type (as determined by
    nbody_X.py), and printed to the standard output.  For example, 
    suppose you are automatically generating 3-body "Angle" interactions using:

    nbody_by_type Angles < old.data > new.data

    The file "new.data" will be identical to "old.data", however the
    "Angles By Type" section will be deleted, and the following lines of
    text will be added to the "Angles" section:

    394 7 5983 5894 5895
    395 7 5984 5895 5896
    396 7 5985 5896 5897
     :  :   :    :    :
    847 9 14827 14848 14849

    The numbers in the first column are counters which assign a ID to 
    every interaction of that type, and start where the original "Angles"
    data left off (New angle ID numbers do not overlap with old ID numbers).
    The text in the second column ("7", "9", ...) matches the text from the 
    first column of the "Angle By Type" section of the input file.

    -------- Example 2 -------

    nbody_by_type.py X \\
                     -atoms atoms.data \\
                     -bonds bonds.data \\
                     -nbody X.data \\
                     -nbodybytype X_by_type.data \\
                     -prefix "SOMESTRING" -suffix "ANOTHERSTRING" \\
                     > new_X.data

    In particular, for Angle interactions:

    nbody_by_type.py Angles \\
                     -atoms atoms.data \\
                     -bonds bonds.data \\
                     -nbody angles.data \\
                     -nbodybytype angles_by_type.data \\
                     > new_Angles.data

    When run this way, nbody_by_type.py behaves exactly the same way
    as in Example 1, however only the lines of text corresponding to
    the new generated interactions are printed, (not the entire data file).
    Also note, that when run this way, nbody_by_type.py does not read the
    LAMMPS data from the standard input.  Instead, it reads each section of
    the data file from a different file indicated by the arguments following
    the "-atoms", "-bonds", "-nbody", and "-nbodybytype" flags.

    "Angles" is a 3-body interaction style.  So when run this way, 
    nbody_by_type.py will create a 5 (=3+2) column file (new_Angles.data).

Note: the atom, bond and other IDs/types in need not be integers.

Note: This program must be distributed with several python modules, including:
        nbody_Angles.py, nbody_Dihedrals.py, and nbody_Impropers.py.  These
      contain bond definitions for angular, dihedral, and improper interactions.
      (In case any new interaction types are ever added to LAMMPS, with only
      a few lines of python it is easy to edit to define new bonded 
      interaction types by supplying new "nbody_X.py" python module.
      Refer to the modules listed above for examples.)

Note: Optional "-prefix" and "-suffix" arguments can be included to decorate
      the interaction IDs (first column).  For example, -prefix "auto_" and
      -suffix "_angle", causes "new_Angles.data" to contain the following text:

    auto_394_angle 7 5983 5894 5895
    auto_395_angle 7 5984 5895 5896
    auto_396_angle 7 5985 5896 5897
          :        :   :    :    :
    auto_847_angle 9 14827 14848 14849

"""


import sys
from extract_lammps_data import *
from nbody_by_type_lib import GenInteractions_str
from ttree_lex import *
from lttree_styles import AtomStyle2ColNames, ColNames2AidAtypeMolid
import os, inspect # <- Needed to import modules in subdirectories (see below)

if sys.version < '2.6':
    raise InputError('Error: Using python '+sys.version+'\n'
                     '       Alas, you must upgrade to a newer version of python (2.6 or later).')
elif sys.version < '2.7':
    sys.stderr.write('--------------------------------------------------------\n'
                     '----------------- WARNING: OLD PYTHON VERSION ----------\n'
                     '  This program is untested on your python version ('+sys.version+').\n'
                     '  PLEASE LET ME KNOW IF THIS PROGRAM CRASHES (and upgrade python).\n'
                     '    -Andrew   2013-10-25\n'
                     '--------------------------------------------------------\n'
                     '--------------------------------------------------------\n')
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict



def GenInteractions_lines(lines_atoms,
                          lines_bonds,
                          lines_nbody,
                          lines_nbodybytype,
                          atom_style,
                          g_bond_pattern,
                          canonical_order, #function to sort atoms and bonds
                          prefix='',
                          suffix='',
                          report_progress = False):

    column_names = AtomStyle2ColNames(atom_style)
    i_atomid, i_atomtype, i_molid = ColNames2AidAtypeMolid(column_names)

    atomids_str = []
    atomtypes_str = []

    for iv in range(0, len(lines_atoms)):
        line = lines_atoms[iv].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) > 0:
            tokens = SplitQuotedString(line)
            if ((len(tokens) <= i_atomid) or (len(tokens) <= i_atomtype)):
                raise(InputError('Error not enough columns on line '+str(iv+1)+' of \"Atoms\" section.'))
            tokens = SplitQuotedString(line)
            atomids_str.append(EscCharStrToChar(tokens[i_atomid]))
            atomtypes_str.append(EscCharStrToChar(tokens[i_atomtype]))

    bondids_str = []
    bondtypes_str = []
    bond_pairs = []
    for ie in range(0, len(lines_bonds)):
        line = lines_bonds[ie].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) > 0:
            tokens = SplitQuotedString(line)
            if len(tokens) < 4:
                raise(InputError('Error not enough columns on line '+str(ie+1)+' of \"Bonds\" section.'))
            bondids_str.append(EscCharStrToChar(tokens[0]))
            bondtypes_str.append(EscCharStrToChar(tokens[1]))
            bond_pairs.append( (EscCharStrToChar(tokens[2]),
                                EscCharStrToChar(tokens[3])) )

    typepattern_to_coefftypes = []

    for i in range(0, len(lines_nbodybytype)):
        line = lines_nbodybytype[i].strip()
        if '#' in line:
            icomment = line.find('#')
            line = (line[:icomment]).strip()
        if len(line) > 0:
            tokens = SplitQuotedString(line)

            if ((len(tokens) != 1 + g_bond_pattern.GetNumVerts()) and
                (len(tokens) != 1 + g_bond_pattern.GetNumVerts() 
                                  + g_bond_pattern.GetNumEdges())):
                raise(InputError('Error: Wrong number of columns in \"By Type\" section of data file.\n'
                                 'Offending line:\n'+
                                 '\"'+line+'\"\n'
                                 'Expected either '+
                                 str(1 + g_bond_pattern.GetNumVerts()) + ' or ' +
                                 str(1 + g_bond_pattern.GetNumVerts() + 
                                     g_bond_pattern.GetNumEdges())
                                 + ' colunms.'))

            coefftype = EscCharStrToChar(tokens[0])
            typepattern = []

            for typestr in tokens[1:]:
                if ((len(typestr) >= 2) and 
                    (typestr[0] == '/') and (typestr[-1] == '/')):
                    regex_str = typestr[1:-1]
                    typepattern.append( re.compile(regex_str) )
                else:
                    typepattern.append(EscCharStrToChar(typestr))

            # If the user neglected to specify the bond types, assume '*'
            if len(tokens) == 1 + g_bond_pattern.GetNumVerts():
                typepattern += ['*'] * g_bond_pattern.GetNumEdges()

            typepattern_to_coefftypes.append([typepattern, coefftype])

    coefftype_to_atomids_str = GenInteractions_str(bond_pairs,
                                                   g_bond_pattern,
                                                   typepattern_to_coefftypes,
                                                   canonical_order,
                                                   atomids_str,
                                                   atomtypes_str,
                                                   bondids_str,
                                                   bondtypes_str,
                                                   report_progress)
    lines_nbody_new = []
    for coefftype, atomids_list in coefftype_to_atomids_str.items():
        for atomids_found in atomids_list:
            n = len(lines_nbody) + len(lines_nbody_new) + 1
            line = prefix+str(n)+suffix+' '+ \
                coefftype+' '+(' '.join(atomids_found))+'\n'
            lines_nbody_new.append(line)

    return lines_nbody_new



def GenInteractions_files(lines_data,
                          src_bond_pattern,
                          fname_atoms,
                          fname_bonds,
                          fname_nbody,
                          fname_nbodybytype,
                          section_name,
                          section_name_bytype,
                          atom_style,
                          prefix='',
                          suffix='',
                          report_progress = False):

    if fname_atoms == None:
        lines_atoms = [line for line in ExtractDataSection(lines_data, 'Atoms')]
    else:
        try:
            f = open(fname_atoms, 'r')
        except:
            sys.stderr.write('Error: Unable to open file \"'+fname_atoms+'\" for reading.\n')
            sys.exit(-1)
        lines_atoms = [line for line in f.readlines() 
                       if ((len(line.strip())>0) and (line.strip()[0] != '#'))]
        f.close()


    if fname_bonds == None:
        lines_bonds = [line for line in ExtractDataSection(lines_data, 'Bonds')]
    else:
        try:
            f = open(fname_bonds, 'r')
        except IOError: 
            sys.stderr.write('Error: Unable to open file \"'+fname_bonds+'\" for reading.\n')
            sys.exit(-1)
        lines_bonds = [line for line in f.readlines() 
                       if ((len(line.strip())>0) and (line.strip()[0] != '#'))]
        f.close()


    if fname_nbody == None:
        lines_nbody = [line for line in ExtractDataSection(lines_data, section_name)]
    else:
        try:
            f = open(fname_nbody, 'r')
            lines_nbody = [line for line in f.readlines() 
                           if ((len(line.strip())>0) and (line.strip()[0] != '#'))]
            f.close()
        except IOError: 
            #sys.stderr.write('    (omitting optional file \"'+fname_nbody+'\")\n')
            lines_nbody = []


    if fname_nbodybytype == None:
        lines_nbodybytype=[line for 
                           line in ExtractDataSection(lines_data,
                                                      section_name_bytype)]

    else:
        try:
            f = open(fname_nbodybytype, 'r')
        except:
            sys.stderr.write('Error: Unable to open file \"'+fname_nbodybytype+'\" for reading.\n')
            sys.exit(-1)
        lines_nbodybytype = [line for line in f.readlines() 
                             if((len(line.strip())>0)and(line.strip()[0]!='#'))]
        f.close()


    try:
        g = __import__(src_bond_pattern) #defines g.bond_pattern, g.canonical_order
    except:
        # If not found, look for it in the "nbody_alternate_symmetry" directory
        #http://stackoverflow.com/questions/279237/import-a-module-from-a-relative-path
        cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"nbody_alternate_symmetry")))
        if cmd_subfolder not in sys.path:
            sys.path.insert(0, cmd_subfolder)
        try:
            g = __import__(src_bond_pattern) #defines g.bond_pattern, g.canonical_order
        except:
            sys.stderr.write('Error: Unable to locate file \"'+src_bond_pattern+'\"\n'
                             '       (Did you mispell the file name?\n'
                             '        Check the \"nbody_alternate_symmetry/\" directory.)\n')
            sys.exit(-1)
        



    return GenInteractions_lines(lines_atoms,
                                 lines_bonds,
                                 lines_nbody,
                                 lines_nbodybytype,
                                 atom_style,
                                 g.bond_pattern,
                                 g.canonical_order,
                                 prefix,
                                 suffix,
                                 report_progress)






if __name__ == "__main__":

    g_program_name = __file__.split('/')[-1]  # = 'nbody_by_type.py'
    g_date_str     = '2014-12-19'
    g_version_str  = '0.18'

    bond_pattern_module_name = ""

    #######  Main Code Below: #######
    sys.stderr.write(g_program_name+' v'+g_version_str+' '+g_date_str+' ')
    if sys.version < '3':
        sys.stderr.write(' (python version < 3)\n')
    else:
        sys.stderr.write('\n')

    try:

        fname_atoms = None
        fname_bonds = None
        fname_nbody = None
        fname_nbodybytype = None
        atom_style = 'full'
        prefix=''
        suffix=''

        argv = [arg for arg in sys.argv]

        if len(argv) == 1:
            raise InputError('Error: Missing argument required.\n'
                             '       The \"'+g_program_name+'\" program requires an argument containing the\n'
                             '       name of a section from a LAMMPS data file storing bonded interactions.\n'
                             '       (For example: "Angles", "Dihedrals", or "Impropers".)\n'
                             #'        Note: The first letter of each section is usually capitalized.)\n'
                             '\n'
                             '--------------- general documentation -------------\n'
                             '\n' + man_page_text + '\n')

        section_name = ''         # (This will be replaced later.)
        section_name_bytype = ''  # (This will be replaced later.)

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
                    raise InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                                     '       text which might appear in the "Atoms" section of a LAMMPS data file.\n')
                fname_atoms = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-bonds':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                                     '       text which might appear in the "Bonds" section of a LAMMPS data file.\n')
                fname_bonds = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-nbody':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a file name\n')

                    #raise InputError('Error: '+argv[i]+' flag should be followed by a file name containing lines of\n'
                    #                 '       text which might appear in the "'+section_name+' section of a LAMMPS data file.\n')
                fname_nbody = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-nbodybytype':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a file name\n')

                    #raise InputError('Error: '+argv[i]+' flag should be followed by a file name containing\n'
                    #                 '       text which might appear in the "'+section_name+' By Type" section\n'
                    #                 '       of a LAMMPS data file.\n')
                fname_nbodybytype = argv[i+1]
                del(argv[i:i+2])

            elif ((argv[i].lower() == '-atom-style') or 
                (argv[i].lower() == '-atom_style')):
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a an atom_style name.\n'
                                     '       (Or single quoted string which includes a space-separated\n'
                                     '       list of column names.)\n')
                atom_style = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-prefix':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a prefix string\n'
                                     '       (a string you want to appear to the left of the integer\n'
                                     '        which counts the bonded interactions you have generated.)\n')
                prefix = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-suffix':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by a suffix string\n'
                                     '       (a string you want to appear to the right of the integer\n'
                                     '        which counts the bonded interactions you have generated.)\n')
                prefix = argv[i+1]
                del(argv[i:i+2])

            elif argv[i].lower() == '-subgraph':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by the name of a python file\n'
                                     '       containing the definition of the subgraph you are searching for\n'
                                     '       and it\'s symmetry properties.\n'
                                     '       (See nbody_Dihedrals.py for example.)\n')
                bond_pattern_module_name = argv[i+1]
                # If the file name ends in ".py", then strip off this suffix.
                # For some reason, the next line does not work:
                #bond_pattern_module_name=bond_pattern_module_name.rstrip('.py')
                # Do this instead
                pc = bond_pattern_module_name.rfind('.py')
                if pc != -1:
                    bond_pattern_module_name = bond_pattern_module_name[0:pc]

                del(argv[i:i+2])

            elif argv[i].lower() == '-section':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by the name of the LAMMPS\n'
                                     '       Data section describing the type of interaction being generated.\n'
                                     '       (For example: \"Angles\", \"Dihedrals\", \"Impropers\", etc...)\n')
                section_name = argv[i+1]
                del(argv[i:i+2])


            elif argv[i].lower() == '-sectionbytype':
                if i+1 >= len(argv):
                    raise InputError('Error: '+argv[i]+' flag should be followed by the name of the\n'

                                     '       write_once(\"???? By Type\") section describing how to create the\n'
                                     '       interactions.  (For example: \"Angles By Type\", \"Dihedrals By Type\",\n'
                                     '        \"Impropers By Type\", etc...  Note that this argument\n'
                                     '        will contain spaces, so surround it with quotes.)\n')
                
                section_name_bytype = argv[i+1]
                del(argv[i:i+2])

            elif argv[i][0] == '-':
                raise InputError('Error('+g_program_name+'):\n'
                                 'Unrecogized command line argument \"'+argv[i]+'\"\n')
            else:
                i += 1

        #if len(argv) == 1:
        #    raise InputError('Error: Missing argument required.\n'
        #                     '       The \"'+g_program_name+'\" program requires an argument containing the\n'
        #                     '       name of a section from a LAMMPS data file storing bonded interactions.\n'
        #                     '       (For example: "Angles", "Dihedrals", or "Impropers".)\n')
        #                     #'        Note: The first letter of each section is usually capitalized.)\n'

        if len(argv) == 1:
            pass
        elif len(argv) == 2:
            section_name = argv[1]
            section_name_bytype = section_name + ' By Type'
            # default bond_pattern_module name
            if bond_pattern_module_name == "":  #<--if not set by user
                bond_pattern_module_name = 'nbody_'+section_name
            del(argv[1:2])
        else:
            # if there are more than 2 remaining arguments,
            problem_args = ['\"'+arg+'\"' for arg in argv[1:]]
            raise InputError('Syntax Error('+g_program_name+'):\n\n'
                             '       Problem with argument list.\n'
                             '       The remaining arguments are:\n\n'
                             '         '+(' '.join(problem_args))+'\n\n'
                             '       (The actual problem may be earlier in the argument list.)\n')

        if ((section_name == '') or
            (section_name_bytype == '') or
            (bond_pattern_module_name == '')):
            raise InputError('Syntax Error('+g_program_name+'):\n\n'
                             '       You have not defined the following arguments:\n'
                             '       -section name\n'
                             '       -sectionbytype namebytype\n'
                             '       -subgraph pythonfile.py\n')

        # ------------ Done parsing argument list ----------

        if (fname_atoms or fname_bonds or fname_nbody or fname_nbodybytype):
            output_full_DATA_file = False
            lines_data = []
        else:
            output_full_DATA_file = True
            lines_data = sys.stdin.readlines()

        # Calculate the interactions and generate a list of lines of text

        lines_new_interactions = \
            GenInteractions_files(lines_data,
                                  bond_pattern_module_name,
                                  fname_atoms,
                                  fname_bonds,
                                  fname_nbody,
                                  fname_nbodybytype,
                                  section_name,
                                  section_name_bytype,
                                  atom_style,
                                  prefix,
                                  suffix,
                                  report_progress=True)

        # Print this text to the standard out.

        # Question: Do we write out the entire DATA file, 
        # or just the portion that was generated by this program?

        if not output_full_DATA_file:
            # ...then only print out the interactions which were generated
            # by this program, omitting any lines from the original data file:

            # (This is the way I usually run this program.)
            for line in lines_new_interactions:
                sys.stdout.write(line)


        else:


            # ...then print out the entire data file, deleting the "By Type" 
            # section, and adding the generated lines of text to the corresponding

            # If present, update the interaction counter at the beginning 
            # of the LAMMPS data file.  (For example, if if 100 new "Angles" 
            # interactions were generated, replace "2 Angles" with "102 Angles")
            # 
            for i in range(0, len(lines_data)):
                line = lines_data[i].strip()
                tokens = SplitQuotedString(line)

                # updating the interaction counter
                if ((len(tokens) == 2) and (tokens[1] == (section_name).lower())):
                    tokens[0] = str(int(tokens[0]) + len(lines_new_interactions))
                    lines_data[i] = ' '.join(tokens) + '\n'

                # stop when you come to a section header
                elif line in lammps_data_sections:
                    #"lammps_data_sections" is defined in "extract_lammps_data.py"
                    break


            # locate the appropriate section of the data file
            # (storing the type of interactions we just created)
            i_nbody_a, i_nbody_b = \
                FindDataSection(lines_data, section_name)

            if i_nbody_a == -1:
                if len(lines_new_interactions) > 0:
                    # If not found, create a new section at the end of the file,
                    # containing a section name followed by the list of lines
                    lines_data += ['\n', section_name+'\n', '\n'] + \
                                   lines_new_interactions + ['\n']
            else:
                # Insert the new lines into the existing section
                lines_data[i_nbody_b:i_nbody_b] = lines_new_interactions

            # Figure out where the "By Type" section is located 
            # (so we skip over it)
            i_bytype_a, i_bytype_b = \
                FindDataSection(lines_data, section_name_bytype)

            in_bytype_section = False
            for i in range(0, len(lines_data)):
                line = lines_data[i].strip()
                # Omit all lines of text in the 'By Type' section (including the 
                # header and commments or blank lines which immediately follow it.)
                if line == section_name_bytype:
                    in_bytype_section = True
                elif i == i_bytype_b:
                    in_bytype_section = False

                if not in_bytype_section:
                    sys.stdout.write(lines_data[i])

    except (ValueError, InputError) as err:
        sys.stderr.write('\n'+str(err)+'\n')
        sys.exit(-1)

