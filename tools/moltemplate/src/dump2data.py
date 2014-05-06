#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dump2data.py

Extract dynamical degrees of freedom from a lammps DUMP file (from the stdin)
and construct a new DATA file (to the stdout).
A reference DATA file is needed (argument).

   basic usage 
./dump2data.py orig_file.data < dump.lammpstrj > new_file.data
   (This extract last frame, uses "full" atom_style.)

   options:
./dump2data.py [-t t -atomstyle style] orig.data < dump.lammpstrj > new.data

"""

# Authors: Andrew Jewett
# License: New BSD License
# Copyright (c) 2014
# All rights reserved.


import sys
from collections import defaultdict
from operator import itemgetter, attrgetter


class InputError(Exception):
    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg


def ErrorLeader(infile, lineno):
    return '\"'+infile+'\", line '+str(lineno)+': '


class MiscSettings(object):
    def __init__(self):
        self.tstart = None
        self.tstop = None
        self.timestep_str = ''
        self.last_frame = False
        self.center_frame = False
        self.output_format = 'data'
        self.input_format = 'dump'
        self.multi = True
        self.skip_interval = 1
        self.scale = None


class AtomStyleSettings(object):
    def __init__(self):
        # The following new member data indicate which columns store
        # LAMMPS-specific information.  
        # The next 6 members store keep track of the different columns 
        # of the "Atoms" section of a LAMMPS data file:
        self.column_names = [] #<--A list of column names (optional)
        self.i_coords=[] #<--A triplet of integers indicating which columns store coordinate data
        #self.ii_coords= [] #<--A list of triplets of column indexes storing coordinate data
        self.ii_vects=[]       #<--A list of triplets of column indexes storing directional data
                               #   (such as dipole or ellipsoid orientations)
        self.i_atomid=None     #<--An integer indicating which column has the atomid
        self.i_atomtype=None   #<--An integer indicating which column has the atomtype
        self.i_molid=None      #<--An integer indicating which column has the molid, if applicable



class DataSettings(AtomStyleSettings):
    def __init__(self):
        AtomStyleSettings.__init__(self)
        self.contents  = ''
        self.file_name = ''



# Atom Styles in LAMMPS as of 2011-7-29
g_style_map = {'angle':    ['atom-ID','molecule-ID','atom-type','x','y','z'],
               'atomic':   ['atom-ID','atom-type','x','y','z'],
               'bond':     ['atom-ID','molecule-ID','atom-type','x','y','z'],
               'charge':   ['atom-ID','atom-type','q','x','y','z'],
               'colloid':  ['atom-ID','atom-type','x','y','z'],
               'dipole':   ['atom-ID','atom-type','q','x','y','z','mux','muy','muz'],
               'electron': ['atom-ID','atom-type','q','spin','eradius','x','y','z'],
               'ellipsoid':['atom-ID','atom-type','x','y','z','quatw','quati','quatj','quatk'],
               'full':     ['atom-ID','molecule-ID','atom-type','q','x','y','z'],
               'granular': ['atom-ID','atom-type','diameter','density','x','y','z'],
               'molecular':['atom-ID','molecule-ID','atom-type','x','y','z'],
               'peri':     ['atom-ID','atom-type','volume','density','x','y','z'],
               'hybrid':   ['atom-ID','atom-type','x','y','z']}




def AtomStyle2ColNames(atom_style_string):

    atom_style_string = atom_style_string.strip()
    if len(atom_style_string) == 0:
        raise InputError('Error(dump2data): Invalid atom_style\n'
                         '       (The atom_style command was followed by an empty string.)\n')
    atom_style_args = atom_style_string.split()
    atom_style      = atom_style_args[0]

    hybrid_args     = atom_style_args[1:]
    if (atom_style not in g_style_map):
        if (len(atom_style_args) >= 2):
            # If the atom_style_string includes at least 2 words, then we 
            # interpret this as a list of the individual column names
            return atom_style_args
        else:
            raise InputError('Error(dump2data): Unrecognized atom_style: \"'+atom_style+'\"\n')

    if (atom_style != 'hybrid'):
        return g_style_map[atom_style]
    else:
        column_names = ['atom-ID','atom-type','x','y','z']
        if (len(hybrid_args)==0):
                raise InputError('Error(dump2data): atom_style hybrid must be followed by a sub_style.\n')
        for sub_style in hybrid_args:
            if (sub_style not in g_style_map):
                raise InputError('Error(dump2data): Unrecognized atom_style: \"'+sub_style+'\"\n')
            for cname in g_style_map[sub_style]:
                if cname not in column_names:
                    column_names.append(cname)

        return column_names


def ColNames2AidAtypeMolid(column_names):
    # Because of the diversity of ways that these 
    # numbers are referred to in the LAMMPS documentation,
    # we have to be flexible and allow the user to refer
    # to these quantities in a variety of ways.
    # Hopefully this covers everything:

    i_atomid = None
    if 'atom-ID' in column_names:
        i_atomid = column_names.index('atom-ID')
    elif 'atom−ID' in column_names: # (− is the character used in the manual)
        i_atomid = column_names.index('atom−ID')
    elif 'atomID' in column_names:
        i_atomid = column_names.index('atomID')
    elif 'atomid' in column_names:
        i_atomid = column_names.index('atomid')
    elif 'id' in column_names:
        i_atomid = column_names.index('id')
    elif 'atom' in column_names:
        i_atomid = column_names.index('atom')
    elif '$atom' in column_names:
        i_atomid = column_names.index('$atom')
    else:
        raise InputError('Error(dump2data): List of column names lacks an \"atom-ID\"\n')

    i_atomtype = None
    if 'atom-type' in column_names:
        i_atomtype = column_names.index('atom-type')
    elif 'atom−type' in column_names: # (− hyphen character used in manual)
        i_atomtype = column_names.index('atom−type')
    elif 'atomtype' in column_names:
        i_atomtype = column_names.index('atomtype')
    elif 'type' in column_names:
        i_atomtype = column_names.index('type')
    elif '@atom' in column_names:
        i_atomtype = column_names.index('@atom')
    else:
        raise InputError('Error(dump2data): List of column names lacks an \"atom-type\"\n')

    i_molid = None
    if 'molecule-ID' in column_names:
        i_molid = column_names.index('molecule-ID')
    elif 'molecule−ID' in column_names: # (− hyphen character used in manual)
        i_molid = column_names.index('molecule−ID')
    elif 'moleculeID' in column_names:
        i_molid = column_names.index('moleculeID')
    elif 'moleculeid' in column_names:
        i_molid = column_names.index('moleculeid')
    elif 'molecule' in column_names:
        i_molid = column_names.index('molecule')
    elif 'molID' in column_names:
        i_molid = column_names.index('molID')
    elif 'molid' in column_names:
        i_molid = column_names.index('molid')
    elif 'mol' in column_names:
        i_molid = column_names.index('mol')
    elif '$mol' in column_names:
        i_molid = column_names.index('$mol')
    else:
        pass # some atom_types do not have a valid molecule-ID

    return i_atomid, i_atomtype, i_molid



def ColNames2Coords(column_names):
    """ Which of the columns correspond to coordinates 
        which must be transformed using rigid-body 
        (affine: rotation + translation) transformations?
        This function outputs a list of lists of triplets of integers.

    """
    i_x = None
    i_y = None
    i_z = None
    if 'x' in column_names:
        i_x = column_names.index('x')
    if 'y' in column_names:
        i_y = column_names.index('y')
    if 'z' in column_names:
        i_z = column_names.index('z')
    if (((i_x != None) != (i_y != None)) or
        ((i_y != None) != (i_z != None)) or
        ((i_z != None) != (i_x != None))):
        raise InputError('Error(dump2data): columns must include \"x\", \"y\", and \"z\".\n')
    return [[i_x, i_y, i_z]]


def ColNames2Vects(column_names):
    """ Which of the columns correspond to coordinates 
        which must be transformed using rotations?
        Some coordinates like dipole moments and 
        ellipsoid orientations should only be rotated
        (not translated).
        This function outputs a list of lists of triplets of integers.

    """
    vects = []
    i_mux = None
    i_muy = None
    i_muz = None
    if 'mux' in column_names:
        i_mux = column_names.index('mux')
    if 'muy' in column_names:
        i_muy = column_names.index('muy')
    if 'muz' in column_names:
        i_muz = column_names.index('muz')
    if (((i_mux != None) != (i_muy != None)) or
        ((i_muy != None) != (i_muz != None)) or
        ((i_muz != None) != (i_mux != None))):
        raise InputError('Error(dump2data): custom atom_style list must define mux, muy, and muz or none.\n')
    if i_mux != None:
        vects.append([i_mux, i_muy, i_muz])
    i_quati = None
    i_quatj = None
    i_quatk = None
    if 'quati' in column_names:
        i_quati = column_names.index('quati')
    if 'quatj' in column_names:
        i_quatj = column_names.index('quatj')
    if 'quatk' in column_names:
        i_quatk = column_names.index('quatk')
    if (((i_quati != None) != (i_quatj != None)) or
        ((i_quatj != None) != (i_quatk != None)) or
        ((i_quatk != None) != (i_quati != None))):
        raise InputError('Error(dump2data): custom atom_style list must define quati, quatj, and quatk or none.\n')
    if i_quati != None:
        vects.append([i_quati, i_quatj, i_quatk])
    return vects





def ParseArgs(argv, 
              misc_settings, 
              data_settings, 
              warning_strings=None):

    # Loop over the remaining arguments not processed yet.
    # These arguments are specific to the lttree.py program
    # and are not understood by this program.
    i = 1
    while i < len(argv):
        #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
        if ((argv[i].lower() == '-atomstyle') or 
            (argv[i].lower() == '-atom_style') or 
            (argv[i].lower() == '-atom-style')):
            in_init = []
            if i+1 >= len(argv):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by a an atom_style name.\n'
                                 '       (Or single quoted string which includes a space-separated\n'
                                 '       list of column names.)\n')
            data_settings.column_names = AtomStyle2ColNames(argv[i+1])
            sys.stderr.write('    \"Atoms\" column format:\n')
            sys.stderr.write('    '+(' '.join(data_settings.column_names))+'\n')

            # ColNames2Coords() and ColNames2Vects() generate lists of
            # triplets of integers, storing the column numbers containing
            # x, y, and z coordinate values, and vx,vy,vz direction vectors.
            data_settings.ii_vects = ColNames2Vects(data_settings.column_names)
            ii_coords = ColNames2Coords(data_settings.column_names)
            # This program assumes that there is only one coordinate triplet
            # (x,y,z) for each atom.  Hence we assume that len(ii_coords)==1
            assert(len(ii_coords) == 1)
            data_settings.i_coords = ii_coords[0]

            # Now figure out which columns correspond to atomid, atomtype, molid
            data_settings.i_atomid, data_settings.i_atomtype, data_settings.i_molid = ColNames2AidAtypeMolid(data_settings.column_names)
            del(argv[i:i+2])

        elif (argv[i].lower() == '-icoord'):
            if i+1 >= len(argv):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by list of integers\n'
                                 '       corresponding to column numbers for coordinates in\n'
                                 '       the \"Atoms\" section of a LAMMPS data file.\n') 
            ilist = argv[i+1].split()
            if (len(ilist) % 3) != 0:
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by list of integers.\n'
                                 '       This is usually a list of 3 intebers, but it can contain more.\n'
                                 '       The number of cooridnate columns must be divisible by 3,\n'
                                 '       (even if the simulation is in 2 dimensions)\n')

            #ii_coords = []
            #for i in range(0, len(ilist)/3):
            #    cols = [ilist[3*i]+1, ilist[3*i+1]+1, ilist[3*i+2]+1]
            #    ii_coords.append(cols)
            #if ((len(ii_coords) != 0) or (len(ii_coords[0]) != 3)):
            #    raise InputError('Error(dump2data): Argument \"'+argv[i]+'\" must be followed by exactly 3 integers.\n')

            data_settings.i_coords = ilist
            if (len(i_coords) != 3):
                raise InputError('Error(dump2data): Argument \"'+argv[i]+'\" must be followed by exactly 3 integers.\n')

            data_settings.i_coords = ii_coords[0]

            del(argv[i:i+2])

        elif (argv[i].lower() == '-ivect'):
            if i+1 >= len(argv):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by list of integers\n'
                                 '       corresponding to column numbers for direction vectors in\n'
                                 '       the \"Atoms\" section of a LAMMPS data file.\n') 
            ilist = argv[i+1].split()
            if (len(ilist) % 3) != 0:
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by list of integers.\n'
                                 '       This is usually a list of 3 intebers, but it can contain more.\n'
                                 '       The number of cooridnate columns must be divisible by 3,\n'
                                 '       (even if the simulation is in 2 dimensions)\n')

            data_settings.ii_vects = []
            for i in range(0, len(ilist)/3):
                cols = [ilist[3*i]+1, ilist[3*i+1]+1, ilist[3*i+2]+1]
                setting.ii_vects.append(cols)
            # This should override any earlier settings as a result of the
            # -atomstyle argument.  So you can specify a custom list of column
            # names using -atomstyle "list of column names", and then afterwards
            # specify which of these columns correspond to direction vectors
            # using the "-ivect" command line argument later on.
            # This way, in theory you should be able to read columns from
            # new custom atom-styles that have not been invented yet.
            # (Although I haven't tested this.)

            del(argv[i:i+2])
        # i_atomid is not really needed for this program, but I load it anyway
        elif ((argv[i].lower() == '-iatomid') or 
              (argv[i].lower() == '-iid') or 
              (argv[i].lower() == '-iatom-id')):
            if ((i+1 >= len(argv)) or (not str.isdigit(argv[i+1]))):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by an integer\n'
                                 '       (>=1) indicating which column in the \"Atoms\" section of a\n'
                                 '       LAMMPS data file contains the atom id number (typically 1).\n'
                                 '       (This argument is unnecessary if you use the -atomstyle argument.)\n')
            i_atomid = int(argv[i+1])-1
            del(argv[i:i+2])
        # i_atomtype is not really needed for this program, but I load it anyway
        elif ((argv[i].lower() == '-iatomtype') or 
              (argv[i].lower() == '-itype') or               
              (argv[i].lower() == '-iatom-type')):
            if ((i+1 >= len(argv)) or (not str.isdigit(argv[i+1]))):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by an integer\n'
                                 '       (>=1) indicating which column in the \"Atoms\" section of a\n'
                                 '       LAMMPS data file contains the atom type.\n'
                                 '       (This argument is unnecessary if you use the -atomstyle argument.)\n')
            i_atomtype = int(argv[i+1])-1
            del(argv[i:i+2])
        # i_molid is not really needed for this program, but I load it anyway
        elif ((argv[i].lower() == '-imolid') or 
              (argv[i].lower() == '-imol') or 
              (argv[i].lower() == '-imol-id') or 
              (argv[i].lower() == '-imoleculeid') or 
              (argv[i].lower() == '-imolecule-id')):
            if ((i+1 >= len(argv)) or (not str.isdigit(argv[i+1]))):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by an integer\n'
                                 '       (>=1) indicating which column in the \"Atoms\" section of a\n'
                                 '       LAMMPS data file contains the molecule id number.\n'
                                 '       (This argument is unnecessary if you use the -atomstyle argument.)\n')
            del(argv[i:i+2])
        # Which frame do we want?
        elif (argv[i].lower() == '-t'):
            if ((i+1 >= len(argv)) or (not str.isdigit(argv[i+1]))):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by an integer indicating\n'
                                 '       the frame you want to extract from the dump file (trajectory).\n'
                                 '       This integer should match the timestep corresponding to the frame\n'
                                 '       whose coordinates you wish to extract.\n')
            misc_settings.timestep_str = argv[i+1]
            del(argv[i:i+2])
            misc_settings.multi = False
            misc_settings.last_frame = False

        elif (argv[i].lower() == '-tstart'):
            if ((i+1 >= len(argv)) or (not str.isdigit(argv[i+1]))):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by an integer indicating\n'
                                 '       the first frame you want to extract from the dump file (trajectory).\n'
                                 '       This integer should match the timestep corresponding to the frame\n'
                                 '       (after which) you wish to extract coordinates.\n')
            misc_settings.tstart = float(argv[i+1])
            del(argv[i:i+2])
            misc_settings.multi = True

        elif (argv[i].lower() == '-tstop'):
            if ((i+1 >= len(argv)) or (not str.isdigit(argv[i+1]))):
                raise InputError('Error(dump2data): '+argv[i]+' flag should be followed by an number indicating\n'
                                 '       the first frame you want to extract from the dump file (trajectory).\n'
                                 '       Frames after this timestep will be ignored.\n')
            misc_settings.tstop = float(argv[i+1])
            del(argv[i:i+2])
            misc_settings.multi = True

        elif (argv[i].lower() == '-center'):
            misc_settings.center_frame = True
            del(argv[i:i+1])

        elif ((argv[i].lower() == '-raw') or (argv[i].lower() == '-rawout')):
            misc_settings.output_format = 'raw'
            del(argv[i:i+1])

        elif (argv[i].lower() == '-rawin'):
            misc_settings.input_format = 'raw'
            misc_settings.multi = False
            del(argv[i:i+1])

        elif ((argv[i].lower() == '-xyz') or (argv[i].lower() == '-xyzout')):
            misc_settings.output_format = 'xyz'
            del(argv[i:i+1])

        elif (argv[i].lower() == '-xyzin'):
            misc_settings.input_format = 'xyz'
            misc_settings.multi = False
            del(argv[i:i+1])

        elif (argv[i].lower() == '-multi'):
            misc_settings.multi = True
            del(argv[i:i+1])

        elif (argv[i].lower() == '-last'):
            misc_settings.last_frame = True
            misc_settings.multi = False
            del(argv[i:i+1])

        elif (argv[i].lower() == '-interval'):
            misc_settings.skip_interval = int(argv[i+1])
            del(argv[i:i+2])

        elif (argv[i].lower() == '-scale'):
            misc_settings.scale = float(argv[i+1])
            del(argv[i:i+2])

        elif ((argv[i][0] == '-') and (__name__ == "__main__")):
            raise InputError('Error(dump2data): Unrecogized command line argument \"'+argv[i]+'\"\n')
        else:
            i += 1

    usage_examples = \
"""    Typical usage:
dump2data.py orig_file.data < dump.lammpstrj > new_file.data
      (This extracts last frame, uses "full" atom_style.)
    Additional options:
dump2data.py -t t -atomstyle style orig.data < dump.lammpstrj > new.data
"""

    #if __name__ == "__main__":

    if (len(argv) > 2):
        # if there are more than 2 remaining arguments,
        #   AND
        # no other function will process the remaining argument list
        # (ie. if __name__ == "__main__")
        #   THEN
        raise InputError('    ----\n'
                         'ERROR(dump2data): You have too many arguments (or unrecognized arguments):\n'
                         '       \"'+(' '.join(argv))+'\"\n'
                         '    ----\n'
                         +usage_examples)
    elif (len(argv) < 2):
        if misc_settings.output_format == 'data':
            raise InputError('    ----\n'
                             'ERROR(dump2data): Problem with argument list:\n'
                             '       Expected a LAMMPS .data file as an argument.\n'
                             '    ----\n'
                             +usage_examples)
    else:
        in_data_file            = open(argv[1], 'r')
        data_settings.file_name = argv[1];
        data_settings.contents  = in_data_file.readlines()
        in_data_file.close()

    #end of if-then statement for "if __name__ == "__main__""

    if len(data_settings.i_coords) == 0:
        if warning_strings != None:
            warning_strings.append('WARNING(dump2data): atom_style unknown. (Use -atomstyle style. Assuming \"full\")')
        warn_atom_style_unspecified = True
        # The default atom_style is "full"
        data_settings.column_names = AtomStyle2ColNames('full')
        ii_coords = ColNames2Coords(data_settings.column_names)
        # This program assumes that there is only one coordinate triplet
        # (x,y,z) for each atom.  Hence we assume that len(ii_coords)==1
        assert(len(ii_coords) == 1)
        data_settings.i_coords = ii_coords[0]
        data_settings.ii_vects = ColNames2Vects(data_settings.column_names)
        data_settings.i_atomid, data_settings.i_atomtype, data_settings.i_molid = ColNames2AidAtypeMolid(data_settings.column_names)

        ### sys.stderr.write('########################################################\n'
        ###                  '##            WARNING: atom_style unspecified         ##\n'
        ###                  '##    --> \"Atoms\" column data has an unknown format.  ##\n'
        ###                  '##              Assuming atom_style = \"full\"          ##\n'
        ###                  '########################################################\n'
        ###                  '## To specify the \"Atoms\" column format you can:      ##\n'
        ###                  '##   1) Use the -atom_style \"STYLE\"  argument         ##\n'
        ###                  '##      where \"STYLE\" is a string indicating a LAMMPS ##\n'
        ###                  '##      atom_style, including hybrid styles.(Standard ##\n' 
        ###                  '##      atom styles defined in 2011 are supported.)   ##\n'
        ###                  '##   2) Use the -atom_style \"COL_LIST\"    argument    ##\n'
        ###                  '##      where \"COL_LIST" is a quoted list of strings  ##\n'
        ###                  '##      indicating the name of each column.           ##\n'
        ###                  '##      Names \"x\",\"y\",\"z\" are interpreted as          ##\n'
        ###                  '##      atomic coordinates. \"mux\",\"muy\",\"muz\"         ##\n'
        ###                  '##      and \"quati\",\"quatj\",\"quatk\" are               ##\n'
        ###                  '##      interpreted as direction vectors.             ##\n'
        ###                  '##   3) Use the -icoord \"cx cy cz...\" argument        ##\n'
        ###                  '##      where \"cx cy cz\" is a list of integers        ##\n'
        ###                  '##      indicating the column numbers for the x,y,z   ##\n'
        ###                  '##      coordinates of each atom.                     ##\n'
        ###                  '##   4) Use the -ivect \"cmux cmuy cmuz...\" argument   ##\n'
        ###                  '##      where \"cmux cmuy cmuz...\" is a list of        ##\n'
        ###                  '##      integers indicating the column numbers for    ##\n'
        ###                  '##      the vector that determines the direction of a ##\n'
        ###                  '##      dipole or ellipsoid (ie. a rotateable vector).##\n'
        ###                  '##      (More than one triplet can be specified. The  ##\n'
        ###                  '##       number of entries must be divisible by 3.)   ##\n'
        ###                  '##   5) Include a                                     ##\n'
        ###                  '##      write(\"in_init.txt\"){atom_style ...}          ##\n'
        ###                  '##      statement in your .ttree file.                ##\n'
        ###                  '########################################################\n')




def GetIntAtomID(pair):
    return int(pair[0])



def WriteFrameToData(out_file,
                     descr_str,
                     misc_settings,
                     data_settings,
                     natoms,
                     coords,
                     coords_ixiyiz,
                     vects,
                     velocities,
                     xlo_str, xhi_str, 
                     ylo_str, yhi_str, 
                     zlo_str, zhi_str,
                     xy_str, xz_str, yz_str):

    """
    Open a data file.  Read the LAMMPS DATA file line by line.
    When the line contains information which is also in the dump file,
    replace that information with information from the dump file.
    (Information from a dump file is stored in the arguments to this function.)
    The resulting file also has LAMMPS DATA format.

    """

    section = ''
    firstline = True
    for line in data_settings.contents:
        line = line.strip()
        ic = line.find('#')
        if ic != -1:
            line = line[:ic]
        if firstline: # Construct a new descriptive header line:
            if descr_str != None:
                line = descr_str
            firstline = False

        if (len(line) > 0):
            # The initial section (section='') is assumed to be 
            # the "LAMMPS Description" section.  This is where the
            # box boundaries are specified.
            if section == '':
                tokens = line.split()
                if ((len(tokens) >= 2) and
                    ((tokens[-2] == 'xlo') and (tokens[-1] == 'xhi')) and
                    ((xlo_str != None) and (xhi_str != None))):
                    tokens[0] = xlo_str
                    tokens[1] = xhi_str
                    line = ' '.join(tokens)
                elif ((len(tokens) >= 2) and
                      ((tokens[-2] == 'ylo') and (tokens[-1] == 'yhi')) and
                      ((ylo_str != None) and (yhi_str != None))):
                    tokens[0] = ylo_str
                    tokens[1] = yhi_str
                    line = ' '.join(tokens)
                elif ((len(tokens) >= 2) and
                      ((tokens[-2] == 'zlo') and (tokens[-1] == 'zhi')) and
                      ((zlo_str != None) and (zhi_str != None))):
                    tokens[0] = zlo_str
                    tokens[1] = zhi_str
                    line = ' '.join(tokens)
                elif ((len(tokens) >= 3) and
                      ((tokens[-3] == 'xy') and
                       (tokens[-2] == 'xz') and
                       (tokens[-1] == 'yz')) and
                      ((xy_str != None) and
                       (xz_str != None) and
                       (yz_str != None))):
                    tokens[0] = xy_str
                    tokens[1] = xz_str
                    tokens[2] = yz_str
                    line = ' '.join(tokens)
            if (line in set(['Masses', 'Velocities', 'Atoms',
                             'Bond Coeffs', 'Angle Coeffs',
                             'Dihedral Coeffs', 'Improper Coeffs',
                             'Bonds', 'Angles', 'Dihedrals', 'Impropers'])):
                section = line
            else:
                if (section == 'Atoms'):
                    tokens = line.split()
                    atomid = tokens[0]
                    if atomid in coords:
                        # Loop over all of the vector degrees of
                        # freedom of the particle, excluding coords
                        # (for example: mu_x, mu_y, mu_z,
                        #            or quat_i, quat_j, quat_k)
                        # In principle, depending on the atom_style,
                        # there could be multiple vectors per atom.
                        for I in range(0,len(data_settings.ii_vects)):
                            vxvyvz = vects[atomid][I]
                            i_vx = data_settings.ii_vects[I][0]
                            i_vy = data_settings.ii_vects[I][1]
                            i_vz = data_settings.ii_vects[I][2]
                            if ((i_vx >= len(tokens)) or
                                (i_vy >= len(tokens)) or
                                (i_vz >= len(tokens))):
                                raise InputError('Error(dump2data): Atom style incompatible with data file.\n'
                                                 '       Specify the atom_style using -atomstyle style.\n')
                            if ((vxvyvz == None) or
                                (type(vxvyvz) is not tuple)):
                                assert(data_settings.column_names[i_vx] not in dump_column_names)
                                raise InputError('Error(dump2data): You have a vector coordinate in your DATA file named \"'+data_settings.column_names[i_vx]+'\"\n'
                                                 '       However there are no columns with this name in your DUMP file\n'
                                                 '       (or the column was not in the expected place).\n'
                                                 '       Hence, the atom styles in the dump and data files do not match.')

                            # Replace the vector components with numbers
                            # from the dump file
                            tokens[i_vx] = vxvyvz[0]
                            tokens[i_vy] = vxvyvz[1]
                            tokens[i_vz] = vxvyvz[2]

                        # Now loop over the coordinates of each atom.
                        #for I in range(0,len(data_settings.ii_coords)):
                        #    xyz = coords[atomid][I]
                        #            THIS LOOP IS SILLY.
                        #            EACH ATOM ONLY HAS ONE SET OF X,Y,Z 
                        #            COORDINATES. COMMENTING OUT THIS LOOP:
                        #    i_x = data_settings.ii_coords[I][0]
                        #    i_y = data_settings.ii_coords[I][1]
                        #    i_z = data_settings.ii_coords[I][2]
                        # USING THIS INSTEAD:

                        xyz = coords[atomid]
                        i_x = data_settings.i_coords[0]
                        i_y = data_settings.i_coords[1]
                        i_z = data_settings.i_coords[2]
                        if ((i_x >= len(tokens)) or
                            (i_y >= len(tokens)) or
                            (i_z >= len(tokens))):
                            raise InputError('Error(dump2data): Atom style incompatible with data file.\n'
                                             '       Specify the atom_style using -atomstyle style.\n')
                        # Replace the coordinates with coordinates from 
                        # the dump file into tokens[i_x]...
                        tokens[i_x] = str(xyz[0])
                        tokens[i_y] = str(xyz[1])
                        tokens[i_z] = str(xyz[2])

                        # Are there there any integer coords 
                        # (ix, iy, iz) in the dump file?
                        if coords_ixiyiz[atomid]:
                            assert(len(coords_ixiyiz[atomid]) == 3)
                            # Integer coords stored in the DATA file too?
                            if len(tokens)==(len(data_settings.column_names)+3):
                                # Then replace the last 3 columns of the
                                # line in the data file with: ix iy iz
                                tokens[-3] = coords_ixiyiz[atomid][0]
                                tokens[-2] = coords_ixiyiz[atomid][1]
                                tokens[-1] = coords_ixiyiz[atomid][2]
                            else:
                                if (not misc_settings.center_frame):
                                    # Append them to the end of the line:
                                    tokens.append(coords_ixiyiz[atomid][0])
                                    tokens.append(coords_ixiyiz[atomid][1])
                                    tokens.append(coords_ixiyiz[atomid][2])

                        # Now finally paste all the tokens together:
                        line = ' '.join(tokens)


                elif (section == 'Velocities'):
                    tokens = line.split()
                    atomid = tokens[0]
                    if atomid in velocities:

                        vxvyvz = velocities[atomid]
                        if len(tokens) < 4:
                            raise InputError('Error(dump2data): Not enough columns in the \"Velocities\" file.\n')
                        # Replace the coordinates with coordinates from 
                        # the dump file into tokens[i_x]...
                        tokens[1] = str(vxvyvz[0])
                        tokens[2] = str(vxvyvz[1])
                        tokens[3] = str(vxvyvz[2])

                        # Now finally paste all the tokens together:
                        line = ' '.join(tokens)


        out_file.write(line+'\n')















if __name__ == "__main__":

    g_program_name = 'dump2data.py'
    g_date_str     = '2013-10-30'
    g_version_str  = 'v0.43'

    #######  Main Code Below: #######
    sys.stderr.write(g_program_name+' '+g_version_str+' '+g_date_str+' ')
    #if sys.version < '3':
    #    sys.stderr.write(' (python version < 3)\n')
    #else:
    sys.stderr.write('\n')

    try:
        data_settings = DataSettings()
        misc_settings = MiscSettings()
        warning_strings = []
        ParseArgs(sys.argv, 
                  misc_settings, 
                  data_settings, 
                  warning_strings)

        # Open the lammps dump file (trajectory file)
        # Skip to the line containing the correct frame/timestep.
        # (this is the last frame by default).
        # Read the "BOX BOUNDS" and the "ATOMS" sections.
        # Store the x,y,z coordinates in the "coords" associative array
        # (indexed by atom id, which could be non-numeric in general).

        section = ''

        #coords = defaultdict(list)
        #coords_ixiyiz = defaultdict(list)
        #vects = defaultdict(list)
        #xlo_str = xhi_str = ylo_str = yhi_str = zlo_str = zhi_str = None
        #xy_str = xz_str = yz_str = None
        #natoms = -1
        #timestep_str = ''

        frame_coords = defaultdict(list)
        frame_coords_ixiyiz = defaultdict(list)
        frame_vects = defaultdict(list)
        frame_velocities = defaultdict(list)
        frame_xlo_str = frame_xhi_str = None
        frame_ylo_str = frame_yhi_str = None
        frame_zlo_str = frame_zhi_str = None
        frame_xy_str  = frame_xz_str  = frame_yz_str = None
        frame_natoms = -1
        frame_timestep_str = ''
        i_atomid = i_atomtype = i_molid = -1
        i_x  = i_y  = i_z  = i_xu  = i_yu  = i_zu  = -1
        i_xs = i_ys = i_zs = i_xsu = i_ysu = i_zsu = -1

        dump_column_names = []

        #num_frames_in = -1
        num_frames_out = 0
        finished_reading_frame = False
        read_last_frame = False

        #in_coord_file = open('traj_nvt.lammpstrj','r')
        #in_coord_file = open('deleteme.lammpstrj','r')
        in_coord_file = sys.stdin

        while True:

            line = in_coord_file.readline()
            if line == '': # if EOF
                if len(frame_coords) > 0:
                    finished_reading_frame = True
                read_last_frame = True

            line = line.strip()
            if (line.find('ITEM:') == 0):
                section = line
                if (section.find('ITEM: ATOMS ') == 0):
                    dump_column_names = line[12:].split()
                    i_atomid, i_atomtype, i_molid = \
                      ColNames2AidAtypeMolid(dump_column_names)
                    #ii_coords = ColNames2Coords(dump_column_names)

                    if 'x' in dump_column_names:
                        i_x = dump_column_names.index('x')
                    elif 'xu' in dump_column_names:
                        i_xu = dump_column_names.index('xu')
                    elif 'xs' in dump_column_names:
                        i_xs = dump_column_names.index('xs')
                    elif 'xsu' in dump_column_names:
                        i_xsu = dump_column_names.index('xsu')
                    else:
                        raise InputError('Error(dump2data): \"ATOMS\" section of dump file lacks a \"x\" column.\n'+
                                         '       (excerpt below)\n' + line)

                    if 'y' in dump_column_names:
                        i_y = dump_column_names.index('y')
                    elif 'yu' in dump_column_names:
                        i_yu = dump_column_names.index('yu')
                    elif 'ys' in dump_column_names:
                        i_ys = dump_column_names.index('ys')
                    elif 'ysu' in dump_column_names:
                        i_ysu = dump_column_names.index('ysu')
                    else:
                        raise InputError('Error(dump2data): \"ATOMS\" section of dump file lacks a \"y\" column.\n'+
                                         '       (excerpt below)\n' + line)

                    if 'z' in dump_column_names:
                        i_z = dump_column_names.index('z')
                    elif 'zu' in dump_column_names:
                        i_zu = dump_column_names.index('zu')
                    elif 'zs' in dump_column_names:
                        i_zs = dump_column_names.index('zs')
                    elif 'zsu' in dump_column_names:
                        i_zsu = dump_column_names.index('zsu')
                    else:
                        raise InputError('Error(dump2data): \"ATOMS\" section of dump file lacks a \"z\" column.\n'+
                                         '       (excerpt below)\n' + line)





                    ii_vects  = ColNames2Vects(dump_column_names)
                    if (len(ii_vects) != len(data_settings.ii_vects)):
                        raise InputError('Error(dump2data): atom styles in data and dump files differ.\n'
                                         '      Some needed columns from the atom_styles are missing in the dump file.')

                    i_ix = i_iy = i_iz = -1
                    if 'ix' in dump_column_names:
                        i_ix = dump_column_names.index('ix')
                    if 'iy' in dump_column_names:
                        i_iy = dump_column_names.index('iy')
                    if 'iz' in dump_column_names:
                        i_iz = dump_column_names.index('iz')


                    i_vx = i_vy = i_vz = -1
                    if 'vx' in dump_column_names:
                        i_vx = dump_column_names.index('vx')
                    if 'vy' in dump_column_names:
                        i_vy = dump_column_names.index('vy')
                    if 'vz' in dump_column_names:
                        i_vz = dump_column_names.index('vz')


                elif (section.find('ITEM: BOX BOUNDS') == 0):
                    avec=[1.0, 0.0, 0.0]
                    bvec=[0.0, 1.0, 0.0]
                    cvec=[0.0, 0.0, 1.0]

                elif (section.find('ITEM: TIMESTEP') == 0):
                    if len(frame_coords) > 0:
                        finished_reading_frame = True

            elif ((len(line) > 0) and (line[0] != '#')):
                if (section.find('ITEM: TIMESTEP') == 0):
                    finished_reading_frame = False
                    frame_timestep_str  = line
                    frame_coords = defaultdict(list)
                    frame_coords_ixiyiz = defaultdict(list)
                    frame_vects  = defaultdict(list)
                    frame_velocities = defaultdict(list)
                    frame_xlo_str = frame_xhi_str = None
                    frame_ylo_str = frame_yhi_str = None
                    frame_zlo_str = frame_zhi_str = None
                    frame_xy_str  = frame_xz_str  = frame_yz_str = None

                elif (section == 'ITEM: NUMBER OF ATOMS'):
                    frame_natoms = int(line)

                elif (section.find('ITEM: BOX BOUNDS') == 0):
                    is_triclinic = (section.find('xy xz yz') == 0)

                    tokens = line.split()
                    if not frame_xlo_str:
                        assert(not frame_xhi_str)
                        frame_xlo_str = tokens[0]
                        frame_xhi_str = tokens[1]
                        avec[0] = float(frame_xhi_str) - float(frame_xlo_str)
                        if (is_triclinic and (len(tokens) > 2)):
                            frame_xy_str = tokens[2]
                            bvec[0] = float(frame_xy_str)
                            #See http://lammps.sandia.gov/doc/Section_howto.html#howto_12
                        #sys.stderr.write('avec='+str(avec)+'\n')

                    elif not frame_ylo_str:
                        assert(not frame_yhi_str)
                        frame_ylo_str = tokens[0]
                        frame_yhi_str = tokens[1]
                        bvec[1] = float(frame_yhi_str) - float(frame_ylo_str)
                        if (is_triclinic and (len(tokens) > 2)):
                            frame_xz_str = tokens[2]
                            cvec[0] = float(frame_xz_str)
                            #See http://lammps.sandia.gov/doc/Section_howto.html#howto_12
                        #sys.stderr.write('bvec='+str(bvec)+'\n')

                    elif not frame_zlo_str:
                        assert(not frame_zhi_str)
                        frame_zlo_str = tokens[0]
                        frame_zhi_str = tokens[1]
                        cvec = [0.0, 0.0, float(frame_zhi_str) - float(frame_zlo_str)]
                        if (is_triclinic and (len(tokens) > 2)):
                            frame_yz_str = tokens[2]
                            cvec[1] = float(frame_yz_str)
                            #See http://lammps.sandia.gov/doc/Section_howto.html#howto_12
                        #sys.stderr.write('cvec='+str(cvec)+'\n')

                elif (section.find('ITEM: ATOMS') == 0):
                    tokens = line.split()
                    atomid = tokens[i_atomid]

                    if ((i_x != -1) and (i_y != -1) and (i_z != -1)):
                        x  = float(tokens[i_x]) #i_x determined above
                        y  = float(tokens[i_y])
                        z  = float(tokens[i_z])

                    elif ((i_xu != -1) and (i_yu != -1) and (i_zu != -1)):
                        x  = float(tokens[i_xu]) #i_x determined above
                        y  = float(tokens[i_yu])
                        z  = float(tokens[i_zu])

                    elif ((i_xs != -1) and (i_ys != -1) and (i_zs != -1)):
                        xs = float(tokens[i_xs]) #i_xs determined above
                        ys = float(tokens[i_ys])
                        zs = float(tokens[i_zs])

                        x = float(xlo_str) + xs*avec[0] + ys*bvec[0] + zs*cvec[0]
                        y = float(ylo_str) + xs*avec[1] + ys*bvec[1] + zs*cvec[1]
                        z = float(zlo_str) + xs*avec[2] + ys*bvec[2] + zs*cvec[2]

                    # avec, bvec, cvec described here:
                    #http://lammps.sandia.gov/doc/Section_howto.html#howto_12

                    elif ((i_xsu != -1) and (i_ysu != -1) and (i_zsu != -1)):
                        xsu = float(tokens[i_xsu]) #i_xs determined above
                        ysu = float(tokens[i_ysu])
                        zsu = float(tokens[i_zsu])

                        x = float(xlo_str) + xsu*avec[0] + ysu*bvec[0] + zsu*cvec[0]
                        y = float(ylo_str) + xsu*avec[1] + ysu*bvec[1] + zsu*cvec[1]
                        z = float(zlo_str) + xsu*avec[2] + ysu*bvec[2] + zsu*cvec[2]

                    # Now deal with ix, iy, iz
                    if i_ix != -1:
                        ix = int(tokens[i_ix])
                        if (misc_settings.center_frame or
                            (misc_settings.output_format != 'data')):
                            #sys.stderr.write('ix = '+str(ix)+', avec='+str(avec)+'\n')
                            x += ix*avec[0]
                            y += ix*avec[1]
                            z += ix*avec[2]
                        else:
                            if atomid not in frame_coords_ixiyiz:
                                frame_coords_ixiyiz[atomid] = ["0", "0", "0"]
                            else:
                                frame_coords_ixiyiz[atomid][0] = str(ix)

                    if i_iy != -1:
                        iy = int(tokens[i_iy])
                        if (misc_settings.center_frame or
                            (misc_settings.output_format != 'data')):
                            #sys.stderr.write('iy = '+str(iy)+', bvec='+str(bvec)+'\n')
                            x += iy*bvec[0]
                            y += iy*bvec[1]
                            z += iy*bvec[2]
                        else:
                            if atomid not in frame_coords_ixiyiz:
                                frame_coords_ixiyiz[atomid] = ["0", "0", "0"]
                            else:
                                frame_coords_ixiyiz[atomid][1] = str(iy)

                    if i_iz != -1:
                        iz = int(tokens[i_iz])
                        if (misc_settings.center_frame or
                            (misc_settings.output_format != 'data')):
                            #sys.stderr.write('iz = '+str(iz)+', cvec='+str(cvec)+'\n')
                            x += iz*cvec[0]
                            y += iz*cvec[1]
                            z += iz*cvec[2]
                        else:
                            if atomid not in frame_coords_ixiyiz:
                                frame_coords_ixiyiz[atomid] = ["0", "0", "0"]
                            else:
                                frame_coords_ixiyiz[atomid][2] = str(iz)

                    #frame_coords[atomid] = [str(x), str(y), str(z)]
                    frame_coords[atomid] = [x, y, z]

                    vx = 0.0
                    vy = 0.0
                    vz = 0.0
                    if i_vx != -1:
                        vx = float(tokens[i_vx])
                    if i_vy != -1:
                        vy = float(tokens[i_vy])
                    if i_vz != -1:
                        vz = float(tokens[i_vz])

                    frame_velocities[atomid] = [vx, vy, vz]

                    # Ugly detail:
                    # There can be multiple "vects" associated with each atom
                    # (for example, dipole moments, ellipsoid directions, etc..)

                    if atomid not in frame_vects:
                        frame_vects[atomid] = [None for I in range(0,len(ii_vects))]

                    for I in range(0, len(ii_vects)):
                        i_vx   = ii_vects[I][0]
                        i_vy   = ii_vects[I][1]
                        i_vz   = ii_vects[I][2]
                        vx_str = tokens[i_vx]
                        vy_str = tokens[i_vy]
                        vz_str = tokens[i_vz]

                        # Now the annoying part:
                        # Which vect is it (mux,muy,muz) or (quati,quatj,quatk)?
                        # The columns could be listed in a different order
                        # in the data file and in the dump file.
                        # Figure out which vector it is in the data file (stored
                        # in the integer "I_data") so that column names match.
                        name_vx = dump_column_names[i_vx]
                        name_vy = dump_column_names[i_vy]
                        name_vz = dump_column_names[i_vz]
                        i_vx_data = 0
                        I_data = -1
                        # This code is ugly and inneficient.
                        # I never want to touch this code again. (Hope it works)
                        while i_vx_data < len(data_settings.column_names):
                            if name_vx == data_settings.column_names[i_vx_data]:
                                I_data = 0
                                while I_data < len(data_settings.ii_vects):
                                    if ii_vects[I] == data_settings.ii_vects[I_data]:
                                        break
                                    I_data += 1

                            if (0<I_data) and (I_data < len(data_settings.ii_vects)):
                                break

                            i_vx_data += 1

                        if (0 <= I_data) and (I_data < len(data_settings.ii_vects)):
                            frame_vects[atomid][I_data] = (vx_str,vy_str,vz_str)
                        else:
                            raise InputError('Error(dump2data): You have a vector coordinate in your dump file named \"'+name_vx+'\"\n'
                                             '       However there are no columns with this name in your data file\n'
                                             '       (or the column was not in the expected place).\n'
                                             '       Hence, the atom styles in the dump and data files do not match.')
                                             

            if finished_reading_frame:

                if misc_settings.scale != None:
                    for atomid in frame_coords:
                        for d in range(0,3):
                            crd = float(frame_coords[atomid][d])
                            frame_coords[atomid][d] = str(crd*misc_settings.scale)

                if len(frame_coords) != frame_natoms:
                    err_msg = 'Number of lines in \"ITEM: ATOMS\" section disagrees with\n' \
                        + '           \"ITEM: NUMBER OF ATOMS\" declared earlier in this file.\n'
                    raise InputError(err_msg)

                if misc_settings.center_frame:
                    cm = [0.0, 0.0, 0.0]
                    for atomid in frame_coords:
                        for d in range(0,3):
                            cm[d] += float(frame_coords[atomid][d])
                    for d in range(0,3):
                        cm[d] /= float(len(frame_coords))
                    for atomid in frame_coords:
                        for d in range(0,3):
                            frame_coords[atomid][d] = "%.7g" % (float(frame_coords[atomid][d]) - cm[d])
                        frame_coords_ixiyiz[atomid] = ["0","0","0"]

                    if misc_settings.output_format != 'data':
                        frame_coords_ixiyiz[atomid] = ["0","0","0"]



                #if (num_frames_in == -1):
                #    if (misc_settings.timestep_str != ''):
                #        if (float(frame_timestep_str) >= 
                #            float(misc_settings.timestep_str)):
                #            num_frames_in = 1
                #        if not misc_settings.multi:
                #            read_last_frame = True
                #    else:
                #        num_frames_in = 1



                # Should we write out the coordinates in this frame?
                write_this_frame = False

                if misc_settings.multi:

                    write_this_frame = True
                    if (misc_settings.tstart and
                        (int(frame_timestep_str) < misc_settings.tstart)):
                        write_this_frame = False
                    if (misc_settings.tstop and
                        (int(frame_timestep_str) > misc_settings.tstop)):
                        write_this_frame = False
                        read_last_frame = True

                    if misc_settings.tstart:
                        tstart = misc_settings.tstart
                    else:
                        tstart = 0

                    if ((int(frame_timestep_str) - tstart)
                        % 
                        misc_settings.skip_interval) != 0:
                        write_this_frame = False

                else:
                    if misc_settings.last_frame:
                        if read_last_frame:
                            write_this_frame = True
                    else:
                        assert(misc_settings.timestep_str)
                        if (int(frame_timestep_str) >= 
                            int(misc_settings.timestep_str)):
                            write_this_frame = True
                            read_last_frame  = True


                if write_this_frame:

                    num_frames_out += 1

                    sys.stderr.write('  (writing frame '+str(num_frames_out)+
                                     ' at timestep '+frame_timestep_str+')\n')


                    # Print the frame
                    # First check which format to output the data:
                    if misc_settings.output_format == 'raw':
                        # Print out the coordinates in simple 3-column text format
                        for atomid, xyz in iter(sorted(frame_coords.items(), key=GetIntAtomID)):
                            if misc_settings.scale == None:
                                sys.stdout.write(str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2])+'\n')
                            else:
                                # Only convert to float and back if misc_settings.scale != None
                                sys.stdout.write(str(misc_settings.scale*float(xyz[0]))+' '+
                                                 str(misc_settings.scale*float(xyz[1]))+' '+
                                                 str(misc_settings.scale*float(xyz[2]))+'\n')
                        sys.stdout.write('\n')

                    elif misc_settings.output_format == 'xyz':
                            # Print out the coordinates in simple 3-column text format
                        sys.stdout.write(str(len(frame_coords))+'\n')
                        descr_str = 'LAMMPS data from timestep '+frame_timestep_str
                        sys.stdout.write(descr_str+'\n')
                        for atomid, xyz in iter(sorted(frame_coords.items(), key=GetIntAtomID)):
                            if misc_settings.scale == None:
                                sys.stdout.write(str(atomid)+' '+
                                                 str(xyz[0])+' '+
                                                 str(xyz[1])+' '+
                                                 str(xyz[2])+'\n')
                            else:
                                # Only convert to float and back if misc_settings.scale != None
                                sys.stdout.write(str(atomid)+' '+
                                                 str(misc_settings.scale*float(xyz[0]))+' '+
                                                 str(misc_settings.scale*float(xyz[1]))+' '+
                                                 str(misc_settings.scale*float(xyz[2]))+'\n')

                    else:
                        # Parse the DATA file specified by the user
                        # and replace appropriate lines or fields with
                        # the corresponding text from the DUMP file.
                        descr_str = 'LAMMPS data from timestep '+frame_timestep_str
                        if misc_settings.multi and (misc_settings.output_format == 'data'):
                            out_file_name = data_settings.file_name + '.'\
                                + str(num_frames_out)
                            sys.stderr.write('  (creating file \"'+out_file_name+'\")\n')
                            out_file = open(out_file_name, 'w')
                        else:
                            out_file = sys.stdout

                        WriteFrameToData(out_file,
                                         descr_str,
                                         misc_settings,
                                         data_settings,
                                         frame_natoms,
                                         frame_coords,
                                         frame_coords_ixiyiz,
                                         frame_vects,
                                         frame_velocities,
                                         frame_xlo_str, frame_xhi_str, 
                                         frame_ylo_str, frame_yhi_str, 
                                         frame_zlo_str, frame_zhi_str,
                                         frame_xy_str, frame_xz_str, frame_yz_str)

                        #if misc_settings.multi:
                        #    out_file.close()


                #if num_frames_in >= 0:
                #    num_frames_in += 1


                if read_last_frame:
                    exit(0)


        for warning_str in warning_strings:
            sys.stderr.write(warning_str+'\n')



    except (ValueError, InputError) as err:
        sys.stderr.write('\n'+str(err)+'\n')
        sys.exit(-1)

