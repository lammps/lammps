#!/usr/bin/env python

"""
   Generate a moltemplate (.lt) file containing a definition of a polymer
   molecule whose monomers are located at the positions specified in
   "coords.raw" (a 3-column text file).  Monomers will be rotated so
   that they point in the direction connecting neighbors (r[i+1]-r[i])
   The user can specify the subunits to use when building the polymer,
   the atoms to to build bonds (and angles, and dihedrals) between monomers
   and the helical pitch of the polymer.  The output of this program is
   a text file in moltemplate (.lt) format containing the sequence of
   moltemplate commands needed to build this polymer molecule(s).  (One must
   then run moltemplate on this file to build the LAMMPS simulation files.)
       Multiple Polymers:
   To make it easier to create polymer melts, multiple polymers can be created
   from coordinates in the same file by using the "-cuts" command line argument.
      Encapsulation:
   If the "-polymer-name PolyName" command line option is given, then these
   moltemplate commands will be nested within the definition of a moltemplate
   object (named "PolyName", in this example. Later in your moltemplate files,
   you must remember to instantiate a copy of this moltemplate object using
   a command like "polymer = new PolyName"  Atoms within this object will
   share the same molecule-ID number.)  If multiple polymers are requested, then
   each of them will have their own polymer object.

"""


g_usage_msg = """
Usage:

   genpoly_lt.py  \\
      [-bond btype a1 a2] \\
      [-helix deltaphi] \\
      [-axis x,y,z] \\
      [-circular yes/no/connected] \\
      [-dir-indices ia ib] \\
      [-angle    atype a1 a2 a3 i1 i2 i3] \\
      [-dihedral dtype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-improper itype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-monomer-name mname] \\
      [-sequence sequence.txt] \\
      [-polymer-name pname] \\
      [-inherits ForceFieldObject] \\
      [-header "import \"monomer.lt\""] \\
      [-cuts cuts.txt] \\
      [-box paddingX,paddingY,paddingZ] \\
      < coords.raw > polymer.lt

"""


import sys
import random
from math import *


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return self.err_msg

    def __repr__(self):
        return str(self)


class GPSettings(object):

    def __init__(self):
        self.direction_orig = [1.0, 0.0, 0.0]
        self.is_circular = False
        self.connect_ends = False
        self.delta_phi = 0.0
        self.header = 'import \"forcefield.lt\"'
        self.name_monomer = 'Monomer'
        self.name_polymer = ''
        self.inherits = ''
        self.name_sequence = []
        self.dir_index_offsets = (-1,1)
        self.cuts = []
        self.box_padding = None
        self.bonds_name = []
        self.bonds_type = []
        self.bonds_atoms = []
        self.bonds_index_offsets = []
        self.angles_name = []
        self.angles_type = []
        self.angles_atoms = []
        self.angles_index_offsets = []
        self.dihedrals_name = []
        self.dihedrals_type = []
        self.dihedrals_atoms = []
        self.dihedrals_index_offsets = []
        self.impropers_name = []
        self.impropers_type = []
        self.impropers_atoms = []
        self.impropers_index_offsets = []

    def ParseArgs(self, argv):
        i = 1
        while i < len(argv):
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
            if argv[i].lower() == '-bond':
                if i + 3 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 4 strings.\n')
                # self.bonds_name.append(argv[i+1])
                self.bonds_type.append(argv[i + 1])
                self.bonds_atoms.append((argv[i + 2],
                                         argv[i + 3]))
                self.bonds_index_offsets.append((0, 1))
                del(argv[i:i + 4])
            elif argv[i].lower() == '-angle':
                if i + 7 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 5 strings and 3 integers.\n')
                # self.angles_name.append(argv[i+1])
                self.angles_type.append(argv[i + 1])
                self.angles_atoms.append((argv[i + 2],
                                          argv[i + 3],
                                          argv[i + 4]))
                self.angles_index_offsets.append((int(argv[i + 5]),
                                                  int(argv[i + 6]),
                                                  int(argv[i + 7])))
                if ((self.angles_index_offsets[-1][0] < 0) or
                        (self.angles_index_offsets[-1][1] < 0) or
                        (self.angles_index_offsets[-1][2] < 0)):
                    raise InputError(
                        'Error: ' + argv[i] + ' indices (i1 i2 i3) must be >= 0\n')
                del(argv[i:i + 8])
            elif argv[i].lower() == '-dihedral':
                if i + 9 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 6 strings and 4 integers.\n')
                # self.dihedrals_name.append(argv[i+1])
                self.dihedrals_type.append(argv[i + 1])
                self.dihedrals_atoms.append((argv[i + 2],
                                             argv[i + 3],
                                             argv[i + 4],
                                             argv[i + 5]))
                self.dihedrals_index_offsets.append((int(argv[i + 6]),
                                                     int(argv[i + 7]),
                                                     int(argv[i + 8]),
                                                     int(argv[i + 9])))
                if ((self.dihedrals_index_offsets[-1][0] < 0) or
                        (self.dihedrals_index_offsets[-1][1] < 0) or
                        (self.dihedrals_index_offsets[-1][2] < 0) or
                        (self.dihedrals_index_offsets[-1][3] < 0)):
                    raise InputError(
                        'Error: ' + argv[i] + ' indices (i1 i2 i3 i4) must be >= 0\n')
                del(argv[i:i + 10])
            elif argv[i].lower() == '-improper':
                if i + 9 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by 6 strings and 4 integers.\n')
                # self.impropers_name.append(argv[i+1])
                self.impropers_type.append(argv[i + 1])
                self.impropers_atoms.append((argv[i + 2],
                                             argv[i + 3],
                                             argv[i + 4],
                                             argv[i + 5]))
                self.impropers_index_offsets.append((int(argv[i + 6]),
                                                     int(argv[i + 7]),
                                                     int(argv[i + 8]),
                                                     int(argv[i + 9])))
                if ((self.impropers_index_offsets[-1][0] < 0) or
                        (self.impropers_index_offsets[-1][1] < 0) or
                        (self.impropers_index_offsets[-1][2] < 0) or
                        (self.impropers_index_offsets[-1][3] < 0)):
                    raise InputError(
                        'Error: ' + argv[i] + ' indices (i1 i2 i3 i4) must be >= 0\n')
                del(argv[i:i + 10])
            elif (argv[i].lower() == '-monomer-name'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string\n')
                self.name_monomer = argv[i + 1]
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-sequence'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')
                try:
                    f = open(argv[i + 1], "r")
                except IOError:
                    raise InputError(
                        'Error: file ' + argv[i + 1] + ' could not be opened for reading\n')
                self.name_sequence = []
                for line_orig in f:
                    line = line_orig.strip()
                    ic = line.find('#')
                    if ic != -1:
                        line = line[:ic]
                    else:
                        line = line.strip()
                    if len(line) > 0:
                        self.name_sequence.append(line)
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-cuts'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a file name\n')
                try:
                    f = open(argv[i + 1], "r")
                except IOError:
                    raise InputError(
                        'Error: file ' + argv[i + 1] + ' could not be opened for reading\n')
                self.name_sequence = []
                for line_orig in f:
                    line = line_orig.strip()
                    ic = line.find('#')
                    if ic != -1:
                        line = line[:ic]
                    else:
                        line = line.strip()
                    if len(line) > 0:
                        try:
                            self.cuts.append(int(line))
                        except ValueError:
                            raise InputError(
                                'Error: file ' + argv[i + 1] + ' should contain only nonnegative integers.\n')
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-polymer-name'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string\n')
                self.name_polymer = argv[i + 1]
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-inherits'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string\n')
                self.inherits = argv[i + 1]
                if self.inherits.find('inherits ') == 0:
                    self.inherits = ' ' + self.inherits
                else:
                    self.inherits = ' inherits ' + self.inherits
                if self.name_polymer == '':
                    self.name_polymer = 'Polymer'  # supply a default name
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-header'):
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a string (usually in quotes)\n')
                self.header = argv[i + 1]
                del(argv[i:i + 2])
            elif argv[i].lower() == '-axis':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed ' +
                                     'by 3 numbers separated by commas (no spaces)\n')
                self.direction_orig = list(map(float, argv[i + 1].split(',')))
                del(argv[i:i + 2])
            elif argv[i].lower() == '-circular':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by an argument\n' +
                                     '       ("yes", "no", or "connected")\n')
                if argv[i + 1].lower() == 'yes':
                    self.connect_ends = True
                    self.is_circular = True
                elif argv[i + 1].lower() == 'connected':
                    self.connect_ends = True
                    self.is_circular = False
                elif argv[i + 1].lower() == 'no':
                    self.connect_ends = False
                    self.is_circular = False
                else:
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by an argument\n' +
                                     '       ("yes", "no", or "connected")\n')
                del(argv[i:i + 2])
            elif argv[i].lower() == '-helix':
                if i + 1 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by a number (angle in degrees)\n')
                self.delta_phi = float(argv[i + 1])
                del(argv[i:i + 2])
            elif (argv[i].lower() == '-dir-indices'):
                if i + 2 >= len(argv):
                    raise InputError(
                        'Error: ' + argv[i] + ' flag should be followed by two integers\n')
                self.dir_index_offsets = (int(argv[i + 1]), int(argv[i + 2]))
                if self.dir_index_offsets[0] == self.dir_index_offsets[1]:
                    raise InputError(
                        'Error: The two numbers following ' + argv[i] + ' must not be equal.\n')
                del(argv[i:i + 3])
            elif (argv[i].lower() == '-box'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed ' +
                                     'by 3 numbers separated by commas (no spaces)\n')
                self.box_padding = list(map(float, argv[i + 1].split(',')))
                if len(self.box_padding) == 1:
                    self.box_padding = self.box_padding * 3
                del(argv[i:i + 2])

            # elif ((argv[i][0] == '-') and (__name__ == '__main__')):
            #
            #    raise InputError('Error('+g_program_name+'):\n'+\
            #        'Unrecogized command line argument \"'+argv[i]+\
            #        '\"\n\n'+\
            #        __doc__)
            else:
                i += 1

        for b in range(0, len(self.bonds_type)):
            if len(self.bonds_type) > 1:
                self.bonds_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.bonds_name.append('genpoly')
        for b in range(0, len(self.angles_type)):
            if len(self.angles_type) > 1:
                self.angles_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.angles_name.append('genpoly')
        for b in range(0, len(self.dihedrals_type)):
            if len(self.dihedrals_type) > 1:
                self.dihedrals_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.dihedrals_name.append('genpoly')
        for b in range(0, len(self.impropers_type)):
            if len(self.impropers_type) > 1:
                self.impropers_name.append('genpoly' + str(b + 1) + '_')
            else:
                self.impropers_name.append('genpoly')


class WrapPeriodic(object):
    """ Wrap() calculates the remainder of i % N.
        It turns out to be convenient to do this multiple times and later
        query whether i/N != 0 in any of them once (by checking bounds_err).

    """
    bounds_err = False

    @classmethod
    def Wrap(obj, i, N):
        if i // N != 0:
            obj.bounds_err = True
        return i % N

    def WrapF(obj, x, L):
        i = floor(x / L)
        if i != 0:
            obj.bounds_err = True
        return x - i * L


class GenPoly(object):
    """
        Read coordinates from a file, and generate a list of \"new\" commands
        in moltemplate format with the position of each monomer located
        at these positions, oriented appropriately, with bonds (and angles,
        dihedrals, etc...) connecting successive monomers together.
        By default (if settings.cuts==False) only a single polymer is created.
        However this class can create multiple polymers of different lengths.
        The list of coordinates for each polymer are saved separately within
        the "self.coords_multi" member.

    """

    def __init__(self):
        self.settings = GPSettings()
        self.coords_multi = []  # a list-of-list-of-lists of numbers Nxnx3
        self.direction_vects = []
        self.box_bounds_min = [0.0, 0.0, 0.0]
        self.box_bounds_max = [0.0, 0.0, 0.0]
        self.N = 0

    def ParseArgs(self, argv):
        # The command above will remove arguments from argv which are
        # understood by GPSettings.ParseArgs(argv).  
        # The remaining arguments will be handled below.
        self.settings.ParseArgs(argv)

    def ReadCoords(self, infile):
        coords = []
        lines = infile.readlines()
        for i in range(0, len(lines)):
            tokens = lines[i].strip().split()
            if (len(tokens) == 3):
                coords.append(list(map(float, tokens)))

        self.N = len(coords)
        if self.N < 2:
            raise InputError(
                "Error: Coordinate file must have at least 2 positions.\n")
        # Now generate self.settings.name_sequence:
        if len(self.settings.name_sequence) != self.N:
            self.settings.name_sequence = [self.settings.name_monomer] * self.N

        self.settings.cuts.append(self.N + 1)
        self.settings.cuts.sort()
        i = 0
        for j in self.settings.cuts:
            self.coords_multi.append(coords[i:j])
            i = j

    def ChooseDirections(self, coords):
        """
        Calculate the direction each monomer subunit should be pointing at:

        """

        self.N = len(coords)
        self.direction_vects = [[0.0, 0.0, 0.0] for i in range(0, self.N + 1)]

        if self.settings.is_circular:
            for i in range(0, self.N):
                # By default, the direction that monomer "i" is pointing is
                # determined by the position of the monomers before and after it
                # (at index i-1, and i+1).  More generally, we allow the user
                # to choose what these offsets are ("dir_index_offsets[")
                ia = WrapPeriodic.Wrap(i + self.settings.dir_index_offsets[0],
                                       self.N)
                ib = WrapPeriodic.Wrap(i + self.settings.dir_index_offsets[1],
                                       self.N)
                for d in range(0, 3):
                    self.direction_vects[i][d] = coords[
                        ib][d] - coords[ia][d]
        else:
            for i in range(1, self.N - 1):
                for d in range(0, 3):
                    self.direction_vects[i][d] = coords[
                        i + self.settings.dir_index_offsets[1]][d] - coords[
                            i + self.settings.dir_index_offsets[0]][d]

            for d in range(0, 3):
                self.direction_vects[0][d] = coords[1][d] - coords[0][d]
                self.direction_vects[
                    self.N - 1][d] = coords[self.N - 1][d] - coords[self.N - 2][d]

        # Optional: normalize the direction vectors

        for i in range(0, self.N):
            direction_len = 0.0
            for d in range(0, 3):
                direction_len += (self.direction_vects[i][d])**2
            direction_len = sqrt(direction_len)
            for d in range(0, 3):
                self.direction_vects[i][d] /= direction_len

        # Special case:  self.direction_vects[-1] is the direction that the original monomer
        # in "monomer.lt" was pointing.  (By default, 1,0,0 <--> the "x"
        # direction)

        self.direction_vects[-1] = self.settings.direction_orig

    def WriteLTFile(self, outfile):
        """ Write an moltemplate (.lt) file containing the definition of
        this polymer object.  (If multiple polymer objects were requested by
        the user (using the -cuts argument), then their definitions will
        appear nested within this object, and each of them will be
        instantiated once when the parent object is instantiated.)

        """

        outfile.write(self.settings.header + "\n\n\n")
        if len(self.coords_multi) == 1:
            self.WritePolymer(outfile,
                              self.settings.name_polymer +
                              self.settings.inherits,
                              self.coords_multi[0])
        else:
            if self.settings.name_polymer != '':
                outfile.write(self.settings.name_polymer + " {\n\n")
            outfile.write('# Definitions of individual polymers to follow\n\n')
            for i in range(0, len(self.coords_multi)):
                self.WritePolymer(outfile,
                                  self.settings.name_polymer + '_sub' + str(i + 1) +
                                  self.settings.inherits,
                                  self.coords_multi[i])
            outfile.write('\n\n'
                          '# Now instantiate all the polymers (once each)\n\n')

            for i in range(0, len(self.coords_multi)):
                outfile.write('polymers[' + str(i) + '] = new ' +
                              self.settings.name_polymer + '_sub' + str(i + 1) + '\n')

            if self.settings.name_polymer != '':
                outfile.write('\n\n'
                              '}  # ' + self.settings.name_polymer + '\n\n')

        if self.settings.box_padding != None:
            for i in range(0, len(self.coords_multi)):
                # calculate the box big enough to collectively enclose
                # all of the coordinates (even multiple coordinate sets)
                self.CalcBoxBoundaries(self.coords_multi[i])
            self.WriteBoxBoundaries(outfile)

    def WritePolymer(self,
                     outfile,
                     name_polymer,
                     coords):
        """ Write a single polymer object to a file.
            This function is invoked by WriteLTFile()

        """
        self.ChooseDirections(coords)

        if name_polymer != '':
            outfile.write(name_polymer + ' {\n'
                          '\n\n\n'
                          'create_var {$mol}\n'
                          '# The line above forces all monomer subunits to share the same molecule-ID\n'
                          '# (Note: Setting the molecule-ID number is optional and is usually ignored.)\n\n\n\n')

        outfile.write("""
# ------------ List of Monomers: ------------
#
# (Note: move(), rot(), and rotvv() commands control the position
#  of each monomer.  (See the moltemplate manual for an explanation
#  of what they do.)  Commands enclosed in push() are cumulative
#  and remain in effect until removed by pop().)



"""
                      )

        outfile.write("push(move(0,0,0))\n")

        for i in range(0, self.N):
            #im1 = i-1
            # if im1 < 0 or self.settings.connect_ends:
            #    if im1 < 0:
            #        im1 += self.N
            outfile.write("pop()\n")
            outfile.write("push(rotvv(" +
                          str(self.direction_vects[i - 1][0]) + "," +
                          str(self.direction_vects[i - 1][1]) + "," +
                          str(self.direction_vects[i - 1][2]) + "," +
                          str(self.direction_vects[i][0]) + "," +
                          str(self.direction_vects[i][1]) + "," +
                          str(self.direction_vects[i][2]) + "))\n")
            # Recall that self.direction_vects[-1] =
            # self.settings.direction_orig  (usually 1,0,0)
            outfile.write("push(move(" +
                          str(coords[i][0]) + "," +
                          str(coords[i][1]) + "," +
                          str(coords[i][2]) + "))\n")

            outfile.write("mon[" + str(i) + "] = new " +
                          self.settings.name_sequence[i] +
                          ".rot(" + str(self.settings.delta_phi * i) + ",1,0,0)\n")

        assert(len(self.settings.bonds_name) ==
               len(self.settings.bonds_type) ==
               len(self.settings.bonds_atoms) ==
               len(self.settings.bonds_index_offsets))
        if len(self.settings.bonds_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Bonds\") {\n")
        WrapPeriodic.bounds_err = False
        for i in range(0, self.N):
            test = False
            for b in range(0, len(self.settings.bonds_type)):
                I = i + self.settings.bonds_index_offsets[b][0]
                J = i + self.settings.bonds_index_offsets[b][1]
                I = WrapPeriodic.Wrap(I, self.N)
                J = WrapPeriodic.Wrap(J, self.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write(
                    "  $bond:" + self.settings.bonds_name[b] + str(i + 1))
                if len(self.settings.bonds_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @bond:" + self.settings.bonds_type[b] + " $atom:mon[" + str(I) + "]/" + self.settings.bonds_atoms[
                              b][0] + " $atom:mon[" + str(J) + "]/" + self.settings.bonds_atoms[b][1] + "\n")
        if len(self.settings.bonds_type) > 0:
            outfile.write("}  # write(\"Data Bonds\") {...\n\n\n")

        assert(len(self.settings.angles_name) ==
               len(self.settings.angles_type) ==
               len(self.settings.angles_atoms) ==
               len(self.settings.angles_index_offsets))
        if len(self.settings.angles_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Angles\") {\n")
        for i in range(0, self.N):
            for b in range(0, len(self.settings.angles_type)):
                I = i + self.settings.angles_index_offsets[b][0]
                J = i + self.settings.angles_index_offsets[b][1]
                K = i + self.settings.angles_index_offsets[b][2]
                I = WrapPeriodic.Wrap(I, self.N)
                J = WrapPeriodic.Wrap(J, self.N)
                K = WrapPeriodic.Wrap(K, self.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write(
                    "  $angle:" + self.settings.angles_name[b] + str(i + 1))
                if len(self.settings.angles_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @angle:" + self.settings.angles_type[b] +
                              " $atom:mon[" + str(I) + "]/" + self.settings.angles_atoms[b][0] +
                              " $atom:mon[" + str(J) + "]/" + self.settings.angles_atoms[b][1] +
                              " $atom:mon[" + str(K) + "]/" + self.settings.angles_atoms[b][2] +
                              "\n")
        if len(self.settings.angles_type) > 0:
            outfile.write("}  # write(\"Data Angles\") {...\n\n\n")

        assert(len(self.settings.dihedrals_name) ==
               len(self.settings.dihedrals_type) ==
               len(self.settings.dihedrals_atoms) ==
               len(self.settings.dihedrals_index_offsets))
        if len(self.settings.dihedrals_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Dihedrals\") {\n")
        for i in range(0, self.N):
            for b in range(0, len(self.settings.dihedrals_type)):
                I = i + self.settings.dihedrals_index_offsets[b][0]
                J = i + self.settings.dihedrals_index_offsets[b][1]
                K = i + self.settings.dihedrals_index_offsets[b][2]
                L = i + self.settings.dihedrals_index_offsets[b][3]
                I = WrapPeriodic.Wrap(I, self.N)
                J = WrapPeriodic.Wrap(J, self.N)
                K = WrapPeriodic.Wrap(K, self.N)
                L = WrapPeriodic.Wrap(L, self.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write("  $dihedral:" +
                              self.settings.dihedrals_name[b] + str(i + 1))
                if len(self.settings.dihedrals_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @dihedral:" + self.settings.dihedrals_type[b] +
                              " $atom:mon[" + str(I) + "]/" + self.settings.dihedrals_atoms[b][0] +
                              " $atom:mon[" + str(J) + "]/" + self.settings.dihedrals_atoms[b][1] +
                              " $atom:mon[" + str(K) + "]/" + self.settings.dihedrals_atoms[b][2] +
                              " $atom:mon[" + str(L) + "]/" + self.settings.dihedrals_atoms[b][3] +
                              "\n")
        if len(self.settings.dihedrals_type) > 0:
            outfile.write("}  # write(\"Data Dihedrals\") {...\n\n\n")

        assert(len(self.settings.impropers_name) ==
               len(self.settings.impropers_type) ==
               len(self.settings.impropers_atoms) ==
               len(self.settings.impropers_index_offsets))
        if len(self.settings.impropers_type) > 0:
            outfile.write("\n"
                          "\n"
                          "write(\"Data Impropers\") {\n")
        for i in range(0, self.N):
            for b in range(0, len(self.settings.impropers_type)):
                I = i + self.settings.impropers_index_offsets[b][0]
                J = i + self.settings.impropers_index_offsets[b][1]
                K = i + self.settings.impropers_index_offsets[b][2]
                L = i + self.settings.impropers_index_offsets[b][3]
                I = WrapPeriodic.Wrap(I, self.N)
                J = WrapPeriodic.Wrap(J, self.N)
                K = WrapPeriodic.Wrap(K, self.N)
                L = WrapPeriodic.Wrap(L, self.N)
                if WrapPeriodic.bounds_err:
                    WrapPeriodic.bounds_err = False
                    if not self.settings.connect_ends:
                        continue
                outfile.write("  $improper:" +
                              self.settings.impropers_name[b] + str(i + 1))
                if len(self.settings.impropers_type) > 1:
                    outfile.write("_" + str(b + 1))
                outfile.write(" @improper:" + self.settings.impropers_type[b] +
                              " $atom:mon[" + str(I) + "]/" + self.settings.impropers_atoms[b][0] +
                              " $atom:mon[" + str(J) + "]/" + self.settings.impropers_atoms[b][1] +
                              " $atom:mon[" + str(K) + "]/" + self.settings.impropers_atoms[b][2] +
                              " $atom:mon[" + str(L) + "]/" + self.settings.impropers_atoms[b][3] +
                              "\n")
        if len(self.settings.impropers_type) > 0:
            outfile.write("}  # write(\"Data Impropers\") {...\n\n\n")

        if name_polymer != '':
            outfile.write("}  # " + name_polymer + "\n\n\n\n")

    def CalcBoxBoundaries(self, coords):
        N = len(coords)
        for i in range(0, N):
            for d in range(0, 3):
                if not self.box_bounds_min:
                    assert(not self.box_bounds_max)
                    self.box_bounds_min = [xd for xd in coords[i]]
                    self.box_bounds_max = [xd for xd in coords[i]]
                else:
                    if coords[i][d] > self.box_bounds_max[d]:
                        self.box_bounds_max[d] = coords[i][d]
                    if coords[i][d] < self.box_bounds_min[d]:
                        self.box_bounds_min[d] = coords[i][d]

    def WriteBoxBoundaries(self, outfile):
        for d in range(0, 3):
            self.box_bounds_min[d] -= self.settings.box_padding[d]
            self.box_bounds_max[d] += self.settings.box_padding[d]
        outfile.write("\n# ---------------- simulation box -----------------\n"
            
                      "# Now define a box big enough to hold a polymer with this (initial) shape\n"
                      "\n\n"
                      "write_once(\"Data Boundary\") {\n"
                      + str(self.box_bounds_min[0]) + "  " +
                      str(self.box_bounds_max[0]) + " xlo xhi\n"
                      + str(self.box_bounds_min[1]) + "  " +
                      str(self.box_bounds_max[1]) + " ylo yhi\n"
                      + str(self.box_bounds_min[2]) + "  " +
                      str(self.box_bounds_max[2]) + " zlo zhi\n"
                      "}\n\n\n")


def main():
    try:
        g_program_name = __file__.split('/')[-1]
        g_version_str = '0.0.5'
        g_date_str = '2017-4-14'
        sys.stderr.write(g_program_name + ' v' +
                         g_version_str + ' ' + g_date_str + '\n')
        argv = [arg for arg in sys.argv]
        infile = sys.stdin
        outfile = sys.stdout
        genpoly = GenPoly()
        genpoly.ParseArgs(argv)
        # Any remain arguments?
        if len(argv) > 1:
            raise InputError('Error(' + g_program_name + '):\n' +
                             'Unrecogized command line argument \"' + argv[1] +
                             '\"\n\n' +
                             g_usage_msg)
        genpoly.ReadCoords(infile)
        genpoly.WriteLTFile(outfile)

    except (ValueError, InputError) as err:
        sys.stderr.write('\n' + str(err) + '\n')
        sys.exit(-1)

    return

if __name__ == '__main__':
    main()
