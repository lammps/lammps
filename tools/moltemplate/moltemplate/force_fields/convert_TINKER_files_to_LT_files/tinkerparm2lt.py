#! /usr/bin/env python

"""
This standalone python script can be used to convert the force-fields 
distributed with TINKER (".prm" files) into moltemplate (".lt") format.

This script has only been tested with the OPLSAA force field (oplsaa.prm).
The full list of force-field files distributed with TINKER is available at:
http://dasher.wustl.edu/tinker/distribution/params

Other TINKER compatible force fields may not work.  (One reason for this is that
some force fields require features, such as point-dipole polarizability, which,
to my knowledge, have not yet been implemented in LAMMPS as of 2017-2-01.)

"""


__author__ = 'Jason Lambert and Andrew Jewett'
# (some additional corrections by Miguel Gonzalez, Yue Chun Chiu and others)
__version__ = '0.3.2'
__date__ = '2018-6-15'


import sys
import os
from sets import Set
from operator import itemgetter


g_program_name = __file__.split('/')[-1]


doc_msg = \
    "Typical Usage:\n\n" + \
    "   " + g_program_name + " -name OPLS < oplsaa.prm > oplsaa.lt\n\n" + \
    "   where \"oplsaa.prm\" is a force-field file downloaded from the TINKER website,\n" + \
    "         \"oplsaa.lt\" is the corresponding file converted into moltemplate (.lt) format.\n" + \
    "   and   \"OPLS\" is the name that future moltemplate users will use to refer to\n" + \
    "         this force-field (optional).\n" + \
    "Optional Arguments\n" + \
    "   -name FORCEFIELDNAME # Give the force-field a name (recommended)\n" + \
    "   -file FILE_NAME      # Read force field parameters from a file\n" + \
    "   -url URL             # Read force field parameters from a file online\n" + \
    "   -atoms \"QUOTED LIST\" # Restrict output to a subset of atom types\n" + \
    "   -hybrid              # Optional LAMMPS \"hybrid\" style compatibility\n" + \
    "   -zeropad N          # Optional zero-padding for bonded interactions\n"


def SplitQuotedString(string,
                      quotes='\'\"',
                      delimiters=' \t\r\f\n',
                      escape='\\',
                      comment_char='#'):
    tokens = []
    token = ''
    reading_token = True
    escaped_state = False
    quote_state = None
    for c in string:

        if (c in comment_char) and (not escaped_state) and (quote_state == None):
            tokens.append(token)
            return tokens

        elif (c in delimiters) and (not escaped_state) and (quote_state == None):
            if reading_token:
                tokens.append(token)
                token = ''
                reading_token = False

        elif c in escape:
            if escaped_state:
                token += c
                reading_token = True
                escaped_state = False
            else:
                escaped_state = True
                # and leave c (the '\' character) out of token
        elif (c in quotes) and (not escaped_state):
            if (quote_state != None):
                if (c == quote_state):
                    quote_state = None
            else:
                quote_state = c
            token += c
            reading_token = True
        else:
            if (c == 'n') and (escaped_state == True):
                c = '\n'
            elif (c == 't') and (escaped_state == True):
                c = '\t'
            elif (c == 'r') and (escaped_state == True):
                c = '\r'
            elif (c == 'f') and (escaped_state == True):
                c = '\f'
            token += c
            reading_token = True
            escaped_state = False

    if len(string) > 0:
        tokens.append(token)
    return tokens


def RemoveOuterQuotes(text, quotes='\"\''):
    if ((len(text) >= 2) and (text[0] in quotes) and (text[-1] == text[0])):
        return text[1:-1]
    else:
        return text



def main():
    try:
        sys.stderr.write(g_program_name + ", version " +
                         __version__ + ", " + __date__ + "\n")
        if sys.version < '2.6':
            raise Exception('Error: Using python ' + sys.version + '\n' +
                            '       Alas, your version of python is too old.\n'
                            '       You must upgrade to a newer version of python (2.6 or later).')
    
        if sys.version < '2.7':
            from ordereddict import OrderedDict
        else:
            from collections import OrderedDict
    
        if sys.version > '3':
            import io
        else:
            import cStringIO
    
        # defaults:
        ffname = "TINKER_FORCE_FIELD"
        type_subset = Set([])
        filename_in = ""
        file_in = sys.stdin
        pair_style_name = "lj/cut/coul/long"
        pair_style_link = "http://lammps.sandia.gov/doc/pair_lj.html"
        bond_style_name = "harmonic"
        bond_style_link = "http://lammps.sandia.gov/doc/bond_harmonic.html"
        angle_style_name = "harmonic"
        angle_style_link = "http://lammps.sandia.gov/doc/angle_harmonic.html"
        dihedral_style_name = "fourier"
        dihedral_style_link = "http://lammps.sandia.gov/doc/dihedral_fourier.html"
        improper_style_name = "harmonic"
        improper_style_link = "http://lammps.sandia.gov/doc/improper_harmonic.html"
        #improper_style_name = "cvff"
        #improper_style_link = "http://lammps.sandia.gov/doc/improper_cvff.html"
        special_bonds_command = "special_bonds lj/coul 0.0 0.0 0.5"
        mixing_style = "geometric"
        use_hybrid = False
        contains_united_atoms = False
        zeropad_ffid = 1
    
        argv = [arg for arg in sys.argv]
    
        i = 1
    
        while i < len(argv):
    
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
    
            if argv[i] == '-atoms':
                if i + 1 >= len(argv):
                    raise Exception('Error: the \"' + argv[i] + '\" argument should be followed by a quoted string\n'
                                    '       which contains a space-delimited list of of a subset of atom types\n'
                                    '       you want to use from the original force-field.\n'
                                    '       Make sure you enclose the entire list in quotes.\n')
                type_subset = Set(argv[i + 1].strip('\"\'').strip().split())
                del argv[i:i + 2]
    
            elif argv[i] == '-name':
                if i + 1 >= len(argv):
                    raise Exception(
                        'Error: ' + argv[i] + ' flag should be followed by the name of the force-field\n')
                ffname = argv[i + 1]
                del argv[i:i + 2]
    
            elif argv[i] in ('-file', '-in-file'):
                if i + 1 >= len(argv):
                    raise Exception(
                        'Error: ' + argv[i] + ' flag should be followed by the name of a force-field file\n')
                filename_in = argv[i + 1]
                try:
                    file_in = open(filename_in, 'r')
                except IOError:
                    sys.stderr.write('Error: Unable to open file\n'
                                     '       \"' + filename_in + '\"\n'
                                     '       for reading.\n')
                    sys.exit(1)
                del argv[i:i + 2]

            elif argv[i] == '-dihedral-style':
                if i + 1 >= len(argv):
                    raise Exception(
                        'Error: ' + argv[i] + ' flag should be followed by either \"opls\" or \"fourier\"\n')
                dihedral_style_name = argv[i + 1]
                if dihedral_style_name == "fourier":
                    dihedral_style_link = "http://lammps.sandia.gov/doc/dihedral_fourier.html"
                if dihedral_style_name == "opls":
                    dihedral_style_link = "http://lammps.sandia.gov/doc/dihedral_opls.html"
                else:
                    raise Exception(
                        'Error: ' + argv[i] + ' ' + dihedral_style_name + ' not supported.\n')
                del argv[i:i + 2]
    
            elif argv[i] in ('-url', '-in-url'):
                import urllib2
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a URL pointing to\n'
                        
                                     ' a TINKER file containing force-field information.\n')
                url = argv[i + 1]
                try:
                    request = urllib2.Request(url)
                    file_in = urllib2.urlopen(request)
                except urllib2.URLError:
                    sys.stdout.write("Error: Unable to open link:\n" + url + "\n")
                    sys.exit(1)
                del argv[i:i + 2]

            elif argv[i] == '-hybrid':
                use_hybrid = True
                del argv[i:i + 1]
    
            elif (argv[i] == '-zeropad' or argv[i] == '-zero-pad'):
                if (i + 1 >= len(argv)) or (argv[i+1][1:] == '-'):
                    raise Exception(
                        'Error: ' + argv[i] + ' flag should be followed by a positive integer\n')
                zeropad_ffid = int(argv[i+1])
                del argv[i:i + 2]
    
            elif argv[i] in ('-help', '--help', '-?', '--?'):
                sys.stderr.write(doc_msg)
                sys.exit(0)
                del argv[i:i + 1]
    
            else:
                i += 1
    
        if len(argv) != 1:
            raise Exception('Error: Unrecongized arguments: ' + ' '.join(argv[1:]) +
                            '\n\n' + doc_msg)
    
        #sys.stderr.write("Reading parameter file...\n")
    
        lines = file_in.readlines()
    
        atom2charge = OrderedDict()  # lookup charge from atom type
        atom2mass = OrderedDict()  # lookup mass from atom type
        atom2vdw_e = OrderedDict()  # lookup Lennard-Jones "epsilon" parameter
        atom2vdw_s = OrderedDict()  # lookup Lennard-Jones "sigma" parameter
        atom2descr = OrderedDict()
        atom2ffid = OrderedDict()   # lookup force-field-ID from atom type
        # force-field-ID is an id number/string used to assign
        # bonds, angles, dihedrals, and impropers.
    
        bonds_by_type = OrderedDict()  # lookup bond parameters by force-field-ID
        angles_by_type = OrderedDict()  # lookup angle parameters by force-field-ID
        dihedrals_by_type = OrderedDict()  # lookup dihedral parameters by force-field-ID
        impropers_by_type = OrderedDict()  # lookup improper parameters by force-field-ID
        lines_ureybrad = []
        lines_biotype = []
    
        for iline in range(0, len(lines)):
            line = lines[iline]
            tokens = SplitQuotedString(line.strip(),
                                       comment_char='#')
    
            if (len(tokens) > 1) and (tokens[0] == 'atom'):
                tokens = map(RemoveOuterQuotes,
                             SplitQuotedString(line.strip(),
                                               comment_char=''))
                if (len(tokens) > 6):
                    if ((len(type_subset) == 0) or (tokens[1] in type_subset)):
                        atom2ffid[tokens[1]] = tokens[2]
                        #atom2mass[tokens[1]] = float(tokens[6])
                        # Some atoms in oplsaa.prm have zero mass. Unfortunately this
                        # causes LAMMPS to crash, even if these atoms are never used,
                        # so I give the mass a non-zero value instead.
                        atom2mass[tokens[1]] = max(float(tokens[6]), 1e-30)
                        atom2descr[tokens[1]] = tokens[4]
                        if tokens[4].find('(UA)') != -1:
                            contains_united_atoms = True
                else:
                    raise Exception('Error: Invalid atom line:\n' + line)
            elif (len(tokens) > 2) and (tokens[0] == 'charge'):
                if ((len(type_subset) == 0) or (tokens[1] in type_subset)):
                    atom2charge[tokens[1]] = float(tokens[2])
            elif (len(tokens) > 3) and (tokens[0] == 'vdw'):
                if ((len(type_subset) == 0) or (tokens[1] in type_subset)):
                    atom2vdw_e[tokens[1]] = float(tokens[3])  # "epsilon"
                    atom2vdw_s[tokens[1]] = float(tokens[2])  # "sigma"
            elif (len(tokens) > 4) and (tokens[0] == 'bond'):
                k = float(tokens[3])
                r0 = float(tokens[4])
                bonds_by_type[tokens[1].rjust(zeropad_ffid,'0'),
                              tokens[2].rjust(zeropad_ffid,'0')] = (k, r0)
            elif (len(tokens) > 5) and (tokens[0] == 'angle'):
                k = float(tokens[4])
                angle0 = float(tokens[5])
                angles_by_type[tokens[1].rjust(zeropad_ffid,'0'),
                               tokens[2].rjust(zeropad_ffid,'0'),
                               tokens[3].rjust(zeropad_ffid,'0')] = (k, angle0)
            elif (len(tokens) > 11) and (tokens[0] == 'torsion'):
                if dihedral_style_name == 'fourier':
                    # http://lammps.sandia.gov/doc/dihedral_fourier.html
                    m = (len(tokens) - 5) / 3
                    K = [0.0 for i in range(0, m)]
                    n = [0.0 for i in range(0, m)]
                    d = [0.0 for i in range(0, m)]
                    for i in range(0, m):
                        K[i] = float(tokens[5 + 3 * i])
                        d[i] = float(tokens[5 + 3 * i + 1])
                        n[i] = float(tokens[5 + 3 * i + 2])
                    dihedrals_by_type[tokens[1].rjust(zeropad_ffid,'0'),
                                      tokens[2].rjust(zeropad_ffid,'0'),
                                      tokens[3].rjust(zeropad_ffid,'0'),
                                      tokens[4].rjust(zeropad_ffid,'0')] = (K, n, d)
                elif dihedral_style_name == 'opls':
                    # http://lammps.sandia.gov/doc/dihedral_opls.html
                    K1 = float(tokens[5])
                    K2 = float(tokens[8])
                    K3 = float(tokens[11])
                    K4 = 0.0
                    if len(tokens) > 14:
                        K4 = float(tokens[14])
                    if ((float(tokens[6]) != 0.0) or (float(tokens[7]) != 1.0) or
                        (float(tokens[9]) not in (180.0, -180.0)) or (float(tokens[10]) != 2.0) or
                        (float(tokens[12]) != 0.0) or (float(tokens[13]) != 3.0) or
                        ((K4 != 0.0) and
                         ((len(tokens) <= 16) or
                          (float(tokens[15]) not in (180.0, -180.0)) or
                          (float(tokens[16]) != 4.0)))):
                        raise Exception("Error: This parameter file is incompatible with -dihedral-style \"" + dihedral_style_name + "\"\n" +
                                        "       (See line number " + str(iline + 1) + " of parameter file.)\n")
                    dihedrals_by_type[tokens[1].rjust(zeropad_ffid,'0'),
                                      tokens[2].rjust(zeropad_ffid,'0'),
                                      tokens[3].rjust(zeropad_ffid,'0'),
                                      tokens[4].rjust(zeropad_ffid,'0')] = (K1, K2, K3, K4)
                else:
                    assert(False)
    
            elif (len(tokens) > 7) and (tokens[0] == 'imptors'):
                k = float(tokens[5])
                angle0 = float(tokens[6])
                multiplicity = float(tokens[7])
                impropers_by_type[tokens[1].rjust(zeropad_ffid,'0'),
                                  tokens[2].rjust(zeropad_ffid,'0'),
                                  tokens[3].rjust(zeropad_ffid,'0'),
                                  tokens[4].rjust(zeropad_ffid,'0')] = (k / multiplicity, angle0)
            elif ((len(tokens) > 0) and (tokens[0] == 'biotype')):
                # I'm not sure what to do with these, so I'll store them for now and
                # append them as comments to the .lt file generated by the program.
                lines_biotype.append(line.rstrip())
            elif ((len(tokens) > 0) and (tokens[0] == 'ureybrad')):
                # I'm not sure what to do with these, so I'll store them for now and
                # append them as comments to the .lt file generated by the program.
                lines_ureybrad.append(line.rstrip())
            elif ((len(tokens) > 1) and (tokens[0] == 'radiusrule')):
                if tokens[1] == 'GEOMETRIC':
                    mixing_style = 'geometric'
                elif tokens[1] == 'ARITHMETIC':
                    mixing_style = 'arithmetic'
                else:
                    raise Exception("Error: Unrecognized mixing style: " + tokens[1] + ", found here:\n" +
                                    line)
            elif ((len(tokens) > 1) and (tokens[0] == 'epsilonrule')):
                if tokens[1] != 'GEOMETRIC':
                    raise Exception("Error: As of 2016-9-21, LAMMPS only supports GEOMETRIC mixing of energies\n" +
                                    "       This force field simply cannot be used with LAMMPS in a general way.\n" +
                                    "       One way around this is to manually change the \"epsilonrule\" back to\n" +
                                    "       GEOMETRIC, and limit the number of atom types considered by this\n" +
                                    "       program by using the -atoms \"LIST OF ATOMS\" argument,\n" +
                                    "       to only include the atoms you care about, and then explicitly\n" +
                                    "       define pair_coeffs for all possible pairs of these atom types.\n" +
                                    "       If this is a popular force-field, then lobby the LAMMPS developers\n" +
                                    "       to consider alternate mixing rules.\n\n" +
                                    "The offending line from the file is line number " + str(iline) + ":\n" +
                                    line + "\n")

        # Zero-pad the atom2ffid values so that they have the same number
        # of digits.  This is usually not necessary, but it can be helpful
        # to remove uncertainty about the meaning of '4*' which could
        # pattern match with '4', '4L', '47', '47L'...  If you replace '4'
        # with '04', '04*' becomes distinguishable from '47*'.
        # This can be useful if you want to augment the force field later,
        # (for example, adding additional atoms to the LOPLSAA variant of OPLSAA)

        for k in atom2ffid.keys():
            atom2ffid[k] = atom2ffid[k].rjust(zeropad_ffid, '0')

            # Horrible hack:  for LOPLSAA, uncomment the next 3 lines:
            #ki = atom2ffid[k].find('L')
            #if ki!=-1:
            #    atom2ffid[k] = atom2ffid[k].rjust(zeropad_ffid + len(atom2ffid[k]) - ki, '0')


        #sys.stderr.write(" done.\n")
        #sys.stderr.write("Converting to moltemplate format...\n")


        system_is_charged = False
        for atom_type in atom2charge:
            if atom2charge[atom_type] != 0.0:
                system_is_charged = True
        
        if system_is_charged:
            pair_style_name = "lj/cut/coul/long"
            pair_style_params = "10.0 10.0"
            kspace_style = "    kspace_style pppm 0.0001\n"
            pair_style_link = "http://lammps.sandia.gov/doc/pair_lj.html"
        else:
            pair_style_name = "lj/cut"
            pair_style_params = "10.0"
            kspace_style = ""
            pair_style_link = "http://lammps.sandia.gov/doc/pair_lj.html"

        pair_style_command = "    pair_style " + ("hybrid " if use_hybrid else "") + \
                             pair_style_name + " " + pair_style_params + "\n"
        
        sys.stdout.write("# This file was generated automatically using:\n")
        sys.stdout.write("# " + g_program_name + " " + " ".join(sys.argv[1:]) + "\n")
        if contains_united_atoms:
            sys.stdout.write("#\n"
                             "# WARNING: Many of these atoms are probably UNITED-ATOM (UA) atoms.\n"
                             "#          The hydrogen atoms have been absorbed into the heavy atoms, and the\n"
                             "#          force-field modified accordingly. Do not mix with ordinary atoms.\n")
        sys.stdout.write("#\n"
                         "# WARNING: The following 1-2, 1-3, and 1-4 weighting parameters were ASSUMED:\n")
        sys.stdout.write("#          " + special_bonds_command + "\n")
        sys.stdout.write(
            "#          (See http://lammps.sandia.gov/doc/special_bonds.html for details)\n")
        if len(lines_ureybrad) > 0:
            sys.stdout.write("#\n"
                             "# WARNING: All Urey-Bradley interactions have been IGNORED including:\n")
            sys.stdout.write("#             ffid1 ffid2 ffid3    K        r0\n# ")
            sys.stdout.write("\n# ".join(lines_ureybrad))
            sys.stdout.write("\n\n")
        sys.stdout.write("\n\n")
        sys.stdout.write(ffname + " {\n\n")
        
        sys.stdout.write("  # Below we will use lammps \"set\" command to assign atom charges\n"
                         "  # by atom type.  http://lammps.sandia.gov/doc/set.html\n\n")
        
        sys.stdout.write("  write_once(\"In Charges\") {\n")
        for atype in atom2mass:
            assert(atype in atom2descr)
            sys.stdout.write("    set type @atom:" + atype + " charge " + str(atom2charge[atype]) +
                             "  # \"" + atom2descr[atype] + "\"\n")
        sys.stdout.write("  } #(end of atom partial charges)\n\n\n")
        
        
        sys.stdout.write("  write_once(\"Data Masses\") {\n")
        for atype in atom2mass:
            sys.stdout.write("    @atom:" + atype + " " + str(atom2mass[atype]) + "\n")
        sys.stdout.write("  } #(end of atom masses)\n\n\n")
        
        
        sys.stdout.write("  # ---------- EQUIVALENCE CATEGORIES for bonded interaction lookup ----------\n"
                         "  #   Each type of atom has a separate ID used for looking up bond parameters\n"
                         "  #   and a separate ID for looking up 3-body angle interaction parameters\n"
                         "  #   and a separate ID for looking up 4-body dihedral interaction parameters\n"
                         "  #   and a separate ID for looking up 4-body improper interaction parameters\n"
                         #"  #   (This is because there are several different types of sp3 carbon atoms\n"
                         #"  #   which have the same torsional properties when within an alkane molecule,\n"
                         #"  #   for example.  If they share the same dihedral-ID, then this frees us\n"
                         #"  #   from being forced define separate dihedral interaction parameters\n"
                         #"  #   for all of them.)\n"
                         "  #   The complete @atom type name includes ALL of these ID numbers.  There's\n"
                         "  #   no need to force the end-user to type the complete name of each atom.\n"
                         "  #   The \"replace\" command used below informs moltemplate that the short\n"
                         "  #   @atom names we have been using abovee are equivalent to the complete\n"
                         "  #   @atom names used below:\n\n")
        
        for atype in atom2ffid:
            ffid = atype + "_ffid" + atom2ffid[atype]
            sys.stdout.write("  replace{ @atom:" + atype +
                             " @atom:" + atype + "_b" + atom2ffid[atype] + "_a" + atom2ffid[atype] + "_d" + atom2ffid[atype] + "_i" + atom2ffid[atype] + " }\n")
        sys.stdout.write("\n\n\n\n")
        
        
        sys.stdout.write("  # --------------- Non-Bonded interactions: ---------------------\n"
                         "  # " + pair_style_link + "\n"
                         "  # Syntax:\n"
                         "  # pair_coeff    AtomType1    AtomType2   " +
                         ("PairStyleName  " if use_hybrid else "") +
                         "parameters...\n\n")
            
        sys.stdout.write("  write_once(\"In Settings\") {\n")
        for atype in atom2vdw_e:
            assert(atype in atom2vdw_s)
            if not (atype in atom2ffid):
                continue
        
            sys.stdout.write("    pair_coeff " +
                             "@atom:" + atype + "_b" + atom2ffid[atype] + "_a" + atom2ffid[
                                 atype] + "_d" + atom2ffid[atype] + "_i" + atom2ffid[atype] + " "
                             "@atom:" + atype + "_b" + atom2ffid[atype] + "_a" + atom2ffid[atype] + "_d" + atom2ffid[atype] + "_i" + atom2ffid[atype] + " " +
                             (pair_style_name if use_hybrid else "") +
                             " " + str(atom2vdw_e[atype]) +
                             " " + str(atom2vdw_s[atype]) + "\n")
        sys.stdout.write("  } #(end of pair_coeffs)\n\n\n\n")
        
        
        sys.stdout.write("  # ------- Bonded Interactions: -------\n"
                         "  # " + bond_style_link + "\n"
                         "  # Syntax:  \n"
                         "  # bond_coeff BondTypeName  " +
                         ("BondStyleName  " if use_hybrid else "") +
                         "parameters...\n\n")
        
        sys.stdout.write("  write_once(\"In Settings\") {\n")
        for btype in bonds_by_type:
            ffid1 = btype[0] if btype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = btype[1] if btype[1] != ("0"*zeropad_ffid) else "X"
            (k, r0) = bonds_by_type[btype]
            sys.stdout.write("    bond_coeff @bond:" + ffid1 + "_" + ffid2 + " " +
                             (bond_style_name if use_hybrid else "") +
                             " " + str(k) + " " + str(r0) + "\n")
        sys.stdout.write("  } #(end of bond_coeffs)\n\n")
        
        sys.stdout.write("  # Rules for assigning bond types by atom type:\n"
                         "  # BondTypeName       AtomType1         AtomType2\n"
                         "  #   (* = wildcard)\n\n")
        
        sys.stdout.write("  write_once(\"Data Bonds By Type\") {\n")
        for btype in bonds_by_type:
            ffid1 = btype[0] if btype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = btype[1] if btype[1] != ("0"*zeropad_ffid) else "X"
            sys.stdout.write("    @bond:" + ffid1 + "_" + ffid2)
            ffid1 = "@atom:*_b" + btype[0] + \
                "*_a*_d*_i*" if btype[0] != ("0"*zeropad_ffid) else "@atom:*"
            ffid2 = "@atom:*_b" + btype[1] + \
                "*_a*_d*_i*" if btype[1] != ("0"*zeropad_ffid) else "@atom:*"
            sys.stdout.write(" " + ffid1 + " " + ffid2 + "\n")
        sys.stdout.write("  } #(end of bonds by type)\n\n\n\n\n")
        
        
        sys.stdout.write("  # ------- Angle Interactions: -------\n"
                         "  # " + angle_style_link + "\n"
                         "  # Syntax:  \n"
                         "  # angle_coeff AngleTypeName  "+
                         ("AngleStyleName  " if use_hybrid else "") +
                         "parameters...\n\n")
        
        sys.stdout.write("  write_once(\"In Settings\") {\n")
        for atype in angles_by_type:
            ffid1 = atype[0] if atype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = atype[1] if atype[1] != ("0"*zeropad_ffid) else "X"
            ffid3 = atype[2] if atype[2] != ("0"*zeropad_ffid) else "X"
            (k, angle0) = angles_by_type[atype]
            sys.stdout.write("    angle_coeff @angle:" + ffid1 + "_" + ffid2 + "_" + ffid3 + " " +
                             (angle_style_name if use_hybrid else "") +
                             " " + str(k) + " " + str(angle0) + "\n")
        sys.stdout.write("  } #(end of angle_coeffs)\n\n")
        
        
        sys.stdout.write("  # Rules for creating angle interactions according to atom type:\n"
                         "  #   AngleTypeName     AtomType1       AtomType2          AtomType3\n"
                         "  #   (* = wildcard)\n\n")
        
        sys.stdout.write("  write_once(\"Data Angles By Type\") {\n")
        for atype in angles_by_type:
            ffid1 = atype[0] if atype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = atype[1] if atype[1] != ("0"*zeropad_ffid) else "X"
            ffid3 = atype[2] if atype[2] != ("0"*zeropad_ffid) else "X"
            sys.stdout.write("    @angle:" + ffid1 + "_" + ffid2 + "_" + ffid3)
            ffid1 = "@atom:*_b*_a" + atype[0] + \
                "*_d*_i*" if atype[0] != ("0"*zeropad_ffid) else "@atom:*"
            ffid2 = "@atom:*_b*_a" + atype[1] + \
                "*_d*_i*" if atype[1] != ("0"*zeropad_ffid) else "@atom:*"
            ffid3 = "@atom:*_b*_a" + atype[2] + \
                "*_d*_i*" if atype[2] != ("0"*zeropad_ffid) else "@atom:*"
            sys.stdout.write(" " + ffid1 + " " + ffid2 + " " + ffid3 + "\n")
        sys.stdout.write("  } #(end of angles by type)\n\n\n\n\n")
        
        
        sys.stdout.write("  # ----------- Dihedral Interactions: ------------\n"
                         "  # " + dihedral_style_link + "\n"
                         "  # Syntax:\n"
                         "  # dihedral_coeff DihedralTypeName  " +
                         ("DihedralStyleName  " if use_hybrid else "") +
                         "parameters...\n\n")
        
        sys.stdout.write("  write_once(\"In Settings\") {\n")
        for dtype in dihedrals_by_type:
            ffid1 = dtype[0] if dtype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = dtype[1] if dtype[1] != ("0"*zeropad_ffid) else "X"
            ffid3 = dtype[2] if dtype[2] != ("0"*zeropad_ffid) else "X"
            ffid4 = dtype[3] if dtype[3] != ("0"*zeropad_ffid) else "X"
            sys.stdout.write("    dihedral_coeff @dihedral:" +
                             ffid1 + "_" + ffid2 + "_" + ffid3 + "_" + ffid4 + " " +
                             (dihedral_style_name if use_hybrid else "") +
                             " ")
            if dihedral_style_name == 'fourier':
                # http://lammps.sandia.gov/doc/dihedral_fourier.html
                (K, n, d) = dihedrals_by_type[dtype]
                m = len(K)
                assert((m == len(n)) and (m == len(d)))
                sys.stdout.write(str(m))
                for i in range(0, m):
                    sys.stdout.write(" " + str(K[i]) +
                                     " " + str(n[i]) + " " + str(d[i]))
                sys.stdout.write("\n")
            elif dihedral_style_name == 'opls':
                # http://lammps.sandia.gov/doc/dihedral_opls.html
                (K1, K2, K3, K4) = dihedrals_by_type[dtype]
                sys.stdout.write(str(K1) + " " + str(K2) + " " +
                                 str(K3) + " " + str(K4) + "\n")
            else:
                assert(False)
        sys.stdout.write("  } #(end of dihedral_coeffs)\n\n")
        
        
        sys.stdout.write("  # Rules for creating dihedral interactions according to atom type:\n"
                         "  #   DihedralTypeName     AtomType1     AtomType2     AtomType3     AtomType4\n"
                         "  #   (* = wildcard)\n\n")
        
        sys.stdout.write("  write_once(\"Data Dihedrals By Type\") {\n")
        for dtype in dihedrals_by_type:
            ffid1 = dtype[0] if dtype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = dtype[1] if dtype[1] != ("0"*zeropad_ffid) else "X"
            ffid3 = dtype[2] if dtype[2] != ("0"*zeropad_ffid) else "X"
            ffid4 = dtype[3] if dtype[3] != ("0"*zeropad_ffid) else "X"
            sys.stdout.write("    @dihedral:" +
                             ffid1 + "_" + ffid2 + "_" +
                             ffid3 + "_" + ffid4)
            ffid1 = "@atom:*_b*_a*_d" + dtype[0] + \
                "*_i*" if dtype[0] != ("0"*zeropad_ffid) else "@atom:*"
            ffid2 = "@atom:*_b*_a*_d" + dtype[1] + \
                "*_i*" if dtype[1] != ("0"*zeropad_ffid) else "@atom:*"
            ffid3 = "@atom:*_b*_a*_d" + dtype[2] + \
                "*_i*" if dtype[2] != ("0"*zeropad_ffid) else "@atom:*"
            ffid4 = "@atom:*_b*_a*_d" + dtype[3] + \
                "*_i*" if dtype[3] != ("0"*zeropad_ffid) else "@atom:*"
            
            sys.stdout.write(" " + ffid1 + " " + ffid2 +
                             " " + ffid3 + " " + ffid4 + "\n")
        sys.stdout.write("  } #(end of dihedrals by type)\n\n\n\n\n")
        
        
        sys.stdout.write("  # ---------- Improper Interactions: ----------\n"
                         "  # " + improper_style_link + "\n"
                         "  # Syntax:\n"
                         "  # improper_coeff ImproperTypeName  " +
                         ("ImproperStyleName  " if use_hybrid else "") +
                         "parameters\n\n")
        
        sys.stdout.write("  write_once(\"In Settings\") {\n")
        for itype in impropers_by_type:
            ffid1 = itype[0] if itype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = itype[1] if itype[1] != ("0"*zeropad_ffid) else "X"
            ffid3 = itype[2] if itype[2] != ("0"*zeropad_ffid) else "X"
            ffid4 = itype[3] if itype[3] != ("0"*zeropad_ffid) else "X"
            (k, angle0) = impropers_by_type[itype]
            sys.stdout.write("    improper_coeff @improper:" +
                             ffid1 + "_" + ffid2 + "_" + ffid3 + "_" + ffid4 + " " +
                             (improper_style_name if use_hybrid else "") +
                             " " + str(k) + " " + str(angle0) + "\n")
        sys.stdout.write("  } #(end of improper_coeffs)\n\n")
        
        
        sys.stdout.write("  # Rules for creating improper interactions according to atom type:\n"
                         "  #   ImproperTypeName   AtomType1    AtomType2       AtomType3       AtomType4\n"
                         "  #   (* = wildcard)\n")
        
        sys.stdout.write("  write_once(\"Data Impropers By Type (opls_imp.py)\") {\n")
        for itype in impropers_by_type:
            ffid1 = itype[0] if itype[0] != ("0"*zeropad_ffid) else "X"
            ffid2 = itype[1] if itype[1] != ("0"*zeropad_ffid) else "X"
            ffid3 = itype[2] if itype[2] != ("0"*zeropad_ffid) else "X"
            ffid4 = itype[3] if itype[3] != ("0"*zeropad_ffid) else "X"
            sys.stdout.write("    @improper:" +
                             ffid1 + "_" + ffid2 + "_" +
                             ffid3 + "_" + ffid4)
            ffid1 = "@atom:*_b*_a*_d*_i" + itype[0]+"*" if itype[0] != ("0"*zeropad_ffid) else "@atom:*"
            ffid2 = "@atom:*_b*_a*_d*_i" + itype[1]+"*" if itype[1] != ("0"*zeropad_ffid) else "@atom:*"
            ffid3 = "@atom:*_b*_a*_d*_i" + itype[2]+"*" if itype[2] != ("0"*zeropad_ffid) else "@atom:*"
            ffid4 = "@atom:*_b*_a*_d*_i" + itype[3]+"*" if itype[3] != ("0"*zeropad_ffid) else "@atom:*"
            sys.stdout.write(" " + ffid1 + " " + ffid2 +
                             " " + ffid3 + " " + ffid4 + "\n")
        sys.stdout.write("  } #(end of impropers by type)\n\n\n\n\n")
        
        sys.stdout.write("  # --------   (descriptive comment)   --------\n")
        sys.stdout.write("  # ---- biologically relevant atom types: ----\n  # ")
        sys.stdout.write("\n  # ".join(lines_biotype))
        sys.stdout.write("\n  # ----------   (end of comment)   ----------\n")
        sys.stdout.write("\n\n\n\n")
        
        
        sys.stdout.write("  # LAMMPS supports many different kinds of bonded and non-bonded\n"
                         "  # interactions which can be selected at run time.  Eventually\n"
                         "  # we must inform LAMMPS which of them we will need.  We specify\n"
                         "  # this in the \"In Init\" section: \n\n")
        
        sys.stdout.write("  write_once(\"In Init\") {\n")
        sys.stdout.write("    units real\n")
        sys.stdout.write("    atom_style full\n")
        sys.stdout.write("    bond_style " +
                         ("hybrid " if use_hybrid else "") +
                         bond_style_name + "\n")
        sys.stdout.write("    angle_style " +
                         ("hybrid " if use_hybrid else "") +
                         angle_style_name + "\n")
        sys.stdout.write("    dihedral_style " +
                         ("hybrid " if use_hybrid else "") +
                         dihedral_style_name + "\n")
        sys.stdout.write("    improper_style " +
                         ("hybrid " if use_hybrid else "") +
                         improper_style_name + "\n")
        sys.stdout.write(pair_style_command)
        sys.stdout.write("    pair_modify mix " + mixing_style + "\n")
        sys.stdout.write("    " + special_bonds_command + "\n")
        sys.stdout.write(kspace_style)
        sys.stdout.write("  } #end of init parameters\n\n")
        
        sys.stdout.write("  # Note: We use \"hybrid\" styles in case the user later wishes to\n"
                         "  #       combine the molecules built using this force-field with other\n"
                         "  #       molecules that use other styles.  (This is not necessarily\n"
                         "  #       a good idea, but LAMMPS and moltemplate both allow it.)\n"
                         "  #       For more information:\n"
                         "  #       http://lammps.sandia.gov/doc/pair_hybrid.html\n"
                         "  #       http://lammps.sandia.gov/doc/bond_hybrid.html\n"
                         "  #       http://lammps.sandia.gov/doc/angle_hybrid.html\n"
                         "  #       http://lammps.sandia.gov/doc/dihedral_hybrid.html\n"
                         "  #       http://lammps.sandia.gov/doc/improper_hybrid.html\n\n\n")
        
        sys.stdout.write("}  # " + ffname + "\n\n")
        
        
        #sys.stderr.write(" done.\n")
        
        if filename_in != "":
            file_in.close()
    
    
    
        
    except Exception as err:
        sys.stderr.write('\n\n' + str(err) + '\n')
        sys.exit(1)


if __name__ == '__main__':
    main()
