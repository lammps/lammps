#! /usr/bin/env python
# Author: Andrew Jewett (jewett.aij at g mail)
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2017, California Institute of Technology
# All rights reserved.

"""
This standalone python script can be used to convert force-field data 
in FRC files (a.k.a. "MSI", "Accelrys", "BIOSYM", "DISCOVERY" files)
...into MOLTEMPLATE/LAMMPS compatible format (.LT files).

Once converted into moltemplate (.LT) format, users can use these files with 
MOLTEMPLATE to prepare LAMMPS simulations of molecules using these force fields 
(without needing any additional software such as msi2lmp).

There are several examples of MSI files in the "tools/msi2lmp/frc_files/"
directory which is distributed with LAMMPS.

Limitations:

Currently (2017-10) this script ignores the "template" information in .FRC files.
When defining a new type of molecule, the user must carefully choose the
complete atom type for each type of atom in the molecule.  In other words,
MOLTEMPLATE will not attempt to determine (from local context) whether 
a carbon atom somewhere in your molecule happens to be an SP3 carbon 
(ie. "c4" in the COMPASS force-field), or an aromatic carbon ("c3a"), 
or something else (for example).  This information is typically contained 
in the "templates" section of these files, and this script currently ignores
that information.  Instead, the user must determine which type of carbon atom
it is manually, for all of the carbon atoms in that kind of molecule.
(This only needs to be done once per molecule definition.
 Once a type of molecule is defined, it can be copied indefinitely.)

"""


__author__ = 'Andrew Jewett'
__version__ = '0.2.1'
__date__ = '2017-10-15'


import sys
import os

from collections import defaultdict, OrderedDict
from operator import itemgetter


g_program_name = __file__.split('/')[-1]


doc_msg = \
    "Typical Usage:\n\n" + \
    "   " + g_program_name + " -name COMPASS < compass_published.frc > compass.lt\n\n" + \
    "   where \"compass_published.frc\" is a force-field file in MSI format.\n" + \
    "         \"comass.lt\" is the corresponding file converted to moltemplate format\n" + \
    "   and   \"COMPASS\" is the name that future moltemplate users will use to refer\n" + \
    "         to this force-field (optional).\n" + \
    "Optional Arguments\n" + \
    "   -name FORCEFIELDNAME # Give the force-field a name\n" + \
    "   -file FILE_NAME      # Read force field parameters from a file\n" + \
    "   -url URL             # Read force field parameters from a file on the web\n" + \
    "   -atoms \"QUOTED LIST\" # Restrict output to a subset of atom types\n" + \
    "  Sometimes an FRC file contains multiple versions.  In that case,\n"+\
    "  you can select between them using these optional arguments:\n"+\
    "   -pair-style \"PAIRSTYLE ARGS\" # LAMMPS pair style and cutoff arg(s)\n" + \
    "   -bond-style BONDSTYLE  # desired LAMMPS bond style (default: \"class2\")\n" + \
    "   -angle-style ANGLESTYLE  # desired LAMMPS angle style\n" + \
    "   -dihedral-style DIHEDRALSTYLE  # desired LAMMPS dihedral style\n" + \
    "   -improper-style IMPROPERSTYLE  # desired LAMMPS improper style\n" + \
    "   -hbond-style \"HBONDTYLE ARGS\" # LAMMPS hydrogen-bond style and args\n"


#   "   -auto                # Consider auto_equivalences in the .frc file \n"+\



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

# It seems like there are no ordered sets in python, (a set that remembers the
# order that you added elements), so I built one by wrapping OrderedDict()

class MyOrderedSet(object):
    def __init__(self, l):
        self.d = OrderedDict()
        for x in l:
            self.d[x] = True
    def __add__(self, x):
        self.d[x] = True
    def __delitem__(self, x):
        del self.d[x]
    def __contains__(self, x):
        return x in self.d
    def __iter__(self):
        self.p = iter(self.d)
        return self
    def __next__(self):
        return next(self.p)
    # the following wrappers might be necessary for python2/3 compatibility:
    def add(self, x):
        self.__add__(x)
    def del_item(self, x):
        self.__del_item__(x)
    def iter(self):
        return self.__iter__()
    def next(self):
        return self.__next__()
    # no need to bother with set unions and intersections


def NSplitQuotedString(string,
                       nmax,
                       quotes,
                       delimiters=' \t\r\f\n',
                       escape='\\',
                       comment_char='#'):
    """
    Split a quoted & commented string into at most "nmax" tokens (if nmax>0),
    where each token is separated by one or more delimeter characters
    in the origingal string, and quoted substrings are not split,
    This function returns a list of strings.  Once the string is split Nmax
    times, any remaining text will be appended to the last entry of the list.
    Comments are stripped from the string before splitting begins.
    """
    tokens = []
    token = ''
    reading_token = True
    escaped_state = False
    quote_state = None
    for c in string:

        if (c in comment_char) and (not escaped_state) and (quote_state == None):
            if len(token) > 0:
                tokens.append(token)
            return tokens

        elif (c in delimiters) and (not escaped_state) and (quote_state == None):
            if reading_token:
                if (nmax == 0) or (len(tokens) < nmax-1):
                    if len(token) > 0:
                        tokens.append(token)
                    token = ''
                    reading_token = False
                else:
                    token += c
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

    if len(token) > 0:
        tokens.append(token)
    return tokens




def SplitQuotedString(string,
                      quotes='\'\"',
                      delimiters=' \t\r\f\n',
                      escape='\\',
                      comment_char='#'):

    return NSplitQuotedString(string,
                              0,
                              quotes,
                              delimiters,
                              escape,
                              comment_char)




def RemoveOuterQuotes(text, quotes='\"\''):
    if ((len(text) >= 2) and (text[0] in quotes) and (text[-1] == text[0])):
        return text[1:-1]
    else:
        return text


def SortByEnds(l_orig):
    """
    Convenient to have a one-line macro for swapping list order if first>last
    """
    l = [x for x in l_orig]
    if l[0] > l[-1]:
        l.reverse()
    return l



#def Repl(tokens, a, b):
#    return [(b if x==a else x) for x in tokens]

def DecodeAName(s):
    if s.find('auto') == 0:
        s = s[4:]
    if s == 'X': # special case: deal with strings like  'X'
        return '*'
    return s


def EncodeAName(s):
    """
    Handle * characters in MSI atom names
    """

    if s.find('auto') == 0:
        s = s[4:]
    # If the atom name begins with *, then it is a wildcard
    if s[:1] == '*': # special case: deal with strings like  *7
        return 'X'   # These have special meaning.  Throw away the integer.
                     # (and replace the * with an X)

    # If the * character occurs later on in the atom name, then it is actually
    # part of the atom's name.  (MSI force fields use many strange characters in
    # atom names.)  Here we change the * to \* to prevent the atom name from
    # being interpreted as a wild card in the rules for generating bonds,
    # angles, dihedrals, and impropers.

    return s.replace('*','star').replace('\'','prime').replace('"','dblpr')
                                 # '*' is reserved for wildcards in moltemplate
                                 # 'star' is a string that is unused in any 
                                 # of the force fields I have seen so far.
                                 # Similarly quote characters (' and ") confuse
                                 # moltemplate, so we replace them with something else.

    # The following approach doesn't work (mistakenly thinks '\*' = wildcard)
    #return s.replace('*','\\*')  # this prevents ttree_lex.MatchesAll()
    #                             # from interpreting the '*' as a wildcard
    

def DetermineAutoPriority(anames):
    """
    Given a list of atom names (including wildcards), generate a number 
    indicating the priority the interaction between these atoms should have:
    Scan through list of strings anames, looking for patterns of the form
    *n
    where n is an integer.
    Make sure this pattern only appears once and return n to the caller.
    (These patterns are used by MSI software when using "auto_equivalences"
     to look up force field parameters for bonded interactions.
     The higher the integer, the lower the priority.
     For details, see "Forcefield based simulations" PDF, Cerius2, p 87)
    Ordinary wildcards ('*' characters not followed by integers) have the
    lowest possible priority.  (Each time a '*' string appears in the
    list of arguments, the priority value increases by HUGE_VAL.)
    """

    # This is terrible code.

    n = -1.0
    num_blank_wildcards = 0
    for a in anames:
        # Sometimes the first atom name contains the prefix 'auto'.  Remove this
        if a.find('auto') == 0:
            a = a[4:]
        if a[:1] == '*':
        #if a[:1] == 'X':
            if len(a) > 1:
                if n == -1.0:
                    n = float(a[1:])
                elif n != float(a[1:]):
                    # Make sure if present, the number appears only once in the list of atom names
                    raise InputError('Error: Inconsistent priority integers in the following interaction:\n'
                                     '      ' + ' '.join(anames) + '\n')
            else:
                num_blank_wildcards += 1

    # A "blank" wildcard (not followed by a number eg '*') has a very low priority
    # Give it a high number, because this corresponds to low priority.  Very confusing
    # For details, see "Forcefield based simulations" PDF, Cerius2, p 87)
    HUGE_VAL = 1.0e5
    return n + num_blank_wildcards*HUGE_VAL




#def DeterminePriority(is_auto,
#                      anames,
#                      version):
#    """
#    Determine the priority of an interaction from 
#    1) whether or not it is an "auto" interaction
#    2) what is the force-field "version" (a number)
#    3) what are the names of the atoms (for auto_equivalences only,
#       some atom "names" are wildcards followed by integers. use the integer)
#    """
#
#    if is_auto:
#        n = DetermineAutoPriority(anames)
#        return (is_auto, n)
#    else:
#        return (is_auto, -version)


def DetermineNumericPriority(is_auto,
                             anames,
                             version):
    """
    Determine the priority of an interaction from 
    2) what is the force-field "version" (a number)
    3) what are the names of the atoms (for auto_equivalences only,
       some atom "names" are wildcards followed by integers. use the integer)
    """

    if is_auto:
        n = DetermineAutoPriority(anames)
        return n  # low priority integers <--> high priority  ()
    else:
        return -float(version)  # later version numbers <--> higher priority
                                # (multiplying by -1 compensates for this)
                                # Note: this means auto interactions always have
                                # lower priority because their numeric priority
                                # will be a positive number.  Otherwise the
                                # numeric priority will be a negative number
                                # (...corresponding to a higher priority
                                #  I don't like this complicated priority system
                                #  but I didn't invent it.  It's not my fault.)


def IsAutoAtom(atom_name):
    return atom_name[-1:] == '_'


#def PossibleAutoAtom(atom_name):
#    """ Auto-equivalences are alternate atom names used in "auto" 
#    interactions. (These are low priority interactions used as a
#    last resort when the interaction parameters could not be located
#    by any other means).  Each atom is given an alternate name which
#    is used in this kind of interaction.  These names typically end
#    '_' followed by an optional integer.  Example "auto" atom names
#    are 'c3m_' and 'c=_3'.  Unfortunately some ordinary atom names
#    also end in an integer preceeded by a _ character. But they
#    never end in a '_' character.  Here we check for both."""
#
#    i = atom_name.rfind('_')
#    if (i == -1) or str.isdigit(atom_name[i:]):
#        return True
#    return False   



def IsAutoInteraction(interaction_name):
    return interaction_name.find('auto') == 0


#def IsAutoInteraction(interaction_name):
#    anames = ExtractAtomNames(interaction_name)
#    for a in anames:
#        if IsAutoAtom(a):
#            return True
#        if not PossibleAutoAtom(a):
#            return False
#    return True




def EncodeInteractionName(anames,
                          is_auto = False):
    if is_auto == False:
        is_auto = False
        # Is the line containing anames from an "_auto" section of 
        # the FRC file?  (I am trying to infer this from whether or 
        # not any of the atom names are followed by the '_' character.)
        for s in anames:  
            if IsAutoAtom(s):
                is_auto = True
    if is_auto:
        priority = DetermineAutoPriority(anames)
        # (If an atom name is a wildcard '*' followed by 
        #  an integer, DetermineAutoPriority() will return 
        #  that integer.  Otherwise it will return '')
        #return str(priority)+'auto'+','.join(anames)
        return 'auto'+','.join(anames)

    return ','.join(anames)



def ExtractANames(interaction_name):
    if IsAutoInteraction(interaction_name):
        return interaction_name[4:].split(',')
    return interaction_name.split(',')



def OOPImproperNameSort(aorig):
    assert(len(aorig) == 4)
    atom_names = map(EncodeAName, aorig)
    if atom_names[0] < atom_names[3]:
        return (atom_names, [0,1,2,3])
    else:
        return ([atom_names[3],
                 atom_names[1],
                 atom_names[2],
                 atom_names[0]],
                [3,1,2,0])


def Class2ImproperNameSort(aorig):
    """
    This function takes a list of 4 strings as an argument representing 4 atom
    names for atoms participating in an "improper" ("wilson-out-of-plane")
    interaction.  This function assumes the second atom is the central ("hub") 
    atom in the interaction, and it sorts the remaining atoms names.
    This function also replaces any occurence of \"*\" with \"X\".
    The new list is returned to the caller, along with the permutation.
    """
    assert(len(aorig) == 4)
    atom_names = [a for a in map(EncodeAName, aorig)]
    z = [x for x in zip([atom_names[0], atom_names[2], atom_names[3]],
                        [0,2,3])]
    z.sort()
    l = [z[0][0], atom_names[1], z[1][0], z[2][0]]
    p = [z[0][1], 1, z[1][1], z[2][1]]
    return (l, p)



def Parity(p):
    """ compute the parity of a permutation 
        (credit: "Weeble") 
    """
    permutation = list(p)
    length = len(permutation)
    elements_seen = [False] * length
    cycles = 0
    for index, already_seen in enumerate(elements_seen):
        if already_seen:
            continue
        cycles += 1
        current = index
        while not elements_seen[current]:
            elements_seen[current] = True
            current = permutation[current]
    return (length-cycles) % 2 == 0



def ImCrossTermID(atom_names):
    """
    # From a list of 4 atom names, corresponding two a pair
    # of angles between atoms# 3,2,1 and 3,2,4,
    # and replaces the list of atoms with a canonical tuple
    # which eliminates order ambiguity.
    # If you swap the first and last atom (#1 and #4), then
    # the -pair- of angles is the same.  Hence if atom #1
    # has a name which is lexicographically less than atom #4,
    # swap atoms 1 and 4.
    """
    if atom_names[0] <= atom_names[3]:
        return (atom_names[0]+','+atom_names[1]+','+
                atom_names[2]+','+atom_names[3])
    else:
        return (atom_names[3]+','+atom_names[1]+','+
                atom_names[2]+','+atom_names[0])




def AtomsMatchPattern(anames, pattern):
    """ 
    Check whether the list of atom names "anames" matches "pattern"
    (Both arguments are lists of strings, but some of the strings 
    in pattern may contain wildcard characters followed by 
    "priority" numbers.  Matches with lower priority numbers are
    given preference whenever multiple distinct matches are found.
    (Note: This function does not check patterns in reverse order.)
    """
    #sys.stderr.write('DEBUG: checking whether '+str(anames)+' matches '+str(pattern)+'\n')
    assert(len(anames) == len(pattern))
    matched = True
    for d in range(0, len(pattern)):
        if (pattern[d] == anames[d]) or (pattern[d][0] == '*'):
            if pattern[d][0] == '*':
                priority = int(pattern[d][1:])
            else:
                priority = 0
        else:
            matched = False
    if matched:
        #sys.stderr.write('DEBUG: '+str(anames)+' matches '+str(pattern)+'\n')
        return priority
    else:
        return None


def LookupBondLength(a1, a2,
                     atom2equiv_bond,
                     bond2r0,
                     atom2auto_bond,
                     bond2r0_auto):
    """ 
    Try to find bond parameters between atoms whose original
    atom names (without equivalences) are a1 and a2.
    Then return both the equilibrium bond length for that bond,
    as well as the equivalent atom names used to lookup that bond.
    (These could be stored in either atom2equiv_bond or atom2auto_bond.)
    If a match was not found, return None.
    """
    return_val = None
    anames = (atom2equiv_bond[a1], atom2equiv_bond[a2])
    bond_name = EncodeInteractionName(SortByEnds(anames))
    if bond_name in bond2r0:
        return_val = (bond2r0[bond_name],
                      [anames[0], anames[1]],
                      False)
    # If no bond between these atoms is defined, 
    # check the bonds in the _auto section(s)
    # This is a lot messier.
    elif ((a1 in atom2auto_bond) and (a2 in atom2auto_bond)):
        anames = [atom2auto_bond[a1], atom2auto_bond[a2]]
        # Because _auto interactions can contain wildcards,
        # there can be multiple entries in bond2r0_auto[]
        # for the same list of atom names, and we have to
        # consider all of them, and pick the one with the
        # most priority (ie. whose priority number is lowest).
        # (Note: The MSI file format uses low priority numbers
        #  to indicate high priority.  Somewhat confusing.
        #  For details, see "Forcefield based simulations" PDF, Cerius2, p 87)
        HUGE_VAL = 2000000000
        best_priority = HUGE_VAL
        pattern = ['','']
        for (pattern[0],pattern[1]), r0 in bond2r0_auto.items():
            priority = AtomsMatchPattern(anames, pattern)
            if (priority != None) and (priority < best_priority):
                best_priority = priority
                return_val = (r0, anames, True)
            # try again with the atom type names in reverse order
            priority = AtomsMatchPattern([anames[1],anames[0]], pattern)
            if ((priority != None) and
                (priority < best_priority)):  #(note: low priority numbers = high priority)
                best_priority = priority
                return_val = (r0, anames, True)
        #if return_val != None:
        #    sys.stderr.write('DEBUG: For atoms '+str((a1,a2))+' ... bond_length, batom_names = '+str(return_val)+'\n')
    return return_val








def LookupBondAngle(a1, a2, a3,
                    atom2equiv_angle,
                    angle2theta0_or,
                    atom2auto_angle,
                    angle2theta0_auto_or):
    """ 
    Try to find angle parameters between atoms whose original atom
    names (without equivalences) are a1, a2, and a3.  Then return
    both the equilibrium rest angle for that 3body interaction
    as well as the equivalent atom names used to look it up. (These
    could be stored in either atom2equiv_angle or atom2auto_angle.)
    If a match was not found, return None.
    """
    return_val = None
    anames = (atom2equiv_angle[a1], atom2equiv_angle[a2], atom2equiv_angle[a3])
    angle_name = EncodeInteractionName(SortByEnds(anames))
    if angle_name in angle2theta0_or:
        return_val = (angle2theta0_or[angle_name],
                      [anames[0], anames[1], anames[2]],
                      False)

    # If no angle between these atoms is defined, 
    # check the angles in the _auto section(s)
    # This is a lot messier.
    elif ((a1 in atom2auto_angle[0]) and
          (a2 in atom2auto_angle[1]) and
          (a3 in atom2auto_angle[2])):

        anames = [atom2auto_angle[0][a1],
                  atom2auto_angle[1][a2],
                  atom2auto_angle[2][a3]]
        #sys.stderr.write('DEBUG: LookupBondAngle(): a1,a2,a3=('+
        #                 a1+','+a2+','+a3+'), anames='+str(anames)+'\n')

        # Because _auto interactions can contain wildcards,
        # there can be multiple entries in angle2theta0_auto_or[]
        # for the same list of atom names, and we have to
        # consider all of them, and pick the one with the
        # most priority (ie. whose priority number is lowest).
        # (Note: The MSI file format uses low priority numbers
        #  to indicate high priority.  Somewhat confusing.)
        HUGE_VAL = 2000000000
        best_priority = HUGE_VAL  # (ie. low priority)
        pattern = ['','','']
        for (pattern[0],pattern[1],pattern[2]), theta0 in angle2theta0_auto_or.items():
            priority = AtomsMatchPattern(anames, pattern)
            if ((priority != None) and
                (priority < best_priority)):  #(note: low priority numbers = high priority)
                best_priority = priority
                return_val = (theta0, anames, True)
            # try again with the atom type names in reverse order
            priority = AtomsMatchPattern([anames[2],anames[1],anames[0]], pattern)
            if (priority != None) and (priority < best_priority):
                best_priority = priority
                return_val = (theta0, anames, True)
        #if return_val != None:
        #    sys.stderr.write('DEBUG: For atoms '+str((a1,a2,a3))+' ... rest_angle, anames = '+str(return_val)+'\n')
    return return_val


                              





def Equivalences2ffids(lines_equivalences,
                       atom_types,
                       atom2equiv_pair,
                       atom2equiv_bond,
                       atom2equiv_angle,
                       atom2equiv_dihedral,
                       atom2equiv_improper):
    """
    This function reads a list of lines containing "equivalences" and
    "auto_equivalences" from an MSI-formatted .FRC file.
    Then, for each atom type, it generates a long string which includes the 
    original atom type name as well as all of the equivalences it belongs to.
    Later on, when it is time to generate angles, dihedrals, or impropers,
    moltemplate will search for patterns contained in these strings to decide
    which type of interaction to generate.
    This function returns a dictionary that converts the original atom type name
    into these strings.
    """
    for line in lines_equivalences:
        #tokens = SplitQuotedString(line.strip(),
        #                           comment_char='!>')

        # skip past both '!' and '>' characters
        ic1 = line.find('!')
        ic = ic1
        ic2 = line.find('>')
        if ic2 != -1 and ic2 < ic1:
            ic = ic2
        if ic != -1:
            line = line[:ic]
        else:
            line = line.rstrip('\n')
        tokens = line.strip().split()
        #sys.stderr.write('DEBUG Equivalences2ffids():\n'
        #                 '      tokens = '+str(tokens)+'\n')
        atype = EncodeAName(tokens[2])
        atom2equiv_pair[atype] = EncodeAName(tokens[3])
        atom2equiv_bond[atype] = EncodeAName(tokens[4])
        atom2equiv_angle[atype] = EncodeAName(tokens[5])
        atom2equiv_dihedral[atype] = EncodeAName(tokens[6])
        atom2equiv_improper[atype] = EncodeAName(tokens[7])

    atom2ffid = OrderedDict()
    for atom in atom_types:
        atom2ffid[atom] = (atom + 
                           ',p'+atom2equiv_pair.get(atom,'') + 
                           ',b'+atom2equiv_bond.get(atom,'') + 
                           ',a'+atom2equiv_angle.get(atom,'') + 
                           ',d'+atom2equiv_dihedral.get(atom,'') + 
                           ',i'+atom2equiv_improper.get(atom,''))
    return atom2ffid






def AutoEquivalences2ffids(lines_equivalences,
                           lines_auto_equivalences,
                           atom_types,
                           atom2equiv_pair,
                           atom2equiv_bond,
                           atom2equiv_angle,
                           atom2equiv_dihedral,
                           atom2equiv_improper,
                           atom2auto_pair,
                           atom2auto_bondincr,
                           atom2auto_bond,
                           atom2auto_angleend,
                           atom2auto_anglecenter,
                           atom2auto_dihedralend,
                           atom2auto_dihedralcenter,
                           atom2auto_improperend,
                           atom2auto_impropercenter):
    """
    This function is a variant of Equivalences2ffids() which also considers
    "auto_equivalences".
    This function returns a dictionary that converts the original atom type name
    into a string that includes that atom's "equivalences",
    as well as its "auto_equivalences".
    moltemplate will search for patterns contained in these strings to decide
    which type of interaction to generate.
    """
    Equivalences2ffids(lines_equivalences,
                       atom_types,
                       atom2equiv_pair,
                       atom2equiv_bond,
                       atom2equiv_angle,
                       atom2equiv_dihedral,
                       atom2equiv_improper)

    # ------ The following lines are for processing "auto_equivalences" -----
    #
    # What is the difference between "equivalences" and "auto_equivalences"?
    #
    # equivalences:
    # Here is an excerpt from the Discover manual describing "equivalences":
    #  "Chemically distinct atoms often differ in some, but not all,
    #   of their forcefield parameters. For example, the bond parameters
    #  for the C-C bonds in ethene and in benzene are quite different,
    #  but the nonbond parameters for the carbon atoms are essentially
    #  the same. Rather than duplicating the nonbond parameters in the
    #  forcefield parameter file, the Discover program uses atom type
    #  equivalences to simplify the problem. In the example, the phenyl
    #  carbon atom type is equivalent to the pure sp2 carbons of ethene
    #  insofar as the nonbond parameters are concerned. The Discover
    #  program recognizes five types of equivalences for each atom
    #  type: nonbond, bond, angle, torsion, and out-of-plane.
    #  Cross terms such as bond-bond terms have the same equivalences
    #  (insofar as atom types are concerned) as the diagonal term of
    #  the topology of all the atoms defining the internal coordinates.
    #  For the bond-bond term, this means that the atom type
    #  equivalences for angles would be used
    #
    # auto_equivalences:
    #   Are similar to equivalences, but apparently with lower priority.
    #   In addition, it seems that, when looking up some of the class2 terms
    #   in the interaction according to atom type using "auto_equivalences"
    #   a distinction is made between end atoms and central atoms.
    #   The parameters for these interactions are also stored in different 
    #   tables in the .frc file, with different comments/tags.
    #   (for example, "cff91_auto" as opposed to "cff91")
    # An excerpt from the Discover manual is somewhat vague:
    #  "A forcefield may include automatic parameters for use when
    #   better-quality explicit parameters are not defined for a
    #   particular bond, angle, torsion, or out-of-plane interaction.
    #   These parameters are intended as temporary patches, to allow
    #   you to begin calculations immediately."

    for line in lines_auto_equivalences:
        #tokens = SplitQuotedString(line.strip(),
        #                           comment_char='!>')

        # skip past both '!' and '>' characters
        ic1 = line.find('!')
        ic = ic1
        ic2 = line.find('>')
        if ic2 != -1 and ic2 < ic1:
            ic = ic2
        if ic != -1:
            line = line[:ic]
        else:
            line = line.rstrip('\n')
        tokens = line.strip().split()
        #sys.stderr.write('DEBUG Equivalences2ffids():\n'
        #                 '      tokens = '+str(tokens)+'\n')
        atype = EncodeAName(tokens[2])
        atom2auto_pair[atype] = EncodeAName(tokens[3])
        atom2auto_bondincr[atype] = EncodeAName(tokens[4])
        atom2auto_bond[atype] = EncodeAName(tokens[5])
        atom2auto_angleend[atype] = EncodeAName(tokens[6])
        atom2auto_anglecenter[atype] = EncodeAName(tokens[7])
        atom2auto_dihedralend[atype] = EncodeAName(tokens[8])
        atom2auto_dihedralcenter[atype] = EncodeAName(tokens[9])
        atom2auto_improperend[atype] = EncodeAName(tokens[10])
        atom2auto_impropercenter[atype] = EncodeAName(tokens[11])

    atom2ffid = OrderedDict()
    for atom in atom_types:
        atom2ffid[atom] = (atom + 
                           ',p'+atom2equiv_pair.get(atom,'') + 
                           ',b'+atom2equiv_bond.get(atom,'') + 
                           ',a'+atom2equiv_angle.get(atom,'') + 
                           ',d'+atom2equiv_dihedral.get(atom,'') + 
                           ',i'+atom2equiv_improper.get(atom,'') + 
                           ',ap'+atom2auto_pair.get(atom,'') +
                           ',aq'+atom2auto_bondincr.get(atom,'') +
                           ',ab'+atom2auto_bond.get(atom,'') +
                           ',aae'+atom2auto_angleend.get(atom,'') + 
                           ',aac'+atom2auto_anglecenter.get(atom,'') + 
                           ',ade'+atom2auto_dihedralend.get(atom,'') + 
                           ',adc'+atom2auto_dihedralcenter.get(atom,'') + 
                           ',aie'+atom2auto_improperend.get(atom,'') + 
                           ',aic'+atom2auto_impropercenter.get(atom,'') +
                           ''
                          )
    return atom2ffid






def main():
    try:
        sys.stderr.write(g_program_name + ", version " +
                         __version__ + ", " + __date__ + "\n")
        if sys.version < '2.6':
            raise InputError('Error: Using python ' + sys.version + '\n' +
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
        ffname = 'BIOSYM_MSI_FORCE_FIELD'
        type_subset = set([])
        filename_in = ''
        file_in = sys.stdin
        #file_in = open('pcff_repaired.frc','r')  #CONTINUEHERE
        include_auto_equivalences = False
        #pair_style_name = 'lj/class2/coul/long'
        #pair_style_params = "10.0 10.0"
        pair_style2docs = {}
        pair_style2args = defaultdict(str)
        pair_style2docs['lj/cut/coul/long'] = 'http://lammps.sandia.gov/doc/pair_lj.html'
        pair_style2args['lj/cut/coul/long'] = '10.0'
        pair_style2docs['lj/class2/coul/long'] = 'http://lammps.sandia.gov/doc/pair_class2.html'
        pair_style2args['lj/class2/coul/long'] = '10.0'
        pair_style2docs['lj/class2/coul/cut'] = 'http://lammps.sandia.gov/doc/pair_class2.html'
        pair_style2args['lj/class2/coul/cut'] = '10.0'

        bond_style2docs = {}
        #bond_style2args = defaultdict(str)
        bond_style2docs['harmonic'] = 'http://lammps.sandia.gov/doc/bond_harmonic.html'
        bond_style2docs['class2'] = 'http://lammps.sandia.gov/doc/bond_class2.html'
        bond_style2docs['morse'] = 'http://lammps.sandia.gov/doc/bond_morse.html'
        bond_symmetry_subgraph = ''   # default

        angle_style2docs = {}
        #angle_style2args = defaultdict(str)
        angle_style2docs['harmonic'] = 'http://lammps.sandia.gov/doc/angle_harmonic.html'
        angle_style2docs['class2'] = 'http://lammps.sandia.gov/doc/angle_class2.html'
        angle_symmetry_subgraph = ''  # default

        dihedral_style2docs = {}
        #dihedral_style2args = defaultdict(str)
        dihedral_style2docs['charmm'] = 'http://lammps.sandia.gov/doc/dihedral_charmm.html'
        dihedral_style2docs['class2'] = 'http://lammps.sandia.gov/doc/dihedral_class2.html'
        dihedral_symmetry_subgraph = ''  # default

        improper_style2docs = {}
        #improper_style2args = defaultdict(str)
        improper_style2docs['cvff'] = 'http://lammps.sandia.gov/doc/improper_cvff.html'
        improper_style2docs['class2'] = 'http://lammps.sandia.gov/doc/improper_class2.html'
        improper_symmetry_subgraph = {}  #'cenJsortIKL'

        pair_mixing_style = 'sixthpower tail yes'

        special_bonds_command = 'special_bonds lj/coul 0.0 0.0 1.0 dihedral yes'
        # Thanks to Paul Saxe for is suggestions
        # http://lammps.sandia.gov/threads/msg11270.html


        kspace_style = 'kspace_style pppm 0.0001'
        pair_styles_selected = set([])
        #pair_style_link = 'http://lammps.sandia.gov/doc/pair_class2.html'
        pair_style_args = {}
        pair_cutoff = '10.0'
        #pair_style_command = "    pair_style hybrid " + \
        #    pair_style_name + " " + pair_style_args + "\n"
        bond_styles_selected = set([])
        #bond_style_link = bond_style2docs[bond_style_name]
        #bond_style_args = ''
        angle_styles_selected = set([])
        #angle_style_link = angle_style2docs[angle_style_name]
        #angle_style_args = ''
        dihedral_styles_selected = set([])
        #dihedral_style_link = dihedral_style2docs[dihedral_style_name]
        #dihedral_style_args = ''
        improper_styles_selected = set([])
        #improper_style_link = improper_style2docs[improper_style_name]
        #improper_style_args = ''
        hbond_style_name = ''
        hbond_style_link = ''
        hbond_style_args = ''

        lines_templates = []
        lines_references = defaultdict(list)
        lines_warnings = []
        
    
        argv = [arg for arg in sys.argv]
    
        i = 1
    
        while i < len(argv):
    
            #sys.stderr.write('argv['+str(i)+'] = \"'+argv[i]+'\"\n')
    
            if argv[i] == '-atoms':
                if i + 1 >= len(argv):
                    raise InputError('Error: the \"' + argv[i] + '\" argument should be followed by a quoted string\n'
                                     '       which contains a space-delimited list of of a subset of atom types\n'
                                     '       you want to use from the original force-field.\n'
                                     '       Make sure you enclose the entire list in quotes.\n')
                type_subset = set(argv[i + 1].strip('\"\'').strip().split())
                del argv[i:i + 2]
    
            elif argv[i] == '-name':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of the force-field\n')
                ffname = argv[i + 1]
                del argv[i:i + 2]
    
            elif argv[i] in ('-file', '-in-file'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of a force-field file\n')
                filename_in = argv[i + 1]
                try:
                    file_in = open(filename_in, 'r')
                except IOError:
                    sys.stderr.write('Error: Unable to open file\n'
                                     '       \"' + filename_in + '\"\n'
                                     '       for reading.\n')
                    sys.exit(1)
                del argv[i:i + 2]
    
            elif argv[i] == '-pair-cutoff':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a number'
                                     '       (the distance cutoff for non-bonded (pair) interactions)\n')
                pair_style_cutoff = argv[i+1]
                del argv[i:i + 2]

            elif argv[i] == '-pair-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by either \"lj/class2/coul/cut\" or \"lj/class2/coul/long\"\n')
                pair_style_list = argv[i + 1].split(',')
                for pair_style in pair_style_list:
                    if pair_style == '9-6':
                        pair_style = 'lj/class2/coul/long'
                    elif pair_style in ('12-6', 'lj', 'LJ'):
                        pair_style = 'lj/cut/coul/long'

                    if  pair_style.find('lj/class2/coul/long') == 0:
                        kspace_style = 'kspace_style pppm 0.0001'
                    elif pair_style.find('lj/cut/coul/long') == 0:
                        kspace_style = 'kspace_style pppm 0.0001'
                    elif pair_style.find('lj/class2/coul/cut') == 0:
                        pass
                        #kspace_style = ''
                    elif pair_style.find('lj/cut') == 0:
                        pass
                        #kspace_style = ''
                    else:
                        raise InputError('Error: ' + argv[i] + ' ' + pair_style + ' not supported.\n'
                                         '          The following pair_styles are supported:\n'
                                         '       lj/class2/coul/cut\n'
                                         '       lj/class2/coul/long\n'
                                         '       lj/cut\n'
                                         '       lj/cut/coul/long\n')
                    pair_styles_selected.add(pair_style)

                del argv[i:i + 2]

            elif argv[i] == '-bond-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible bond_style.\n')
                bond_styles = argv[i + 1].split(',')
                for bond_style in bond_styles:
                    bond_styles_selected.add(bond_style)
                #bond_style2args[bond_style] = argv[i + 1].split()[1:]
                #if bond_style_name.find('harmonic') == 0:
                #    pass
                #    #bond_style_link = 'http://lammps.sandia.gov/doc/bond_harmonic.html'
                #elif bond_style_name.find('morse') == 0:
                #    pass
                #    #bond_style_link = 'http://lammps.sandia.gov/doc/bond_morse.html'
                #elif bond_style_name.find('class2') == 0:
                #    pass
                #    #bond_style_link = 'http://lammps.sandia.gov/doc/bond_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\", \"class2\", or \"morse\".\n')
                del argv[i:i + 2]

            elif argv[i] == '-angle-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible angle_style.\n')
                angle_styles = argv[i + 1].split(',')
                for angle_style in angle_styles:
                    angle_styles_selected.add(angle_style)
                #if angle_style_name.find('harmonic') == 0:
                #    pass
                #    #angle_style_link = 'http://lammps.sandia.gov/doc/angle_harmonic.html'
                #elif angle_style_name.find('class2') == 0:
                #    pass
                #    #angle_style_link = 'http://lammps.sandia.gov/doc/angle_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\" or \"class2\"\n')
                del argv[i:i + 2]

            elif argv[i] == '-dihedral-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible dihedral_style.\n')
                dihedral_styles = argv[i + 1].split(',')
                for dihedral_style in dihedral_styles:
                    dihedral_styles_selected.add(dihedral_style)
                #if dihedral_style_name.find('charmm') == 0:
                #    pass
                #    #dihedral_style_link = 'http://lammps.sandia.gov/doc/dihedral_charmm.html'
                #elif dihedral_style_name.find('class2') == 0:
                #    pass
                #    #dihedral_style_link = 'http://lammps.sandia.gov/doc/dihedral_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\" or \"class2\"\n')
                del argv[i:i + 2]

            elif argv[i] == '-improper-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       a compatible impropoer_style.\n')
                improper_styles = argv[i + 1].split(',')
                for improper_style in improper_styles:
                    improper_styles_selected.add(improper_style)
                #if impropoer_style_name.find('harmonic') == 0:
                #    pass
                #    #impropoer_style_link = 'http://lammps.sandia.gov/doc/impropoer_harmonic.html'
                #elif impropoer_style_name.find('class2') == 0:
                #    pass
                #    #impropoer_style_link = 'http://lammps.sandia.gov/doc/impropoer_class2.html'
                #else:
                #    raise InputError('Error: ' + argv[i] + ' must be followed by either:\n'
                #                     '       \"harmonic\" or \"class2\"\n')
                del argv[i:i + 2]

            elif argv[i] == '-hbond-style':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' ' + hbond_style_name + '\n'
                                     '       should be followed by a compatible pair_style.\n')
                hbond_style_name = argv[i + 1]
                hbond_style_link = 'http://lammps.sandia.gov/doc/pair_hbond_dreiding.html'
                if hbond_style_name.find('none') == 0:
                    hbond_style_name = ''
                    hbond_style_args = ''
                elif hbond_style_name.find('hbond/dreiding/lj') == 0:
                    n = len('hbond/dreiding/lj')
                    hbond_style_args = hbond_style_name[n+1:]
                    hbond_style_name = hbond_style_name[:n]
                elif hbond_style_name.find('hbond/dreiding/morse') == 0:
                    n = len('hbond/dreiding/morse')
                    hbond_style_args = hbond_style_name[n+1:]
                    hbond_style_name = hbond_style_name[:n]
                else:
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by either\n'
                                     '       \"hbond/dreiding/lj\" or \"hbond/dreiding/morse"\n')
                del argv[i:i + 2]

            elif argv[i] in ('-url', '-in-url'):
                import urllib2
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a URL pointing to\n'
                                     '       a file containing force-field information in msi/frc format.\n')
                url = argv[i + 1]
                try:
                    request = urllib2.Request(url)
                    file_in = urllib2.urlopen(request)
                except urllib2.URLError:
                    sys.stdout.write("Error: Unable to open link:\n" + url + "\n")
                    sys.exit(1)
                del argv[i:i + 2]
    
            elif argv[i] == '-auto':
                include_auto_equivalences = True
                del argv[i:i + 1]
    
            elif argv[i] in ('-help', '--help', '-?', '--?'):
                sys.stderr.write(doc_msg)
                sys.exit(0)
                del argv[i:i + 1]
    
            else:
                i += 1
    
        if len(argv) != 1:
            raise InputError('Error: Unrecongized arguments: ' + ' '.join(argv[1:]) +
                             '\n\n' + doc_msg)

        # Default styles:
        if len(bond_styles_selected) == 0:
            bond_styles_selected.add('class2')
        if len(angle_styles_selected) == 0:
            angle_styles_selected.add('class2')
        if len(dihedral_styles_selected) == 0:
            dihedral_styles_selected.add('class2')
        if len(improper_styles_selected) == 0:
            improper_styles_selected.add('class2')
        if len(pair_styles_selected) == 0:
            pair_styles_selected.add('lj/class2/coul/long')

        #sys.stderr.write("Reading parameter file...\n")

        lines = file_in.readlines()
        atom2charge = OrderedDict()  # lookup charge from atom type
        atom2mass = OrderedDict()  # lookup mass from atom type
        # equivalences lookup
        atom2ffid = OrderedDict()  # lookup "force-field-ID" a string containing
                                   # equivalences to lookup bonded interactions
        atom2equiv_pair = OrderedDict() # lookup the equivalent symbol used for
                                        # looking up pair interactions
        atom2equiv_bond = OrderedDict()
        atom2equiv_angle = OrderedDict()
        atom2equiv_dihedral = OrderedDict()
        atom2equiv_improper = OrderedDict()
        # inverse equivalences lookup
        equiv_pair2atom = defaultdict(set)
        equiv_bond2atom = defaultdict(set)
        equiv_angle2atom = defaultdict(set)
        equiv_dihedral2atom = defaultdict(set)
        equiv_improper2atom = defaultdict(set)

        # auto equivalences lookup
        atom2auto_pair = OrderedDict()
        atom2auto_bondincr = OrderedDict()
        atom2auto_bond = OrderedDict()
        atom2auto_angleend = OrderedDict()
        atom2auto_anglecenter = OrderedDict()
        atom2auto_dihedralend = OrderedDict()
        atom2auto_dihedralcenter = OrderedDict()
        atom2auto_improperend = OrderedDict()
        atom2auto_impropercenter = OrderedDict()
        # inverse auto equivalences lookup
        auto_pair2atom = defaultdict(set)
        auto_bondincr2atom = defaultdict(set)
        auto_bond2atom = defaultdict(set)
        auto_angleend2atom = defaultdict(set)
        auto_anglecenter2atom = defaultdict(set)
        auto_dihedralend2atom = defaultdict(set)
        auto_dihedralcenter2atom = defaultdict(set)
        auto_improperend2atom = defaultdict(set)
        auto_impropercenter2atom = defaultdict(set)


        atom2element = OrderedDict()  # Optional:
                                      # which element (eg 'C', 'O') ? (Note this
                                      # is different from atom type: 'C1', 'Oh')
        atom2numbonds = OrderedDict() # Optional: how many bonds emanate from
        atom2descr = OrderedDict()    # Optional: a brief description
        atom2ver = OrderedDict()  # atoms introduced in different versions of ff
        atom2ref = OrderedDict()  # reference to paper where atom introduced
        lines_equivalences = []      # equivalences for force-field lookup
        lines_auto_equivalences = [] # auto_equivalences have lower priority

        pair2params = OrderedDict()
        pair2style = OrderedDict()
        pair_styles = set([])
        pair2ver = OrderedDict()
        pair2ref = OrderedDict()

        bond2chargepair = OrderedDict()      # a.k.a "bond increments"
        charge_pair_priority = OrderedDict() # priority in case multiple entries
                                             # exist for the same pair of atoms
        charge_pair_ver = OrderedDict()      # which version of the force field?
        charge_pair_ref = OrderedDict()      # paper introducing this chargepair

        bond2params = OrderedDict()  # store a tuple with the 2-body bond
                                     # interaction type, and its parameters
                                     # for every type of bond
        bond2priority = OrderedDict() # What is the priority of this interaction?
        bond2style = OrderedDict()    # What LAMMPS bond style (formula)
                                      # is used for a given interaction?
        bond_styles = set([])         # Contains all bond styles used.
        bond2ver = OrderedDict()
        bond2ref = OrderedDict()
        bond2r0 = OrderedDict()
        bond2r0_auto = OrderedDict()

        angle2params = OrderedDict() # store a tuple with the 3-body angle
                                     # interaction type, and its parameters
                                     # for every type of angle

        angle2params_or = OrderedDict()
        # http://lammps.sandia.gov/doc/angle_class2.html
        #angle2class2_a = OrderedDict()  # params for the "a" class2 terms
        angle2class2_bb = OrderedDict() # params for the "bb" class2 terms
        angle2class2_bb_or = OrderedDict()
        angle2class2_ba = OrderedDict() # params for the "ba" class2 terms
        angle2class2_ba_or = OrderedDict()
        angle2priority = OrderedDict()  # What is the priority of this interaction?
        angle2priority_or = OrderedDict()
        angle_is_secondary_or = OrderedDict()
        angle2style = OrderedDict()    # What LAMMPS angle style (formula)
                                       # is used for a given interaction?
        angle2style_or = OrderedDict()
        angle_styles = set([])         # Contains all angle styles used.
        angle2ref = OrderedDict()
        angle2ver = OrderedDict()
        angle2ref_or = OrderedDict()
        angle2ver_or = OrderedDict()
        angle2ver_bb = OrderedDict()
        angle2ver_bb_or = OrderedDict()
        angle2ref_bb = OrderedDict()
        angle2ref_bb_or = OrderedDict()
        angle2ver_ba = OrderedDict()
        angle2ver_ba_or = OrderedDict()
        angle2ref_ba = OrderedDict()
        angle2ref_ba_or = OrderedDict()
        angle2theta0_or = OrderedDict()
        angle2theta0_auto_or = OrderedDict()

        # http://lammps.sandia.gov/doc/dihedral_class2.html
        dihedral2params = OrderedDict() # store a tuple with the 4-body dihedral
                                        # interaction type, and its parameters
                                        # for every type of dihedral
        dihedral2params_or = OrderedDict()
        #dihedral2class2_d = OrderedDict() # params for the "d" class2 term
        dihedral2class2_mbt = OrderedDict() # params for the "mbt" class2 term
        dihedral2class2_mbt_or = OrderedDict()
        dihedral2class2_ebt = OrderedDict() # params for the "ebt" class2 term
        dihedral2class2_ebt_or = OrderedDict()
        #dihedral2sym_ebt = OrderedDict()
        dihedral2class2_at = OrderedDict() # params for the "at" class2 term
        dihedral2class2_at_or = OrderedDict()
        #dihedral2sym_at = OrderedDict()
        dihedral2class2_aat = OrderedDict() # params for the "aat" class2 term
        dihedral2class2_aat_or = OrderedDict()
        #dihedral2sym_aat = OrderedDict()
        dihedral2class2_bb13 = OrderedDict() # params for the "bb13" class2 term
        dihedral2class2_bb13_or = OrderedDict()
        #dihedral2sym_bb13 = OrderedDict()
        dihedral2priority = OrderedDict()  # What is the priority of this interaction?
        dihedral2priority_or = OrderedDict()
        dihedral_is_secondary_or = OrderedDict()
        dihedral2style = OrderedDict()    # What LAMMPS dihedral style (formula)
                                          # is used for a given interaction?
        dihedral2style_or = OrderedDict()
        dihedral_styles = set([])         # Contains all dihedral styles used.
        dihedral2ref = OrderedDict()
        dihedral2ver = OrderedDict()
        dihedral2ver_or = OrderedDict()
        dihedral2ref_or = OrderedDict()
        dihedral2ver_mbt = OrderedDict()
        dihedral2ver_mbt_or = OrderedDict()
        dihedral2ref_mbt = OrderedDict()
        dihedral2ref_mbt_or = OrderedDict()
        dihedral2ver_ebt = OrderedDict()
        dihedral2ver_ebt_or = OrderedDict()
        dihedral2ref_ebt = OrderedDict()
        dihedral2ref_ebt_or = OrderedDict()
        dihedral2ver_at = OrderedDict()
        dihedral2ver_at_or = OrderedDict()
        dihedral2ref_at = OrderedDict()
        dihedral2ref_at_or = OrderedDict()
        dihedral2ver_aat = OrderedDict()
        dihedral2ver_aat_or = OrderedDict()
        dihedral2ref_aat = OrderedDict()
        dihedral2ref_aat_or = OrderedDict()
        dihedral2ver_bb13 = OrderedDict()
        dihedral2ver_bb13_or = OrderedDict()
        dihedral2ref_bb13 = OrderedDict()
        dihedral2ref_bb13_or = OrderedDict()


        # http://lammps.sandia.gov/doc/improper_class2.html
        improper2params = OrderedDict() # store a tuple with the 4-body improper
                                        # interaction type, and its parameters
                                        # for every type of imporpoer
        improper2params_or = OrderedDict()
        improper2class2_aa = OrderedDict() # params for the "aa" class2 term
        improper2class2_aa_or = OrderedDict()

        improper2cross = defaultdict(dict)
                           # improper2cross[imp_name][atoms] stores the 
                           # coefficient (K) for the angle-angle ("aa") 
                           # improper interactions between a pair of 
                           # neighboring 3-body angles (in the .FRC file).
                           # "imp_name" is the name of the improper interaction
                           #   (which is a concatination of the central atom and
                           #   the 3 surrounding leaf atoms (which are sorted))
                           # "atoms" indicates, for that K value, the list of
                           #   leaf atoms for that K value as they appear in the
                           #   corresponding line of the .frc file (however the
                           #   and last atom names are swapped if the first
                           #   atom name is lexicographically > the last, to
                           #   eliminate redundancy and ambiguity.)

        improper2sym = defaultdict(set)
                           # improper2sym[imp_name] indicates which subset of
                           # leaf atoms (from 0 to 2) are equivalent and can
                           # tolerate having their order rearranged without
                           # effecting the energy.  Later on this will be used
                           # to reduce the number of improper interactions that
                           # will be generated by moltemplate.

        improper2priority = OrderedDict() # What is the priority of this interaction?
        improper2priority_or = OrderedDict()
        improper_is_secondary_or = OrderedDict()
        improper2style = OrderedDict()    # What LAMMPS improper style (formula)
                                          # is used for a given interaction?
        improper2style_or = OrderedDict()
        improper_styles = set([])         # Contains all improper styles used.
        improper2ver = OrderedDict()
        improper2ver_or = OrderedDict()
        improper2ref = OrderedDict()
        improper2ref_or = OrderedDict()
        improper2ver_aa = OrderedDict()
        improper2ver_aa_or = OrderedDict()
        improper2ref_aa = OrderedDict()
        improper2ref_aa_or = OrderedDict()


        # Warn users if force field contains terms which cannot yet
        # be simulated with LAMMPS (as of 2017-10-13)
        display_OOP_OOP_warning = False
        display_torsion_torsion_1_warning = False


        """
         --- these next few lines of code appear to be unnecessary.
         --- I'll probably delete this code in a later version
        hbond2params = OrderedDict()    # lookup hbond parameters and atom types
        hbond2donors = OrderedDict()    # according to the identifier in the 2nd
        hbond2acceptors = OrderedDict() #  column of the "#hbond_definition"
        hbond2hydrogens = OrderedDict() # section of an .frc file.
        """

        allowed_section_names = set(['#define',
                                     # sections used in all MSI force-fields
                                     '#atom_types',
                                     '#equivalence',
                                     '#auto_equivalence',
                                     '#nonbond(9-6)',
                                     '#nonbond(12-6)',
                                     '#quadratic_bond',
                                     '#quartic_bond',
                                     '#morse_bond',
                                     '#quadratic_angle',
                                     '#quartic_angle',
                                     '#bond-bond',
                                     '#bond-angle',
                                     '#torsion_1',
                                     '#torsion_3',
                                     '#middle_bond-torsion_3',
                                     '#end_bond-torsion_3',
                                     '#angle-torsion_3',
                                     '#angle-angle-torsion_1',#(class2 dihedral)
                                     '#bond-bond_1_3', #(a class2 dihedral term)
                                     '#out_of_plane',
                                     '#wilson_out_of_plane',
                                     '#angle-angle',   #(a class2 improper term)
                                     '#out_of_plane-out_of_plane', # UNSUPPORTED
                                     '#torsion-torsion_1',         # UNSUPPORTED
                                     '#bond_increments',
                                     '#hbond_definition',          # irrelevant?
                                     '#templates',
                                     '#reference',
                                     '#end'
                                     ])

        icol_type = icol_mass = icol_elem = icol_nbonds = icol_comment = icol_ver = icol_ref = -1

        section_name = ''
        section_is_auto = False

        sys.stderr.write("parsing file pass1: look for atom types and equivalences...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue
            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if tokens[0] in allowed_section_names:
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    tokens_after_section_name = tokens[1:]
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif not tokens[0] in ('#version',
                                       '#define'):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')
            elif (len(tokens) == 8) and (section_name == '#equivalence'):
                if line.lstrip().find('!') == 0:
                    continue
                lines_equivalences.append(line)
            elif (len(tokens) == 12) and (section_name == '#auto_equivalence'):
                if line.lstrip().find('!') == 0:
                    continue
                lines_auto_equivalences.append(line)
            elif (len(tokens) > 0) and (section_name == '#atom_types'):
                # Different FRC files put this information in different
                # columns.  Column order is stored in the !Ver comment line:
                if line.lstrip().find('!Ver') == 0:
                    tokens = line.strip().split()
                    for i in range(0, len(tokens)):
                        if tokens[i].lower() == 'type':
                            icol_type = i
                        elif tokens[i].lower() == 'mass':
                            icol_mass = i
                        elif tokens[i].lower() == 'element':
                            icol_elem = i
                        elif tokens[i].lower() == 'connections':
                            icol_nbonds = i
                        elif tokens[i].lower() == 'comment':
                            icol_comment = i
                        elif tokens[i].lower() == '!ver':   #(version of ff)
                            icol_ver = i
                        elif tokens[i].lower() == 'ref':
                            icol_ref = i
                    assert(icol_ver == 0)

                    if -1 in (icol_type, icol_mass):
                        raise InputError('Error: Invalid #atom_types section.\n'
                                         '       The meaning of each column cannot be determined.\n'
                                         '       This file needs a valid "!Ver..." comment.\n')
                    if icol_comment == -1:
                        icol_comment = max(icol_type, icol_mass,
                                           icol_elem, icol_nbonds) + 1

                    sys.stderr.write('icol_ver = '+str(icol_ver)+'\n')
                    sys.stderr.write('icol_ref = '+str(icol_ref)+'\n')
                    sys.stderr.write('icol_mass = '+str(icol_mass)+'\n')
                    sys.stderr.write('icol_nelem = '+str(icol_elem)+'\n')
                    sys.stderr.write('icol_nbonds = '+str(icol_nbonds)+'\n')
                    sys.stderr.write('icol_comment = '+str(icol_comment)+'\n')
                    continue

                tokens = map(RemoveOuterQuotes,
                             NSplitQuotedString(line.strip(),
                                                icol_comment+1,
                                                quotes='',
                                                comment_char='>'))
                tokens = list(tokens)

                if (len(tokens) > 4):
                    if ((len(type_subset) == 0) or (tokens[1] in type_subset)):
                        aname = EncodeAName(tokens[icol_type])
                        atom2mass[aname] = str(max(float(tokens[icol_mass]), 1.0e-06))
                        # Some atoms in cvff.prm have zero mass. Unfortunately this
                        # causes LAMMPS to crash, even if these atoms are never used,
                        # so I give the mass a non-zero value instead.

                        if icol_elem != -1:
                            atom2element[aname] = tokens[icol_elem]
                        if icol_nbonds != -1:
                            atom2numbonds[aname] = int(tokens[icol_nbonds])
                        atom2descr[aname] = tokens[icol_comment]
                        atom2ver[aname] = tokens[icol_ver]
                        atom2ref[aname] = tokens[icol_ref]

                elif len(tokens) > 0:
                    raise InputError('Error: Invalid atom line: (line#'+str(iline)+')\n' +
                                     '\"'+line.strip()+'\"')

        atom_types = [x for x in atom2mass]

        # Now construct the lookup tables and inverse tables
        # we will need to understand the remainder of the file:
        if not include_auto_equivalences:
            atom2ffid = Equivalences2ffids(lines_equivalences,
                                           atom_types,
                                           atom2equiv_pair,
                                           atom2equiv_bond,
                                           atom2equiv_angle,
                                           atom2equiv_dihedral,
                                           atom2equiv_improper)
        else:
            atom2ffid = AutoEquivalences2ffids(lines_equivalences,
                                               lines_auto_equivalences,
                                               atom_types,
                                               atom2equiv_pair,
                                               atom2equiv_bond,
                                               atom2equiv_angle,
                                               atom2equiv_dihedral,
                                               atom2equiv_improper,
                                               atom2auto_pair,
                                               atom2auto_bondincr,
                                               atom2auto_bond,
                                               atom2auto_angleend,
                                               atom2auto_anglecenter,
                                               atom2auto_dihedralend,
                                               atom2auto_dihedralcenter,
                                               atom2auto_improperend,
                                               atom2auto_impropercenter)

        for a,e in atom2equiv_pair.items():
            equiv_pair2atom[e].add(a)
        for a,e in atom2equiv_bond.items():
            equiv_bond2atom[e].add(a)
        for a,e in atom2equiv_angle.items():
            equiv_angle2atom[e].add(a)
        for a,e in atom2equiv_dihedral.items():
            equiv_dihedral2atom[e].add(a)
        for a,e in atom2equiv_improper.items():
            equiv_improper2atom[e].add(a)

        # the inverse lookup for '*' matches all atom types
        for a in atom_types:
            #equiv_pair2atom['*'].add(EncodeAName(a))
            equiv_pair2atom['X'].add(EncodeAName(a))
            #equiv_bond2atom['*'].add(EncodeAName(a))
            equiv_bond2atom['X'].add(EncodeAName(a))
            #equiv_angle2atom['*'].add(EncodeAName(a))
            equiv_angle2atom['X'].add(EncodeAName(a))
            #equiv_dihedral2atom['*'].add(EncodeAName(a))
            equiv_dihedral2atom['X'].add(EncodeAName(a))
            #equiv_improper2atom['*'].add(EncodeAName(a))
            equiv_improper2atom['X'].add(EncodeAName(a))

        for a,e in atom2auto_pair.items():
            auto_pair2atom[e].add(a)
        for a,e in atom2auto_bondincr.items():
            auto_bondincr2atom[e].add(a)
        for a,e in atom2auto_bond.items():
            auto_bond2atom[e].add(a)
        for a,e in atom2auto_angleend.items():
            auto_angleend2atom[e].add(a)
            #auto_angle[0][e].add(a)
            #auto_angle[2][e].add(a)
        for a,e in atom2auto_anglecenter.items():
            auto_anglecenter2atom[e].add(a)
            #auto_angle[1][e].add(a)
        for a,e in atom2auto_dihedralend.items():
            auto_dihedralend2atom[e].add(a)
            #auto_dihedral2atom[0][e].add(a)
            #auto_dihedral2atom[3][e].add(a)
        for a,e in atom2auto_dihedralcenter.items():
            auto_dihedralcenter2atom[e].add(a)
            #auto_dihedral2atom[1][e].add(a)
            #auto_dihedral2atom[2][e].add(a)
        for a,e in atom2auto_improperend.items():
            auto_improperend2atom[e].add(a)
        for a,e in atom2auto_impropercenter.items():
            auto_impropercenter2atom[e].add(a)

        # the inverse lookup for '*' matches all atom types
        for a in atom_types:
            #auto_pair2atom['*'].add(EncodeAName(a))
            auto_pair2atom['X'].add(EncodeAName(a))
            #auto_bondincr2atom['*'].add(EncodeAName(a))
            auto_bondincr2atom['X'].add(EncodeAName(a))
            #auto_bond2atom['*'].add(EncodeAName(a))
            auto_bond2atom['X'].add(EncodeAName(a))
            #auto_angleend2atom['*'].add(EncodeAName(a))
            auto_angleend2atom['X'].add(EncodeAName(a))
            #auto_anglecenter2atom['*'].add(EncodeAName(a))
            auto_anglecenter2atom['X'].add(EncodeAName(a))
            #auto_dihedralend2atom['*'].add(EncodeAName(a))
            auto_dihedralend2atom['X'].add(EncodeAName(a))
            #auto_dihedralcenter2atom['*'].add(EncodeAName(a))
            auto_dihedralcenter2atom['X'].add(EncodeAName(a))
            #auto_improperend2atom['*'].add(EncodeAName(a))
            auto_improperend2atom['X'].add(EncodeAName(a))
            #auto_impropercenter2atom['*'].add(EncodeAName(a))
            auto_impropercenter2atom['X'].add(EncodeAName(a))








        sys.stderr.write("parsing file pass2: look for bonds, bond_increments and nonbonded (pair) interactions...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue
            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if (tokens[0] in allowed_section_names):
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    tokens_after_section_name = tokens[1:]
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif (not tokens[0] in ('#version','#define')):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')


            elif ((len(tokens) > 4) and (section_name == '#nonbond(12-6)')
                  and (pair_styles_selected & set(['lj','lj/cut','lj/cut/coul/long',
                                                   'lj/cut/coul/cut','lj/cut/coul/debye',
                                                   'lj/cut/coul/dsf','lj/cut/coul/msm',
                                                   '12-6','nonbond(12-6)']))):

                if line.lstrip().find('!') == 0:
                    continue
                atom_name = EncodeAName(tokens[2])
                pair2ver[atom_name] = tokens[0]
                pair2ref[atom_name] = tokens[1]
                A = float(tokens[3])
                B = float(tokens[4])
                epsilon = B*B/(4*A)
                sigma = pow(B/A, 1.0/6)
                if sigma == 0.0:
                    sigma = 1.0   #(non-zero to avoid nan error later)
                pair_styles.add('lj/cut/coul/long')
                pair_style_args['lj/cut/coul/long'] = pair_cutoff
                pair2style[atom_name] = 'lj/cut/coul/long'
                pair2params[atom_name] = (str(epsilon)+' '+str(sigma))
                pair_mixing_style = 'geometric tail yes'
                #if pair_style_name.find('lj/cut') == 0:
                #    pair2params[atom_name] = (str(epsilon)+' '+str(sigma))
                #    pair_mixing_style = 'geometric tail yes'


            elif ((len(tokens) > 4) and (section_name == '#nonbond(9-6)')
                  and (pair_styles_selected &
                       set(['class2', '9-6', 'nonbond(9-6)',
                            'lj/class2/coul/long']))):
                if line.lstrip().find('!') == 0:
                    continue
                atom_name = EncodeAName(tokens[2])
                pair2ver[atom_name] = tokens[0]
                pair2ref[atom_name] = tokens[1]
                sigma = tokens[3]
                epsilon = tokens[4]
                pair_styles.add('lj/class2/coul/long')
                pair_style_args['lj/class2/coul/long'] = pair_cutoff
                pair2style[atom_name] = 'lj/class2/coul/long'
                pair2params[atom_name] = (epsilon+' '+sigma)
                pair_mixing_style = 'sixthpower tail yes'
                #if pair_style_name.find('lj/class2') == 0:
                #    pair2params[atom_name] = (epsilon+' '+sigma)
                #    pair_mixing_style = 'sixthpower tail yes'


            elif (len(tokens) == 6) and (section_name == '#bond_increments'):
                if line.lstrip().find('!') == 0:
                    continue
                aorig = [a for a in map(EncodeAName, tokens[2:4])]
                delta_q = tokens[4:6]
                atom_names = [a for a in aorig]
                # swap the order of the atoms?
                order_reversed = aorig[0] > aorig[-1]
                if order_reversed:
                    delta_q.reverse()
                    atom_names.reverse()
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                charge_pair_ver[bond_name] = tokens[0]
                charge_pair_ref[bond_name] = tokens[1]
                charge_pair_priority[bond_name] = \
                    (section_is_auto,
                     DetermineNumericPriority(section_is_auto,
                                              tokens[2:4],
                                              float(charge_pair_ver[bond_name])))
                bond2chargepair[bond_name] = (delta_q[0] + ' ' + delta_q[1])


            elif ((len(tokens) > 5) and (section_name == '#quadratic_bond')
                  and (bond_styles_selected & set(['harmonic','quadratic','quadratic_bond']))):
                if line.lstrip().find('!') == 0:
                    continue
                bond_styles.add('harmonic')
                atom_names = SortByEnds(map(EncodeAName, tokens[2:4]))
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                bond2ver[bond_name] = tokens[0]
                bond2ref[bond_name] = tokens[1]
                bond2priority[bond_name] = \
                    (section_is_auto,
                     DetermineNumericPriority(section_is_auto,
                                              tokens[2:4],
                                              float(bond2ver[bond_name])))
                r0 = tokens[4]
                k = tokens[5]
                if not section_is_auto:
                    bond2r0[bond_name] = r0
                    sys.stderr.write('bond2r0['+bond_name+'] = ' + str(r0) + '\n')
                else:
                    bond2r0_auto[(atom_names[0], atom_names[1])] = r0
                    sys.stderr.write('bond2r0_auto['+str(atom_names)+'] = ' + str(r0) + '\n')
                bond2style[bond_name] = 'harmonic'
                bond2params[bond_name] = (k+' '+r0)


            elif ((len(tokens) > 6) and (section_name == '#morse_bond')
                  and (bond_styles_selected & set(['morse','morse_bond']))):
                if line.lstrip().find('!') == 0:
                    continue
                bond_styles.add('morse')
                atom_names = SortByEnds(map(EncodeAName, tokens[2:4]))
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                bond2ver[bond_name] = tokens[0]
                bond2ref[bond_name] = tokens[1]
                bond2priority[bond_name] = \
                    (section_is_auto,
                     DetermineNumericPriority(section_is_auto,
                                              tokens[2:4],
                                              float(bond2ver[bond_name])))
                r0 = tokens[4]
                D = tokens[5]
                alpha = tokens[6]
                sys.stderr.write('DEBUG: morse: atom_names = '+str(atom_names)+'\n')
                if not section_is_auto:
                    bond2r0[bond_name] = r0
                    sys.stderr.write('bond2r0['+bond_name+'] = ' + str(r0) + '\n')
                else:
                    bond2r0_auto[(atom_names[0], atom_names[1])] = r0
                    sys.stderr.write('bond2r0_auto['+str(atom_names)+'] = ' + str(r0) + '\n')
                bond2style[bond_name] = 'morse'
                bond2params[bond_name] = (D+' '+alpha+' '+r0)



            elif ((len(tokens) > 7) and (section_name == '#quartic_bond')
                  and (bond_styles_selected & set(['class2','quartic','quartic_bond']))):
                if line.lstrip().find('!') == 0:
                    continue
                bond_styles.add('class2')
                atom_names = SortByEnds(map(EncodeAName, tokens[2:4]))
                bond_name = EncodeInteractionName(atom_names, section_is_auto)
                bond2ver[bond_name] = tokens[0]
                bond2ref[bond_name] = tokens[1]
                bond2priority[bond_name] = \
                    (section_is_auto,
                     DetermineNumericPriority(section_is_auto,
                                              tokens[2:4],
                                              float(bond2ver[bond_name])))
                r0 = tokens[4]
                if not section_is_auto:
                    bond2r0[bond_name] = r0
                    sys.stderr.write('bond2r0['+bond_name+'] = ' + str(r0) + '\n')
                else:
                    bond2r0_auto[(atom_names[0], atom_names[1])] = r0
                    sys.stderr.write('bond2r0_auto['+str(atom_names)+'] = ' + str(r0) + '\n')
                K2 = tokens[5]
                K3 = tokens[6]
                K4 = tokens[7]
                bond2style[bond_name] = 'class2'
                bond2params[bond_name] = (r0+' '+K2+' '+K3+' '+K4)





        sys.stderr.write("parsing file pass3: look for (3-body) angle interactions...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue
            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if (tokens[0] in allowed_section_names):
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    tokens_after_section_name = tokens[1:]
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif (not tokens[0] in ('#version','#define')):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')







            elif (len(tokens) > 6) and (section_name == '#quadratic_angle'):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = SortByEnds(map(EncodeAName, tokens[2:5]))
                angle_name = EncodeInteractionName(atom_names, section_is_auto)

                angle2ver[angle_name] = tokens[0]
                angle2ref[angle_name] = tokens[1]
                angle2priority_or[angle_name] = \
                    DetermineNumericPriority(section_is_auto,
                                             tokens[2:5],
                                             float(angle2ver[angle_name]))
                angle_is_secondary_or[angle_name] = False
                angle2priority[angle_name] = \
                    (section_is_auto,
                     angle_is_secondary_or[angle_name],
                     angle2priority_or[angle_name])
                theta0 = tokens[5]
                k = tokens[6]
                if not section_is_auto:
                    angle2theta0_or[angle_name] = theta0
                    sys.stderr.write('angle2theta0_or['+angle_name+'] = ' + str(theta0) + '\n')
                else:
                    angle2theta0_auto_or[(atom_names[0], atom_names[1], atom_names[2])] = theta0
                    sys.stderr.write('angle2theta0_auto_or['+str(atom_names)+'] = ' + str(theta0) + '\n')
                if (angle_styles_selected & set(['harmonic',
                                                 'quadratic',
                                                 'quadratic_angle'])):
                    angle_styles.add('harmonic')
                    angle2style[angle_name] = 'harmonic'
                    angle2params[angle_name] = (k+' '+theta0)
                elif (angle_styles_selected & set(['class2',
                                                   'quartic',
                                                   'quartic_angle'])):
                    # Then this is a special case of the class2 angle where
                    # the (theta-theta0)^3 and (theta-theta0)^4 terms = 0
                    angle_styles.add('class2')
                    angle2style_or[angle_name] = 'class2'
                    angle2params_or[angle_name] = (theta0+' '+k+' 0 0')



            elif ((len(tokens) > 8) and (section_name == '#quartic_angle')
                  and (angle_styles_selected & set(['class2','quartic','quartic_angle']))):
                if line.lstrip().find('!') == 0:
                    continue
                angle_styles.add('class2')
                atom_names = SortByEnds(map(EncodeAName, tokens[2:5]))
                ang_name_orig = EncodeInteractionName(atom_names, section_is_auto)
                version = tokens[0]
                reference = tokens[1]
                angle2ver_or[ang_name_orig] = version
                angle2ref_or[ang_name_orig] = reference
                angle2priority_or[ang_name_orig] = \
                    DetermineNumericPriority(section_is_auto,
                                             tokens[2:5],
                                             float(angle2ver_or[ang_name_orig]))
                angle_is_secondary_or[ang_name_orig] = False
                #angle2priority[ang_name_orig] = \
                #    (section_is_auto,
                #     angle_is_secondary_or[ang_name_orig],
                #     angle2priority_or[ang_name_orig])
                theta0 = tokens[5]
                if not section_is_auto:
                    angle2theta0_or[ang_name_orig] = theta0
                    sys.stderr.write('angle2theta0_or['+ang_name_orig+'] = ' + str(theta0) + '\n')
                else:
                    angle2theta0_auto_or[(atom_names[0], atom_names[1], atom_names[2])] = theta0
                    sys.stderr.write('angle2theta0_auto_or['+str(atom_names)+'] = ' + str(theta0) + '\n')
                K2 = tokens[6]
                K3 = tokens[7]
                K4 = tokens[8]
                angle2style_or[ang_name_orig] = 'class2'
                angle2params_or[ang_name_orig] = [theta0, K2, K3, K4]
                if not ang_name_orig in angle2class2_bb_or:
                    angle2class2_bb_or[ang_name_orig] = '0.0'          # default value
                    angle2ver_bb_or[ang_name_orig] = version           # default value
                    angle2ref_bb_or[ang_name_orig] = reference         # default value
                if not ang_name_orig in angle2class2_ba_or:
                    angle2class2_ba_or[ang_name_orig] = ['0.0', '0.0'] # default value
                    angle2ver_ba_or[ang_name_orig] = version           # default value
                    angle2ref_ba_or[ang_name_orig] = reference         # default value

            elif ((len(tokens) > 5) and
                  (section_name in ('#bond-bond', '#bond-angle')) and
                  (angle_styles_selected &
                   set(['class2', 'quartic', 'quartic_angle']))):
                if line.lstrip().find('!') == 0:
                    continue
                version = tokens[0]
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = [a for a in map(EncodeAName, tokens[2:5])]
                atom_names = SortByEnds(aorig)
                ang_name_orig = EncodeInteractionName(atom_names, section_is_auto)
                K = ['', '']
                K[0] = tokens[5]
                K[1] = K[0]
                if len(tokens) > 6:
                    K[1] = tokens[6]
                order_reversed = aorig[0] > aorig[-1]
                if order_reversed:
                    K.reverse()
                if (section_name == '#bond-bond'):
                    angle2class2_bb_or[ang_name_orig] = K[0]
                    angle2ver_bb_or[ang_name_orig] = version
                    angle2ref_bb_or[ang_name_orig] = reference
                elif (section_name == '#bond-angle'):
                    angle2class2_ba_or[ang_name_orig] = [k for k in K]
                    angle2ver_ba_or[ang_name_orig] = version
                    angle2ref_ba_or[ang_name_orig] = reference
                if not ang_name_orig in angle2params_or:
                    angle_is_secondary_or[ang_name_orig] = True   #only cross terms have been defined so far
                    angle2params_or[ang_name_orig] = ['0.0', '0.0', '0.0', '0.0']  # default value
                    angle2ver_or[ang_name_orig] = version
                    angle2ref_or[ang_name_orig] = reference
                    angle2priority_or[ang_name_orig] = 0.0









        sys.stderr.write("parsing file pass4: look for dihedrals(torsions) and impropers(out_of_plane)...")

        for iline in range(0, len(lines)):
            line = lines[iline]
            sys.stderr.write('line=\"' + line.strip() + '\"\n')
            tokens = SplitQuotedString(line.strip(),
                                       quotes='',
                                       comment_char='>')
            #sys.stderr.write('tokens = ' + str(tokens) + '\n')
            if line.lstrip().find('!') == 0 and tokens[0] != '!Ver':
                continue


            if line.lstrip(' ').find('#') == 0:
                #sys.stderr.write('allowed_section_names = ' +
                #                 str(allowed_section_names) + '\n')
                if (tokens[0] in allowed_section_names):
                    section_name = tokens[0]
                    section_is_auto = tokens[-1].endswith('_auto')
                    tokens_after_section_name = tokens[1:]
                    sys.stderr.write(' encountered section \"'+tokens[0]+'\"\n')
                    continue
                elif (not tokens[0] in ('#version','#define')):
                    raise InputError('Error: Line# '+str(iline) +'\n'
                                     '       Unrecognized section name:\n'
                                     '       \"' + tokens[0] + '\"\n')




            elif (len(tokens) > 8) and (section_name == '#torsion_1'):
                if line.lstrip().find('!') == 0:
                    continue
                atom_names = SortByEnds(map(EncodeAName, tokens[2:6]))
                dihedral_name = EncodeInteractionName(atom_names, section_is_auto)
                dihedral2ver[dihedral_name] = tokens[0]
                dihedral2ref[dihedral_name] = tokens[1]
                dihedral2priority_or[dihedral_name] = \
                    DetermineNumericPriority(section_is_auto,
                                             tokens[2:6],
                                             float(dihedral2ver[dihedral_name]))
                dihedral_is_secondary_or[dihedral_name] = False
                dihedral2priority[dihedral_name] = \
                    (section_is_auto,
                     dihedral_is_secondary_or[dihedral_name],
                     dihedral2priority_or[dihedral_name])
                K = tokens[6]
                n = tokens[7]
                d = tokens[8]

                w = '0.0'  #ignore: this is only used by the CHARMM force field

                if (dihedral_styles_selected & set(['charmm','torsion_1'])):
                    dihedral_styles.add('charmm')
                    dihedral2style[dihedral_name] = 'charmm'
                    #dihedral2params_or[dihedral_name] = [K,n,d,w]
                    dihedral2params[dihedral_name] = (K+' '+n+' '+d+' '+w)
                elif (dihedral_styles_selected & set(['class2','torsion_3'])):
                    # Then this is a special case of the class2 angle
                    # lacking the higher terms in the Fourier series
                    dihedral_styles.add('class2')
                    dihedral2style[dihedral_name] = 'class2'
                    dihedral2params_or[dihedral_name] = [K,d,0,0,0,0]




            elif ((len(tokens) > 7) and (section_name == '#torsion_3')
                  and (dihedral_styles_selected & set(['class2','torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                dihedral_styles.add('class2')
                atom_names = SortByEnds(map(EncodeAName, tokens[2:6]))
                dih_name_orig = EncodeInteractionName(atom_names, section_is_auto)
                version = tokens[0]
                reference = tokens[1]
                dihedral2priority_or[dih_name_orig] = \
                    DetermineNumericPriority(section_is_auto,
                                             tokens[2:6],
                                             float(version))
                dihedral_is_secondary_or[dih_name_orig] = False
                #dihedral2priority[dih_name_orig] = \
                #    (section_is_auto,
                #     dihedral_is_secondary_or[dih_name_orig],
                #     dihedral2priority_or[dih_name_orig])
                V1 = tokens[6]
                phi0_1 = tokens[7]
                V2 = phi0_2 = V3 = phi0_3 = '0.0'
                if len(tokens) > 9:
                    V2 = tokens[8]
                    phi0_2 = tokens[9]
                if len(tokens) > 11:
                    V3 = tokens[10]
                    phi0_3 = tokens[11]
                dihedral2style_or[dih_name_orig] = 'class2'
                dihedral2ver_or[dih_name_orig] = version
                dihedral2ref_or[dih_name_orig] = reference
                dihedral2params_or[dih_name_orig] = [V1, phi0_1, V2, phi0_2, V3, phi0_3]
                # default values for cross terms:
                if not dih_name_orig in dihedral2class2_mbt_or:
                    dihedral2class2_mbt_or[dih_name_orig] = ['0.0','0.0','0.0']  # default value
                    dihedral2ver_mbt_or[dih_name_orig] = version
                    dihedral2ref_mbt_or[dih_name_orig] = reference
                if not dih_name_orig in dihedral2class2_ebt_or:
                    dihedral2class2_ebt_or[dih_name_orig] = [['0.0','0.0','0.0'],['0.0','0.0','0.0']] # default value
                    dihedral2ver_ebt_or[dih_name_orig] = version
                    dihedral2ref_ebt_or[dih_name_orig] = reference
                if not dih_name_orig in dihedral2class2_bb13_or:
                    dihedral2class2_bb13_or[dih_name_orig] = '0.0' # default value
                    dihedral2ver_bb13_or[dih_name_orig] = version
                    dihedral2ref_bb13_or[dih_name_orig] = reference
                if not dih_name_orig in dihedral2class2_at_or:
                    dihedral2class2_at_or[dih_name_orig] = [['0.0','0.0','0.0'],['0.0','0.0','0.0']] # default value
                    dihedral2ver_at_or[dih_name_orig] = version
                    dihedral2ref_at_or[dih_name_orig] = reference
                if not dih_name_orig in dihedral2class2_aat_or:
                    dihedral2class2_aat_or[dih_name_orig] = '0.0' # default value
                    dihedral2ver_aat_or[dih_name_orig] = version
                    dihedral2ref_aat_or[dih_name_orig] = reference





            elif ((len(tokens) > 6) and (section_name == '#middle_bond-torsion_3')
                  and (dihedral_styles_selected & set(['class2','torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                dihedral_styles.add('class2')
                version = tokens[0]
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = [a for a in map(EncodeAName, tokens[2:6])]
                atom_names = SortByEnds(aorig)

                Fmbt = [tokens[6], '0.0', '0.0']
                if len(tokens) > 7:
                    Fmbt[1] = tokens[7]
                if len(tokens) > 8:
                    Fmbt[2] = tokens[8]

                dih_name_orig = EncodeInteractionName(atom_names, section_is_auto)
                                                         
                #sys.stderr.write('DEBUG: (a2,a3) = '+str((a2,a3))+', '
                #                 ' (b1,b2) = '+str(batoms)+'\n')
                dihedral2style[dih_name_orig] = 'class2'
                dihedral2class2_mbt_or[dih_name_orig] = [F for F in Fmbt]
                dihedral2ver_mbt_or[dih_name_orig] = version
                dihedral2ref_mbt_or[dih_name_orig] = reference
                if not dih_name_orig in dihedral2params_or:
                    dihedral_is_secondary_or[dih_name_orig] = True   #only cross terms have been defined so far
                    dihedral2params_or[dih_name_orig] = ['0.0', '0.0', '0.0', '0.0', '0.0', '0.0']
                    dihedral2ver_or[dih_name_orig] = version
                    dihedral2ref_or[dih_name_orig] = reference
                    dihedral2priority_or[dih_name_orig] = 0.0




            elif ((len(tokens) > 6) and
                  (section_name in ('#end_bond-torsion_3',
                                    '#bond-bond_1_3')) and
                  (dihedral_styles_selected &
                   set(['class2', 'torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                dihedral_styles.add('class2')
                version = tokens[0]
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = [a for a in map(EncodeAName, tokens[2:6])]
                atom_names = SortByEnds(aorig)

                dih_name_orig = EncodeInteractionName(atom_names, section_is_auto)

                dihedral2style[dih_name_orig] = 'class2'
                if section_name == '#end_bond-torsion_3':
                    Febt = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
                    Febt[0][0] = tokens[6]
                    if len(tokens) > 7:
                        Febt[0][1] = tokens[7]
                    if len(tokens) > 8:
                        Febt[0][2] = tokens[8]
                    Febt[1][0] = Febt[0][0]
                    Febt[1][1] = Febt[0][1]
                    Febt[1][2] = Febt[0][2]
                    if len(tokens) > 9:
                        Febt[1][0] = tokens[9]
                    if len(tokens) > 10:
                        Febt[1][1] = tokens[10]
                    if len(tokens) > 11:
                        Febt[1][2] = tokens[11]
                    order_reversed = aorig[0] > aorig[-1]
                    if order_reversed:
                            Febt.reverse()
                    dihedral2class2_ebt_or[dih_name_orig] = [ [F_ij for F_ij in F_i] for F_i in Febt] #deep copy of Febt[][]
                    dihedral2ver_ebt_or[dih_name_orig] = version
                    dihedral2ref_ebt_or[dih_name_orig] = reference

                elif section_name == '#bond-bond_1_3':
                    Kbb13 = tokens[6]
                    #dihedral2ver_bb13[dih_name_orig] = version
                    dihedral2class2_bb13_or[dih_name_orig] = Kbb13
                    dihedral2ver_bb13_or[dih_name_orig] = version
                    dihedral2ref_bb13_or[dih_name_orig] = reference
                else:
                    assert(False)
                if not dih_name_orig in dihedral2params_or:
                    dihedral_is_secondary_or[dih_name_orig] = True   #only cross terms have been defined so far
                    dihedral2params_or[dih_name_orig] = ['0.0', '0.0', '0.0', '0.0', '0.0', '0.0']
                    dihedral2ver_or[dih_name_orig] = version
                    dihedral2ref_or[dih_name_orig] = reference
                    dihedral2priority_or[dih_name_orig] = 0.0










            elif ((len(tokens) > 6) and
                  (section_name in ('#angle-torsion_3',
                                    '#angle-angle-torsion_1')) and
                  (dihedral_styles_selected &
                   set(['class2', 'torsion_3']))):
                if line.lstrip().find('!') == 0:
                    continue
                dihedral_styles.add('class2')
                version = tokens[0]
                reference = tokens[1]
                if line.lstrip().find('!') == 0:
                    continue
                aorig = [a for a in map(EncodeAName, tokens[2:6])]
                atom_names = SortByEnds(aorig)

                dih_name_orig = EncodeInteractionName(atom_names, section_is_auto)

                dihedral2style[dih_name_orig] = 'class2'

                if section_name == '#angle-torsion_3':
                    Fat = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
                    Fat[0][0] = tokens[6]
                    if len(tokens) > 7:
                        Fat[0][1] = tokens[7]
                    if len(tokens) > 8:
                        Fat[0][2] = tokens[8]
                    Fat[1][0] = Fat[0][0]
                    Fat[1][1] = Fat[0][1]
                    Fat[1][2] = Fat[0][2]
                    if len(tokens) > 9:
                        Fat[1][0] = tokens[9]
                    if len(tokens) > 10:
                        Fat[1][1] = tokens[10]
                    if len(tokens) > 11:
                        Fat[1][2] = tokens[11]
                    order_reversed = aorig[0] > aorig[-1]
                    if order_reversed:
                        Fat.reverse()
                        Fat[0].reverse()
                        Fat[1].reverse()
                    dihedral2class2_at_or[dih_name_orig] = [ [F_ij for F_ij in F_i] for F_i in Fat] #deep copy of Fat
                    dihedral2ver_at_or[dih_name_orig] = version
                    dihedral2ref_at_or[dih_name_orig] = reference
                elif section_name == '#angle-angle-torsion_1':
                    Kaat = tokens[6]
                    dihedral2class2_aat_or[dih_name_orig] = Kaat
                    dihedral2ver_aat_or[dih_name_orig] = version
                    dihedral2ref_aat_or[dih_name_orig] = reference
                else:
                    assert(False)

                if not dih_name_orig in dihedral2params_or:
                    dihedral_is_secondary_or[dih_name_orig] = True   #only cross terms have been defined so far
                    dihedral2params_or[dih_name_orig] = ['0.0', '0.0', '0.0', '0.0', '0.0', '0.0'] # default value
                    dihedral2ver_or[dih_name_orig] = version
                    dihedral2ref_or[dih_name_orig] = reference
                    dihedral2priority_or[dih_name_orig] = 0.0










            elif ((len(tokens) > 8) and (section_name == '#out_of_plane')
                   and (improper_styles_selected & set(['cvff','out_of_plane']))):
                if line.lstrip().find('!') == 0:
                    continue
                improper_styles.add('cvff')
                aorig = [a for a in map(EncodeAName, tokens[2:6])]
                atom_names,_ignore  = OOPImproperNameSort(tokens[2:6])
                improper_name = EncodeInteractionName(atom_names, section_is_auto)
                imsym = improper_symmetry_subgraph[improper_name] = 'cenJflipIL'
                subgraph2impname['cenJflipIL'].add(improper_name) CONTINUEHERE
                improper2ver[imsym][improper_name] = tokens[0]
                improper2ref[imsym][improper_name] = tokens[1]
                improper2priority_or[imsym][improper_name] = \
                     DetermineNumericPriority(section_is_auto,
                                              tokens[2:6],
                                              float(improper2ver[imsym][improper_name]))
                improper_is_secondary_or[imsym][imp_name_orig] = False
                improper2priority[imsym][improper_name] = \
                    (section_is_auto,
                     improper_is_secondary_or[imsym][imp_name_orig],
                     improper2priority_or[imsym][improper_name])
                K = tokens[6]
                n = tokens[7]
                chi0 = tokens[8]
                improper2style[imsym][improper_name] = 'cvff'
                improper2params[imsym][improper_name] = (Kchi+' '+n+' '+chi0)
                #if improper_style_name == 'cvff':
                #    improper2params[improper_name] = (Kchi+' '+n+' '+chi0)
                #    improper_symmetry_subgraph[improper_name] = 'cenJswapIL'


            elif ((len(tokens) > 7) and (section_name == '#wilson_out_of_plane')
                  and (improper_styles_selected and set(['class2','wilson_out_of_plane']))):
                if line.lstrip().find('!') == 0:
                    continue
                improper_styles.add('class2')
                sys.stderr.write('tokens = ' + str(tokens) + '\n')

                version = tokens[0]
                reference = tokens[1]
                aorig = [a for a in map(EncodeAName, tokens[2:6])]

                # To avoid redundancy, it is necessary to order the atoms
                # in the interaction so that two equivalent ways of ordering
                # the atoms in an improper interaction do not get misinterpreted
                # as two different types of improper interactions.  So we sort
                # the 3 "leaf" atoms surrounding the central "hub" by name.

                atom_names, permutation = Class2ImproperNameSort(tokens[2:6])

                # This will effect the formula for the energy.
                # (specifically the "chi0" parameter)
                # When we lookup the various cross-term interactions for that
                # same improper interaction, we will be sure to sort them
                # in the same way to make sure those interactions are
                # associated with the same improper interaction.

                imp_name_orig = EncodeInteractionName(atom_names, section_is_auto)
                #improper_symmetry_subgraph_or[improper_name] = 'dihedrals_nosym'  (<--no)
                imsym = improper_symmetry_subgraph_or[imp_name_orig] = 'cenJsortIKL'
                improper2ver_or[imsym][imp_name_orig] = version
                improper2ref_or[imsym][imp_name_orig] = reference
                improper2priority_or[imsym][imp_name_orig] = \
                     DetermineNumericPriority(section_is_auto,
                                              tokens[2:6],
                                              float(improper2ver_or[imp_name_orig]))
                improper_is_secondary_or[imsym][imp_name_orig] = False
                #improper2priority[imp_name_orig] = \
                #    (section_is_auto,
                #     improper_is_secondary_or[imp_name_orig],
                #     improper2priority_or[imp_name_orig])
                K = tokens[6]
                chi0 = tokens[7]

                if Parity(permutation) != 0:
                    # Each time the order of a pair of atoms is swapped in
                    # the interaction, all 3 of the "X" (chi) angles change sign
                    # The formula for the ordinary term in the improper
                    # interaction is Ei = K*((Xijkl + Xkjli + Xljik)/3 - chi0)^2
                    # This formula is invariant if we change the sign of all
                    # Xijkl, Xkjli, Xljik, chi0
                    # Hence, we can account for a change in atom order by
                    # changing the sign of the "chi0" parameter.
                    # We calculate the "Parity" of the permutation (ie whether
                    # the permutation has an even or odd number of swaps)
                    # and multiply chi0 by -1 for each swap.
                    # It's not clear if this is necessary since in practice
                    # the "chi0" parameter is usually zero.

                    chi0 = str(-1.0*float(chi0))  # same as ('-' + chi0)

                improper2style_or[imsym][imp_name_orig] = 'class2'
                improper2params_or[imsym][imp_name_orig] = [K, chi0]
                #improper2params[imp_name_orig] = K + ' ' + chi0
                # default values for cross terms:
                if not imp_name_orig in improper2class2_aa_or:
                    improper2class2_aa_or[imsym][imp_name_orig] = '0.0' #(default)
                    improper2ver_aa_or[imsym][imp_name_orig] = version
                    improper2ref_aa_or[imsym][imp_name_orig] = reference
                    # Initially, set all of the angle-angle cross terms to zero
                    # Start with the first cross term between aorig[0],aorig[1],aorig[2] & aorig[2],aorig[1],aorig[3]
                    improper2cross[imp_name_orig][ImCrossTermID([aorig[0],aorig[1],aorig[2],aorig[3]])] = '0.0'
                    # ...then cyclically permute the 3 "leaf" atoms (aorig[0], aorig[2], aorig[3]) around the "hub" atom (aorig[1])
                    improper2cross[imp_name_orig][ImCrossTermID([aorig[2],aorig[1],aorig[3],aorig[0]])] = '0.0'
                    improper2cross[imp_name_orig][ImCrossTermID([aorig[3],aorig[1],aorig[0],aorig[2]])] = '0.0'

            elif ((len(tokens) > 6) and (section_name == '#angle-angle')
                  and (improper_styles_selected and set(['class2','wilson_out_of_plane']))):
                if line.lstrip().find('!') == 0:
                    continue
                improper_styles.add('class2')
                version = tokens[0]
                reference = tokens[1]
                aorig = [a for a in map(EncodeAName, tokens[2:6])]
                atom_names, permutation = Class2ImproperNameSort(tokens[2:6])
                imp_name_orig = EncodeInteractionName(atom_names, section_is_auto)
                imsym = improper_symmetry_subgraph_or[imp_name_orig] = 'cenJsortIKL'
                improper2ver_aa_or[imsym][imp_name_orig] = version
                improper2ref_aa_or[imsym][imp_name_orig] = reference
                K = tokens[6]
                improper2style_or[imsym][imp_name_orig] = 'class2'
                if not imp_name_orig in improper2params_or:
                    improper_is_secondary_or[imsym][imp_name_orig] = True   #only cross terms have been defined so far
                    improper2params_or[imsym][imp_name_orig] = ['0.0', '0.0']
                    improper2ver_or[imsym][imp_name_orig] = version
                    improper2ref_or[imsym][imp_name_orig] = reference
                    improper2priority_or[imsym][imp_name_orig] = 0.0
                if not imp_name_orig in improper2cross:
                    # then initialize all of the cross terms to zero
                    improper2cross[imp_name_orig][ImCrossTermID([aorig[0],aorig[1],aorig[2],aorig[3]])] = '0.0'
                    # ...then cyclically permute the 3 "leaf" atoms (aorig[0], aorig[2], aorig[3]) around the "hub" atom (aorig[1])
                    improper2cross[imp_name_orig][ImCrossTermID([aorig[2],aorig[1],aorig[3],aorig[0]])] = '0.0'
                    improper2cross[imp_name_orig][ImCrossTermID([aorig[3],aorig[1],aorig[0],aorig[2]])] = '0.0'
                #improper2class2_aa_or[imp_name_orig] = K   (not needed)
                improper2cross[imp_name_orig][ImCrossTermID(aorig)] = K

            elif (len(tokens) > 0) and (section_name == '#out_of_plane-out_of_plane'):
                if line.lstrip().find('!') == 0:
                    continue
                display_OOP_OOP_warning = True

            elif (len(tokens) > 0) and (section_name == '#torsion-torsion_1'):
                if line.lstrip().find('!') == 0:
                    continue
                display_torsion_torsion_1_warning = True

            elif section_name == '#templates':
                #if line.lstrip().find('!') == 0:
                #    continue
                lines_templates.append(line)

            elif section_name == '#reference':
                if line.lstrip().find('!') == 0:
                    continue
                if len(tokens_after_section_name) > 0:
                    ref_number = int(tokens_after_section_name[0])
                if len(line.strip()) > 0:
                    lines_references[ref_number].append(line)



            """
             --- these next few lines of code appear to be unnecessary.
             --- I'll probably delete this code in a later version
            elif (len(tokens) > 3) and (section_name == '#hbond_definition'):
                hbondID = tokens[1]
                if tokens[2] == 'distance':
                    hbond2distance[hbondID] = tokens[3]
                if tokens[2] == 'angle':
                    hbond2angle[hbondID] = tokens[3]
                if tokens[2] == 'donors':
                    hbond2donors[hbondID] = map(EncodeAName, tokens[2:])
                if tokens[2] == 'acceptors':
                    hbond2acceptors[hbondID] = map(EncodeAname(),tokens[2:])
            """


        if display_OOP_OOP_warning:
            lines_warnings.append('###########################################################\n'
                             '# WARNING\n'
                             '#      ALL \"out-of-plane_out-of_plane\" INTERACTIONS ARE IGNORED.\n'
                             '#      CHECK THAT THESE TERMS ARE NEGLEGIBLY SMALL.\n'
                             '#      \"out-of-plane_out-of_plane\" interactions are not yet supported in LAMMPS\n'
                                  '#      (...as of 2017-10-13)  There is no way that moltemplate can produce\n'
                             '#      LAMMPS compatible parameter files for these interactions.\n'
                             '###########################################################\n')

        if display_torsion_torsion_1_warning:
            lines_warnings.append('###########################################################\n'
                             '# WARNING\n'
                             '#      ALL \"torsion_torsion_1\" INTERACTIONS ARE IGNORED.\n'
                             '#      CHECK THAT THESE TERMS ARE NEGLEGIBLY SMALL.\n'
                             '#      \"torsion_torsion_1\" interactions are not yet supported in LAMMPS\n'
                                  '#      (...as of 2017-10-13)  There is no way that moltemplate can produce\n'
                             '#      LAMMPS compatible parameter files for these interactions.\n'
                             '###########################################################\n')


        sys.stderr.write(' done.\n'
                         'building lookup tables...')







        """
         --- these next few lines of code appear to be unnecessary.
         --- I'll probably delete them eventually
        if len(hbond2params) > 0:
            sys.stdout.write('\n\n  write_once("In Settings") {\n')
            if hbond_style == 'hbond/dreiding/lj':
                for hbondID, angle in hbond2angle:
                    hbond2params[hbondID] =  hbond2distance[hbondID]+' '+hbond2angle[hbondID]  ##<--this is not correct
            for hbondID, params in hbond2params:
                for donor in hbond2donors[hbondID]:
                    for acceptor in hbond2acceptors[hbondID]:
                        for hydrogen in hbond2hydrogens[hbondID]:
                            sys.stdout.write('pair_coeff @atom:'+donor+' @atom:'+acceptor+' '+hbond_style+' @atom:'+hydrogen+' i '+params+'\n')
            sys.stdout.write('  }   # (DREIDING style H-bond parameters)\n\n\n')
        """






        sys.stderr.write(" done.\n")
        sys.stderr.write("Trying all combinations of atom types...")







        ##################### POST-PROCESSING ########################





        for ang_name_orig in angle2params_or:

            is_auto = (ang_name_orig.find('auto_') == 0)

            atom_names = ExtractANames(ang_name_orig)

            num_angles = 0

            atom_combos = [set([]), set([]), set([])]

            # We must consider every possible combination of atom types
            # which satisfy BOTH angle_equivalences and bond_equivalences.
            # ...AND we must consider BOTH regular AND auto equivalences.
            # For each combination generate a separate @angle interaction.
            # (I fear this will make the resulting .LT file large.)

            # Use different auto equivalence lookup tables for different
            # atoms in the interaction. (ie the "center" and "end" atoms)
            auto_angle2atom = [auto_angleend2atom,
                               auto_anglecenter2atom,
                               auto_angleend2atom]

            for i in range(0, 3):
                angle_atom_name = atom_names[i]
                sys.stderr.write('DEBUG: angle_atom_name = '+angle_atom_name+'\n')
                if not is_auto:
                    assert(angle_atom_name[-1] != '_')
                    # assume regular equivalences when looking up atom types
                    sys.stderr.write('DEBUG: equiv_angle2atom['+angle_atom_name+'] = '+
                                     str(equiv_angle2atom[angle_atom_name])+'\n')
                    for a in equiv_angle2atom[angle_atom_name]:
                        atom_combos[i].add(a)
                else:
                    #assert((angle_atom_name[-1] == '_') or (angle_atom_name[0] == '*'))  (<--some exceptions. don't assert this)

                    # assume "auto" equivalences when looking up atom types
                    sys.stderr.write('DEBUG: auto_angle2atom['+str(i)+']['+angle_atom_name+'] = \n'
                                     '       '+str(equiv_angle2atom[i][angle_atom_name])+'\n')
                    for a in auto_angle2atom[i][angle_atom_name]:
                        atom_combos[i].add(a)

            found_at_least_one = False
            #for a1 in atom_combos[0]:
            for a1 in sorted(list(atom_combos[0])):
                #for a2 in atom_combos[1]:
                for a2 in sorted(list(atom_combos[1])):
                    #sys.stderr.write('atom2auto_bond = '+str(atom2auto_bond)+'\n')
                    bond_data1 = LookupBondLength(a1, a2,
                                                  atom2equiv_bond,
                                                  bond2r0,
                                                  atom2auto_bond,
                                                  bond2r0_auto)
                    if bond_data1 == None: # Save time by continuing only if a
                        continue           # bond was defined between a1 and a2

                    #for a3 in atom_combos[2]:
                    for a3 in sorted(list(atom_combos[2])):
                        bond_data2 = LookupBondLength(a2, a3,
                                                      atom2equiv_bond,
                                                      bond2r0,
                                                      atom2auto_bond,
                                                      bond2r0_auto)
                        if bond_data2 == None:
                            continue

                        #bond lengths:
                        r0s = [0.0, 0.0]
                        #equivalent atom names used to lookup the bonds:
                        batoms = [['', ''], ['', '']]
                        #were "auto" equivalences needed to lookup the bond length?
                        b_is_auto = [False, False]
                        r0s[0], batoms[0], b_is_auto[0] = bond_data1
                        r0s[1], batoms[1], b_is_auto[1] = bond_data2
                        order_reversed = aorig[0] > aorig[-1]
                        if order_reversed:
                            batoms.reverse()
                            batoms[0].reverse()
                            batoms[1].reverse()
                            b_is_auto.reverse()
                        ang_name_full = (ang_name_orig + ',' + 
                                         EncodeInteractionName(batoms[0], b_is_auto[0]) + ',' +
                                         EncodeInteractionName(batoms[1], b_is_auto[1]))


                        #sys.stderr.write('DEBUG: (a1,a2,a3) = '+str((a1,a2,a3))+', '
                        #                 ' (b11,b12,b21,b22) = '+str(batoms)+'\n')
                        angle2ref_or[ang_name_full] = reference
                        angle2style_or[ang_name_full] = 'class2'
                        theta0_K_params = angle2params_or[ang_name_orig]
                        angle2params[ang_name_full] = ' '.join(theta0_K_params)
                        if ang_name_orig in angle2class2_bb_or:
                            Kbb = angle2class2_bb_or[ang_name_orig]
                            assert(ang_name_orig in angle2ver_bb_or)
                            assert(ang_name_orig in angle2ref_bb_or)
                        else:    #(use default values)
                            Kbb = '0.0'
                            angle2class2_bb_or[ang_name_orig] = Kbb
                            angle2ver_bb_or[ang_name_orig] = angle2ver_or[ang_name_orig]
                            angle2ref_bb_or[ang_name_orig] = angle2ref_or[ang_name_orig]
                        angle2class2_bb[ang_name_full] = (Kbb+' '+r0s[0]+' '+r0s[1])
                        angle2priority_bb = \
                            DetermineNumericPriority(is_auto,
                                                     batoms[0] + batoms[1],
                                                     float(angle2ver_bb_or[ang_name_orig]))
                        angle2ver_bb[ang_name_full] = angle2ver_bb_or[ang_name_orig]
                        angle2ref_bb[ang_name_full] = angle2ref_bb_or[ang_name_orig]

                        if ang_name_orig in angle2class2_ba_or:
                            Kba = angle2class2_ba_or[ang_name_orig]
                            assert(ang_name_orig in angle2ver_ba_or)
                            assert(ang_name_orig in angle2ref_ba_or)
                        else:    #(use default values)
                            Kba = ['0.0', '0.0']
                            angle2class2_ba_or[ang_name_orig] = Kba
                            angle2ver_ba_or[ang_name_orig] = angle2ver_or[ang_name_orig]
                            angle2ref_ba_or[ang_name_orig] = angle2ref_or[ang_name_orig]
                        angle2class2_ba[ang_name_full] = (Kba[0]+' '+Kba[1]+' '+r0s[0]+' '+r0s[1])
                        angle2sym_ba = (Kba[0] == Kba[1])
                        angle2priority_ba = \
                            DetermineNumericPriority(is_auto,
                                                     batoms[0] + batoms[1],
                                                     angle2ver_ba_or[ang_name_orig])
                        angle2ver_ba[ang_name_full] = angle2ver_ba_or[ang_name_orig]
                        angle2ref_ba[ang_name_full] = angle2ref_ba_or[ang_name_orig]

                        version = max((angle2ver_or[ang_name_orig],
                                       angle2ver_bb_or[ang_name_orig],
                                       angle2ver_ba_or[ang_name_orig]))
                        angle2ver[ang_name_full] = version
                        angle2ref[ang_name_full] = angle2ref_or[ang_name_orig]
                        angle2style[ang_name_full] = 'class2'
                        angle2priority[ang_name_full] = \
                            (is_auto,
                             angle_is_secondary_or[ang_name_orig],
                             angle2priority_or[ang_name_orig],
                             angle2priority_bb,
                             angle2priority_ba)

                        if num_angles < len(angle2params):
                            sys.stderr.write('DEBUG: '+section_name[1:]+' r0 ('+ang_name_full+') = ('+r0s[0]+', '+r0s[1]+')\n')
                            sys.stderr.write('DEBUG: len(angle2class2_bb) = '+str(len(angle2class2_bb))+'\n')
                            sys.stderr.write('DEBUG: '+section_name[1:]+' r0 ('+ang_name_full+') = ('+r0s[0]+', '+r0s[1]+')\n')
                            #sys.stderr.write('DEBUG: len(angle2class2_ba) = '+str(len(angle2class2_ba))+'\n')
                            num_angles = len(angle2params)

                        if ((not angle2sym_ba)
                            and
                            (atom_names[0] == atom_names[2])):
                            raise InputError('Error: Unsupported angle interaction: \"@angle:'+str(ang_name_orig)+'\"\n'
                                             '       This interaction has symmetric atom names:\n'
                                             ', '.join(atom_names)+'\n'
                                             '       and yet it lacks symmetry in the corresponding force field parameters.\n'
                                             '       (If this is not a mistake in the .frc file, then explain\n'
                                             '        why to andrew so he can fix this.)\n')


                        found_at_least_one = True


            if not found_at_least_one:
                lines_warnings.append('# WARNING: Undefined bond length (r0) in angle: ' +
                                      ' '.join(atom_names)+'\n')
                # Then we were unable to define cross terms for this interaction
                # because at least one of the bond lengths could not be determined.
                # This usually occurs because most of the .FRC files which are
                # in circulation are incomplete.  We have to handle this gracefully.
                ang_name_full = (ang_name_orig + ',X,X,X,X,X,X')
                version = angle2ver_or[ang_name_orig]
                reference = angle2ref_or[ang_name_orig]
                angle2ref[ang_name_full] = reference
                angle2ver[ang_name_full] = version
                angle2style[ang_name_full] = 'class2'
                angle2params[ang_name_full] = ' '.join(angle2params_or[ang_name_orig])
                # substitute zeros for all the cross term interactions
                angle2priority[ang_name_full] = angle2priority_or[ang_name_orig]
                angle2class2_bb[ang_name_full] = '0.0 1.0 1.0'
                angle2ref_bb[ang_name_full] = reference
                angle2ver_bb[ang_name_full] = version
                angle2class2_ba[ang_name_full] = '0.0 0.0 1.0 1.0'
                angle2ref_ba[ang_name_full] = reference
                angle2ver_ba[ang_name_full] = version
            #sys.stderr.write('bond_names = ' + str(bond_names) + '\n')





        ############ POST-PROCESSING DIHEDRALS ###########



        for dih_name_orig in dihedral2params_or:
            #assert(dih_name_orig in dihedral2class2_mbt_or)
            #assert(dih_name_orig in dihedral2class2_ebt_or)
            #assert(dih_name_orig in dihedral2class2_bb13_or)
            #assert(dih_name_orig in dihedral2class2_at_or)
            #assert(dih_name_orig in dihedral2class2_aat_or)

            is_auto = (dih_name_orig.find('auto_') == 0)

            atom_names = ExtractANames(dih_name_orig) 

            num_dihedrals = 0

            atom_combos = [set([]), set([]), set([]), set([])]

            # We must consider every possible combination of atom types
            # which satisfy all three:
            #   dihedral_equivalences
            #   bond_equivalences
            #   angle_equivalences
            # ...AND we must consider BOTH regular AND auto equivalences.
            # For each combination generate a separate @dihedral interaction.
            # (I fear this will make the resulting .LT file large.)

            # Use different auto equivalence lookup tables for different
            # atoms in the interaction. (ie the "center" and "end" atoms)
            auto_dihedral2atom = [auto_dihedralend2atom,
                                  auto_dihedralcenter2atom,
                                  auto_dihedralcenter2atom,
                                  auto_dihedralend2atom]

            for i in range(0, 4):
                dihedral_atom_name = atom_names[i]
                sys.stderr.write('DEBUG: dihedral_atom_name = '+dihedral_atom_name+'\n')
                if not is_auto:
                    assert(dihedral_atom_name[-1] != '_')
                    # assume regular equivalences when looking up atom types
                    sys.stderr.write('DEBUG: equiv_dihedral2atom['+dihedral_atom_name+'] = '+
                                     str(equiv_dihedral2atom[dihedral_atom_name])+'\n')
                    for a in equiv_dihedral2atom[dihedral_atom_name]:
                        atom_combos[i].add(a)
                else:
                    assert((dihedral_atom_name[-1] == '_') or (ange_atom_name[0] == '*'))
                    # assume "auto" equivalences when looking up atom types
                    sys.stderr.write('DEBUG: auto_dihedral2atom['+str(i)+']['+dihedral_atom_name+'] = \n'
                                     '       '+str(equiv_dihedral2atom[i][dihedral_atom_name])+'\n')
                    for a in auto_dihedral2atom[i][dihedral_atom_name]:
                        atom_combos[i].add(a)

            found_at_least_one = False

            #for a1 in atom_combos[0]:
            for a1 in sorted(list(atom_combos[0])):

                #for a2 in atom_combos[1]:
                for a2 in sorted(list(atom_combos[1])):

                    #sys.stderr.write('atom2auto_bond = '+str(atom2auto_bond)+'\n')
                    bond_data12 = LookupBondLength(a1, a2,
                                                   atom2equiv_bond,
                                                   bond2r0,
                                                   atom2auto_bond,
                                                   bond2r0_auto)
                    if bond_data12 == None:
                        # Save time by only continuing if a bond was
                        # found between a1 and a2
                        continue
                    #for a3 in atom_combos[2]:
                    for a3 in sorted(list(atom_combos[2])):
                        bond_data23 = LookupBondLength(a2, a3,
                                                       atom2equiv_bond,
                                                       bond2r0,
                                                       atom2auto_bond,
                                                       bond2r0_auto)
                        if bond_data23 == None:
                            # Save time by only continuing if a bond was
                            # found between a2 and a3
                            continue

                        angle_data123 = LookupBondAngle(a1, a2, a3,
                                                        atom2equiv_angle,
                                                        angle2theta0_or,
                                                        [atom2auto_angleend,
                                                         atom2auto_anglecenter,
                                                         atom2auto_anglecenter],
                                                        angle2theta0_auto_or)
                        if angle_data123 == None:
                            # Save time by only continuing if an angle was
                            # found between a1, a2, a3
                            continue


                        #for a4 in atom_combos[3]:
                        for a4 in sorted(list(atom_combos[3])):
                            bond_data34 = LookupBondLength(a3, a4,
                                                           atom2equiv_bond,
                                                           bond2r0,
                                                           atom2auto_bond,
                                                           bond2r0_auto)
                            if bond_data34 == None:
                                # Save time by only continuing if a bond was
                                # found between a3 and a4
                                continue

                            #rest bond lengths:
                            r0s = [0.0, 0.0, 0,0]
                            #equivalent atom names used to lookup the bonds:
                            batoms = [['', ''], ['', ''], ['','']]
                            #are these bond interactions "auto" interactions?
                            #were "auto" equivalences needed to lookup the bond length?
                            b_is_auto = [False, False, False]
                            r0s[0], batoms[0], b_is_auto[0] = bond_data12
                            r0s[1], batoms[1], b_is_auto[1] = bond_data23
                            r0s[2], batoms[2], b_is_auto[2] = bond_data34

                            angle_data234 = LookupBondAngle(a2, a3, a4,
                                                            atom2equiv_angle,
                                                            angle2theta0_or,
                                                            [atom2auto_angleend,
                                                             atom2auto_anglecenter,
                                                             atom2auto_anglecenter],
                                                            angle2theta0_auto_or)
                            if angle_data234 == None:
                                # Save time by only continuing if an angle was
                                # found between a2, a3, a4
                                continue

                            #rest angles:
                            theta0s = [0.0, 0.0]
                            #equivalent atom names used to lookup angles:
                            aatoms = [['', '',''], ['', '','']]
                            #were "auto" equivalences needed to lookup the bond-angle?
                            a_is_auto = [False, False]
                            theta0s[0], aatoms[0], a_is_auto[0] = angle_data123
                            theta0s[1], aatoms[1], a_is_auto[1] = angle_data234
                            order_reversed = aorig[0] > aorig[-1]
                            if order_reversed:
                                batoms.reverse()
                                batoms[0].reverse()
                                batoms[1].reverse()
                                batoms[2].reverse()
                                b_is_auto.reverse()
                                theta0s.reverse()
                                aatoms.reverse()
                                aatoms[0].reverse()
                                aatoms[1].reverse()
                                a_is_auto.reverse()

                            #if is_auto:
                            dih_name_full = (dih_name_orig + ',' + 
                                             EncodeInteractionName(batoms[0], b_is_auto[0]) + ',' +
                                             EncodeInteractionName(batoms[1], b_is_auto[1]) + ',' +
                                             EncodeInteractionName(batoms[2], b_is_auto[2]) + ',' +
                                             EncodeInteractionName(aatoms[0], a_is_auto[0]) + ',' +
                                             EncodeInteractionName(aatoms[1], a_is_auto[1]))
                            #else:
                            #    assert(batoms[0][1] == batoms[1][0])
                            #    assert(batoms[1][1] == batoms[2][0])
                            #    assert(aatoms[0][1] == aatoms[1][0])
                            #    assert(aatoms[0][2] == aatoms[1][1])
                            #    dih_name_full = dih_name_orig + ',' + \
                            #        EncodeInteractionName([batoms[0][0], batoms[0][1]
                            #                               batoms[2][0], batoms[2][1],
                            #                               aatoms[0][0], aatoms[0][1],
                            #                               aatoms[0][2], aatoms[1][0]],
                            #                               False)

                            ########### Fourier terms ###########
                            #if dih_name_orig in dihedral2param_or:
                            V_phi0_params = dihedral2params_or[dih_name_orig]
                            dihedral2params[dih_name_full] = ' '.join(V_phi0_params)
                            #else:
                            #    dihedral2params[dih_name_full] = '0.0 0.0 0.0 0.0 0.0 0.0'

                            ########### "mbt", "ebt", and "aat" terms ###########
                            # "mbt" terms:
                            if dih_name_orig in dihedral2class2_mbt_or:
                                Fmbt = dihedral2class2_mbt_or[dih_name_orig]
                            else:
                                Fmbt = ['0.0', '0.0', '0.0']
                                dihedral2class2_mbt_or[dih_name_orig] = Fmbt
                                dihedral2ver_mbt_or[dih_name_orig] = dihedral2ver_or[dih_name_orig]
                                dihedral2ref_mbt_or[dih_name_orig] = dihedral2ref_or[dih_name_orig]
                            dihedral2class2_mbt[dih_name_full] = \
                               (Fmbt[0]+' '+Fmbt[1]+' '+Fmbt[2]+' '+r0s[1])
                            dihedral2priority_mbt = \
                                DetermineNumericPriority(is_auto,
                                                         batoms[1],
                                                         float(dihedral2ver_mbt_or[dih_name_orig]))
                            dihedral2ver_mbt[dih_name_full] = dihedral2ver_mbt_or[dih_name_orig]
                            dihedral2ref_mbt[dih_name_full] = dihedral2ref_mbt_or[dih_name_orig]

                            # "ebt" terms:
                            if dih_name_orig in dihedral2class2_ebt_or:
                                Febt = dihedral2class2_ebt_or[dih_name_orig]
                                dihedral2sym_ebt = ((Febt[0][0] == Febt[1][0]) and
                                                    (Febt[0][1] == Febt[1][1]) and
                                                    (Febt[0][2] == Febt[1][2]))
                                                    #and (r0s[0] == r0s[2]))
                            else:
                                Febt = [['0.0','0.0','0.0'], ['0.0','0.0','0.0']]
                                dihedral2class2_ebt_or[dih_name_orig] = Febt
                                dihedral2ver_ebt_or[dih_name_orig] = dihedral2ver_or[dih_name_orig]
                                dihedral2ref_ebt_or[dih_name_orig] = dihedral2ref_or[dih_name_orig]
                                dihedral2sym_ebt = True
                            dihedral2class2_ebt[dih_name_full]= (Febt[0][0] + ' ' +
                                                                 Febt[0][1] + ' ' +
                                                                 Febt[0][2] + ' ' +
                                                                 Febt[1][0] + ' ' +
                                                                 Febt[1][1] + ' ' +
                                                                 Febt[1][2] + ' ' +
                                                                 r0s[0]+' '+r0s[2])

                            dihedral2priority_ebt = \
                                DetermineNumericPriority(is_auto,
                                                         batoms[0] + batoms[2],
                                                         float(dihedral2ver_ebt_or[dih_name_orig]))
                            dihedral2ver_ebt[dih_name_full] = dihedral2ver_ebt_or[dih_name_orig]
                            dihedral2ref_ebt[dih_name_full] = dihedral2ref_ebt_or[dih_name_orig]

                            #(Note:  large atom_priority number <==> low priority
                            # Only one of the atom priority numbers should be > 0)

                            # "bb13" terms:
                            if dih_name_orig in dihedral2class2_bb13_or:
                                Kbb13 = dihedral2class2_bb13_or[dih_name_orig]
                                #dihedral2sym_bb13 = (r0s[0] == r0s[2])
                                dihedral2sym_bb13 = True
                            else:
                                Kbb13 = '0.0'
                                dihedral2class2_bb13_or[dih_name_orig] = Kbb13
                                dihedral2ver_bb13_or[dih_name_orig] = dihedral2ver_or[dih_name_orig]
                                dihedral2ref_bb13_or[dih_name_orig] = dihedral2ref_or[dih_name_orig]
                                dihedral2sym_bb13 = True

                            dihedral2class2_bb13[dih_name_full] = (Kbb13+' '+r0s[0]+' '+r0s[2])
                            dihedral2priority_bb13 = \
                                DetermineNumericPriority(is_auto,
                                                         batoms[0] + batoms[2],
                                                         float(dihedral2ver_bb13_or[dih_name_orig]))
                            dihedral2ver_bb13[dih_name_full] = dihedral2ver_bb13_or[dih_name_orig]
                            dihedral2ref_bb13[dih_name_full] = dihedral2ref_bb13_or[dih_name_orig]


                            ########### "at" and "aat" terms ###########
                            # "at" terms:
                            if dih_name_orig in dihedral2class2_at_or:
                                Fat = dihedral2class2_at_or[dih_name_orig]
                                dihedral2sym_at = ((Fat[0][0] == Fat[1][0]) and
                                                   (Fat[0][1] == Fat[1][1]) and
                                                   (Fat[0][2] == Fat[1][2]))
                                                   #and (theta0[0] == theta0[1]))
                            else:
                                Fat = [['0.0','0.0','0.0'], ['0.0','0.0','0.0']]
                                dihedral2class2_at_or[dih_name_orig] = Fat
                                dihedral2ver_at_or[dih_name_orig] = dihedral2ver_or[dih_name_orig]
                                dihedral2ref_at_or[dih_name_orig] = dihedral2ref_or[dih_name_orig]
                                dihedral2sym_at = True
                            dihedral2class2_at[dih_name_full] = \
                                (Fat[0][0] + ' ' +
                                 Fat[0][1] + ' ' +
                                 Fat[0][2] + ' ' +
                                 Fat[1][0] + ' ' +
                                 Fat[1][1] + ' ' +
                                 Fat[1][2] + ' ' +
                                 theta0s[0] + ' ' +
                                 theta0s[1])                            
                            dihedral2priority_at = \
                                DetermineNumericPriority(is_auto,
                                                         aatoms[0] + aatoms[1],
                                                         float(dihedral2ver_at_or[dih_name_orig]))
                            dihedral2ver_at[dih_name_full] = dihedral2ver_at_or[dih_name_orig]
                            dihedral2ref_at[dih_name_full] = dihedral2ref_at_or[dih_name_orig]


                            # "aat" terms:
                            if dih_name_orig in dihedral2class2_aat_or:
                                Kaat = dihedral2class2_aat_or[dih_name_orig]
                                #dihedral2sym_aat = (theta0[0] == theta0[1])
                                dihedral2sym_aat = True
                            else:
                                Kaat = '0.0'
                                dihedral2class2_aat_or[dih_name_orig] = Kaat
                                dihedral2ver_aat_or[dih_name_orig] = dihedral2ver_or[dih_name_orig]
                                dihedral2ref_aat_or[dih_name_orig] = dihedral2ref_or[dih_name_orig]
                                dihedral2sym_aat = True
                            dihedral2class2_aat[dih_name_full] = \
                                        (Kaat+' '+theta0s[0]+' '+theta0s[1])
                            dihedral2priority_aat = \
                                DetermineNumericPriority(is_auto,
                                                         aatoms[0] + aatoms[1],
                                                         float(dihedral2ver_aat_or[dih_name_orig]))
                            dihedral2ver_aat[dih_name_full] = dihedral2ver_aat_or[dih_name_orig]
                            dihedral2ref_aat[dih_name_full] = dihedral2ref_aat_or[dih_name_orig]

                            if len(dihedral2params) > num_dihedrals:
                                sys.stderr.write('DEBUG: dihedral['+dih_name_full+']:\n'
                                                 '(r12,r23,r34) = ('
                                                 +r0s[0]+','+r0s[1]+','+r0s[2]+') \n'
                                                 '(theta123,theta234) = ('
                                                 +theta0s[0]+','+theta0s[1]+') \n')
                                sys.stderr.write('DEBUG: num_dihedrals = len(dihedral2params) = '
                                                 +str(len(dihedral2params))+'\n')
                            version = max((dihedral2ver_or[dih_name_orig],
                                           dihedral2ver_mbt_or[dih_name_orig],
                                           dihedral2ver_ebt_or[dih_name_orig],
                                           dihedral2ver_bb13_or[dih_name_orig],
                                           dihedral2ver_at_or[dih_name_orig],
                                           dihedral2ver_aat_or[dih_name_orig]))

                            dihedral2style[dih_name_full] = 'class2'
                            dihedral2ver[dih_name_full] = version
                            dihedral2ref[dih_name_full] = dihedral2ref_or[dih_name_orig]
                            dihedral2priority[dih_name_full] = \
                                (is_auto,
                                 dihedral_is_secondary_or[dih_name_orig],
                                 dihedral2priority_or[dih_name_orig],
                                 dihedral2priority_mbt,
                                 dihedral2priority_ebt,
                                 dihedral2priority_bb13,
                                 dihedral2priority_at,
                                 dihedral2priority_aat)

                            num_dihedrals = len(dihedral2params)

                            if ((not (dihedral2sym_ebt and
                                      #dihedral2sym_mbt and
                                      # (note: symmetry doesn't make sense for mbt)
                                      dihedral2sym_at and
                                      dihedral2sym_aat and
                                      dihedral2sym_bb13))
                                and
                                ((atom_names[0] == atom_names[3]) and
                                 (atom_names[1] == atom_names[2]))):
                                raise InputError('Error: Unsupported dihedral interaction: \"@dihedral:'+str(dih_name_orig)+'\"\n'
                                                 '       This interaction has symmetric atom names:\n'+
                                                 ', '.join(atom_names)+'\n'+
                                                 '       and yet it lacks symmetry in the corresponding force field parameters.\n'+
                                                 '       (If this is not a mistake in the .frc file, then explain\n'+
                                                 '        why to andrew so he can fix this.)\n')

                            found_at_least_one = True


            #sys.stderr.write('DEBUG: number of interactions = '+str(len(dihedral2class2_bb))+'\n')
            if not found_at_least_one:
                lines_warnings.append('# WARNING: Undefined bond length (r0) or rest angle (theta0) in dihedral: ' +
                                      #'#          the dihedral interaction between: ' +
                                      ' '.join(atom_names)+'\n')
                # Then we were unable to define cross terms for this interaction because
                # at least one of the bond lengths or bond angles could not be determined.
                # This usually occurs because most of the .FRC files which are
                # in circulation are incomplete.  We have to handle this gracefully.
                dih_name_full = (dih_name_orig + ',X,X,X,X,X,X,X,X,X,X,X,X')
                reference = dihedral2ref_or[dih_name_orig]
                version = dihedral2ver_or[dih_name_orig]
                dihedral2ref[dih_name_full] = reference
                dihedral2ver[dih_name_full] = version
                dihedral2style[dih_name_full] = 'class2'
                dihedral2priority[dih_name_full] = dihedral2priority_or[dih_name_orig]
                dihedral2params[dih_name_full] = ' '.join(dihedral2params_or[dih_name_orig])
                # substitute zeros for all the cross term interactions

                dihedral2class2_mbt[dih_name_full] = '0.0 0.0 0.0 1.0'
                dihedral2ref_mbt[dih_name_full] = reference
                dihedral2ver_mbt[dih_name_full] = version

                dihedral2class2_ebt[dih_name_full] = '0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0'
                dihedral2ref_ebt[dih_name_full] = reference
                dihedral2ver_ebt[dih_name_full] = version

                dihedral2class2_bb13[dih_name_full] = '0.0 1.0 1.0'
                dihedral2ref_bb13[dih_name_full] = reference
                dihedral2ver_bb13[dih_name_full] = version

                dihedral2class2_at[dih_name_full] = '0.0 0.0 0.0 0.0 0.0 0.0 120.0 120.0'
                dihedral2ref_at[dih_name_full] = reference
                dihedral2ver_at[dih_name_full] = version

                dihedral2class2_aat[dih_name_full] = '0.0 120.0 120.0'
                dihedral2ref_aat[dih_name_full] = reference
                dihedral2ver_aat[dih_name_full] = version








        ############ POST-PROCESSING IMPROPERS ###########




        imsym = 'cenJsortIKL'
        for imp_name_orig in improper2cross[imsym]:

            if improper2style_or[imsym][imp_name_orig] != 'class2':
                continue

            assert(imp_name_orig in improper2params_or[imsym])
            assert(imp_name_orig in improper2class2_aa_or[imsym])

            is_auto = (imp_name_orig.find('auto') == 0)

            atom_names = ExtractANames(imp_name_orig) 

            num_impropers = 0

            atom_combos = [set([]), set([]), set([]), set([])]

            # We must consider every possible combination of atom types
            # which satisfy both:
            #   improper_equivalences
            #   angle_equivalences
            # ...AND we must consider BOTH regular AND auto equivalences.
            # For each combination generate a separate @improper interaction.
            # (I fear this will make the resulting .LT file large.)

            # Use different auto equivalence lookup tables for different
            # atoms in the interaction. (ie the "center" and "end" atoms)
                           
            auto_improper2atom = [auto_improperend2atom,
                                  auto_impropercenter2atom,
                                  auto_improperend2atom,
                                  auto_improperend2atom]

            for i in range(0, 4):
                improper_atom_name = atom_names[i]
                sys.stderr.write('DEBUG: improper_atom_name = '+improper_atom_name+'\n')
                if not is_auto:
                    assert(improper_atom_name[-1] != '_')
                    # assume regular equivalences when looking up atom types
                    sys.stderr.write('DEBUG: equiv_improper2atom['+improper_atom_name+'] = '+
                                     str(equiv_improper2atom[improper_atom_name])+'\n')
                    for a in equiv_improper2atom[improper_atom_name]:
                        atom_combos[i].add(a)
                else:
                    assert((improper_atom_name[-1] == '_') or (improper_atom_name[0] == 'X'))
                    # assume "auto" equivalences when looking up atom types
                    sys.stderr.write('DEBUG: auto_improper2atom['+str(i)+']['+improper_atom_name+'] = \n'
                                     '       '+str(auto_improper2atom[i][improper_atom_name])+'\n')
                    for a in auto_improper2atom[i][improper_atom_name]:
                        atom_combos[i].add(a)

            is_auto = IsAutoInteraction(imp_name_orig)   # is this an "auto" interaction?

            atom_names = ExtractANames(imp_name_orig)            # names of all 4 atoms
            lnames = [atom_names[0], atom_names[2], atom_names[3]]  # names of "leaf" atoms

            #M1     = improper2cross[imp_name_orig][ 2 ]
            #M2     = improper2cross[imp_name_orig][ 0 ]
            #M3     = improper2cross[imp_name_orig][ 3 ]

            #try:
            M1 = improper2cross[imp_name_orig][ImCrossTermID([atom_names[0],
                                                              atom_names[1],
                                                              atom_names[2],
                                                              atom_names[3]])]
            #except KeyError:
            #    M1 = '0.0'

            #try:
            M2 = improper2cross[imp_name_orig][ImCrossTermID([atom_names[2],
                                                              atom_names[1],
                                                              atom_names[0],
                                                              atom_names[3]])]
            #except KeyError:
            #    M2 = '0.0'

            #try:
            M3 = improper2cross[imp_name_orig][ImCrossTermID([atom_names[0],
                                                              atom_names[1],
                                                              atom_names[3],
                                                              atom_names[2]])]
            #except KeyError:
            #    M3 = '0.0'






            # ###### Symmetry: ######
            # Unfortunately, it's time to wade into the messy issue of symmetry.
            #    We desire a way to detect whether an improper interaction
            # between 4 atoms is invariant with respect to atom reordering
            # of the 3 peripheral "leaf" atoms which surround the central atom.
            # In principle, any rearrangement of atoms would require a separate
            # class2 improper interaction.  However, in some cases, when the
            # parameters for these rearrangements are symmetric, we can detect
            # that and warn moltemplate that it is not necessary to generate new
            # improper interactions for every conceivable permutation of these
            # atoms.  Figuring out when it is safe to do that is a headache.
            #   (...but it's necessary.  Otherwise each junction in the molecule
            #   will generate 3*2*1=6 improper interactions which are usually
            #   redundant.  This will slow down the simulation significantly
            #   and may make it difficult to compare the resulting LAMMPS 
            #   input files with those generated by other tools like msi2lmp.)
            #
            # To make this easier, I store the parameters in arrays which 
            # are arranged in a more symmetric way
            M = [0.0, 0.0, 0.0]
            theta0 = [0.0, 0.0, 0.0]
            # noti3[i] = the sorted tuple of integers from the 
            #            set {0,1,2} which remain after deleting i
            noti3 = ((1,2), (0,2), (0,1))
            i_neigh = [ ([0,2,3][ noti3[i][0] ],   # neighbor leaves of ith leaf
                         [0,2,3][ noti3[i][1] ]) for i in range(0,3)]
            for i in range(0, 3):
                # You will notice the pattern "[0,2,3][i]" appears often in the
                # code below because for class 2 force-fields, the second atom
                # (with index 1) is the central atom ("hub" atom), and the three
                # that surround it ("leaf" atoms) have indices 0,2,3.  I want
                # to skip over the central atoms and loop over the leaf atoms
                imTermID = ImCrossTermID([atom_names[ i_neigh[i][0] ],
                                          atom_names[ 1 ],
                                          atom_names[ [0,2,3][i] ],
                                          atom_names[ i_neigh[i][1] ]])
                M[i] = float(improper2cross[imp_name_orig][imTermID])
                ##i_leaf = [0,2,3][i]
                ##M[i] = float(improper2cross[imp_name_orig][ i_leaf ])
                #angle_name_l = SortByEnds([atom_names[i_neigh[i][0]],
                #                           atom_names[ 1 ],
                #                           atom_names[i_neigh[i][1]]])
                #angle_name = EncodeInteractionName(angle_name_l, is_auto)
                #theta0[i] = float(angle2theta0_or[angle_name])

            for i in range(0, 3):
                if (M[ noti3[i][0] ] == M[ noti3[i][1] ]):
                    #and (theta0[ noti3[i][0] ] == theta0[ noti3[i][1] ])):
                    # Then it is safe to swap the order of these two atoms in
                    # the list of atoms when looking up force-field parameters
                    improper2sym[imp_name_orig].add(i_neigh[i][0])
                    improper2sym[imp_name_orig].add(i_neigh[i][1])
                    # Later, I can use these to decide whether or not I need to
                    # change the default script with symmetry rules. (I'm hoping
                    # that "cenJsortIKL.py" should work in most cases.)
                    # CONTINUEHERE: FIGURE OUT WHETHER TO WORRY ABOUT improper2sym
                else:
                    if atom_names[i_neigh[i][0]] == atom_names[i_neigh[i][1]]:
                        raise InputError('Error: Unsupported improper interaction: \"@improper:'+str(imp_name_orig)+'\"\n'
                                         '       This interaction has matching atom aliases:\n'
                                         '       (@atom:'+str(atom_names[i_neigh[i][0]])+
                                         ', @atom:'+str(atom_names[i_neigh[i][1]])+')\n'
                                         '       and yet it lacks symmetry in the corresponding force field parameters.\n'
                                         '       (If this is not a mistake in the .frc file, then ask andrew to\n'
                                         '       fix this limitation.)\n')


            found_at_least_one = False
            for a1 in sorted(list(atom_combos[0])):
                for a2 in sorted(list(atom_combos[1])):
                    sys.stderr.write('DEBUG: improper '+imp_name_orig+' substitutions: '+a1+','+a2+',...\n')
                    for a3 in sorted(list(atom_combos[2])):
                        #(Note: sorting "atom_combos" makes it faster and easier
                        # to follow the loop's progress. This nested loop can be very slow.)
                        theta0s = ['0.0', '0.0', '0.0']
                        aatoms = [['', '',''], ['', '',''], ['', '', '']]
                        #were "auto" equivalences needed to lookup the bond-angle?
                        a_is_auto = [False, False, False]
                        # Collect information from the different terms in a class2 improper:
                        # http://lammps.sandia.gov/doc/improper_class2.html

                        # Loop over the neighbors of the central atom in each improper
                        # interaction and collect all the Mi and Ti parameters. Collect 
                        # them in the order they appear in the formula for the Eaa
                        # term as it appears in the documentation for improper_style class2:
                        # 
                        #    http://lammps.sandia.gov/doc/improper_class2.html
                        #
                        # Eaa = M1 (Tijk - T0)(Tkjl - T2) +   #common leaf node: k (index 2)
                        #       M2 (Tijk - T0)(Tijl - T1) +   #common leaf node: i (index 0)
                        #       M3 (Tijl - T1)(Tkjl - T2)     #common leaf node: l (index 3)
                        # (I'm trying to match the variable names used in this web page
                        #  I wish the author had chosen the M1,M2,M3, T1,T2,T3 order in more
                        #  symmetric way, or at least in a way that makes more sense to me.)

                        #angle_name_l = SortByEnds([atom_names[0], atom_names[1], atom_names[2]])
                        #angle_name = EncodeInteractionName(angle_name_l, is_auto)
                        #theta01 = angle2theta0_or[angle_name]
                        angle_data = LookupBondAngle(a1, a2, a3,
                                                     atom2equiv_angle,
                                                     angle2theta0_or,
                                                     [atom2auto_improperend,
                                                      atom2auto_impropercenter,
                                                      atom2auto_improperend],
                                                     angle2theta0_auto_or)
                        if angle_data == None:
                            # Save time by only continuing if an angle was
                                # found between a1, a2, a3
                            continue
                        theta0s[0], aatoms[0], a_is_auto[0] = angle_data


                        for a4 in sorted(list(atom_combos[3])):
                            theta0s[1] = theta0s[2] = '0.0'
                            aatoms[1] = aatoms[2] = ['', '','']

                            #angle_name_l = SortByEnds(aatoms[0])
                            #angle_name = EncodeInteractionName(angle_name_l[0], is_auto)

                            #theta02 = angle2theta0_or[angle_name]
                            angle_data = LookupBondAngle(a1, a2, a4,
                                                         atom2equiv_angle,
                                                         angle2theta0_or,
                                                         [atom2auto_improperend,
                                                          atom2auto_impropercenter,
                                                          atom2auto_improperend],
                                                         angle2theta0_auto_or)
                            if angle_data == None:
                                # Save time by only continuing if an angle was
                                # found between a1, a2, a4
                                continue
                            theta0s[1], aatoms[1], a_is_auto[1] = angle_data

                            #angle_name_l = SortByEnds(aatoms[1])
                            #angle_name = EncodeInteractionName(angle_name_l, is_auto)


                            #theta03 = angle2theta0_or[angle_name]
                            angle_data = LookupBondAngle(a3, a2, a4,
                                                         atom2equiv_angle,
                                                         angle2theta0_or,
                                                         [atom2auto_improperend,
                                                          atom2auto_impropercenter,
                                                          atom2auto_improperend],
                                                         angle2theta0_auto_or)
                            if angle_data == None:
                                # Save time by only continuing if an angle was
                                # found between a3, a2, a4
                                continue
                            theta0s[2], aatoms[2], a_is_auto[2] = angle_data


                            # The following asserts checks that the two theta0s
                            # are defined whenever the corresponding M is defined.
                            # (Note: The order is LAMMPS-implementation specific.
                            #  See http://lammps.sandia.gov/doc/improper_class2.html)
                            assert((float(theta0s[0]) != 0) or (float(M1) == 0))
                            assert((float(theta0s[2]) != 0) or (float(M1) == 0))
                            assert((float(theta0s[0]) != 0) or (float(M2) == 0))
                            assert((float(theta0s[1]) != 0) or (float(M2) == 0))
                            assert((float(theta0s[1]) != 0) or (float(M3) == 0))
                            assert((float(theta0s[2]) != 0) or (float(M3) == 0))

                            #angle_name_l = SortByEnds(aatoms[2])
                            #angle_name = EncodeInteractionName(angle_name_l, is_auto)


                            imp_name_full = (imp_name_orig + ',' + 
                                             EncodeInteractionName(aatoms[0], a_is_auto[0]) + ',' +
                                             EncodeInteractionName(aatoms[1], a_is_auto[1]) + ',' +
                                             EncodeInteractionName(aatoms[2], a_is_auto[2]))

                            #if imp_name_orig in improper2params_or[imsym][imp_name_orig]:
                            improper2params[imsym][imp_name_full] = ' '.join(improper2params_or[imsym][imp_name_orig])
                            #else:
                            #    improper2params[imsym][imp_name_full] = '0.0 0.0'

                            #if imp_name_orig in improper2cross:
                            improper2class2_aa[imsym][imp_name_full] = \
                                (str(M1)+' '+str(M2)+' '+str(M3)+' '+
                                 str(theta0s[0])+' '+str(theta0s[1])+' '+str(theta0s[2]))
                            #else:
                            #    improper2class2_aa[imsym][imp_name_full] = '0.0 0.0 0.0 0.0 0.0 0.0'
                            #    improper2ver_aa_or[imsym][imp_name_orig] = improper2ver_or[imsym][imp_name_orig]
                            #    improper2ref_aa_or[imsym][imp_name_orig] = improper2ref_or[imsym][imp_name_orig]

                        improper2priority_aa = \
                            DetermineNumericPriority(is_auto,
                                                     aatoms[0] + aatoms[1] + aatoms[2],
                                                     float(improper2ver_aa_or[imsym][imp_name_orig]))
                        improper2ver_aa[imsym][imp_name_full] = improper2ver_aa_or[imsym][imp_name_orig]
                        improper2ref_aa[imsym][imp_name_full] = improper2ref_aa_or[imsym][imp_name_orig]


                        version = max((improper2ver_or[imsym][imp_name_orig],
                                       improper2ver_aa_or[imsym][imp_name_orig]))
                        improper2style[imsym][imp_name_full] = 'class2'
                        improper2ref[imsym][imp_name_full] = improper2ref_or[imsym][imp_name_orig]
                        improper2ver[imsym][imp_name_full] = version
                        improper2priority[imsym][imp_name_full] = \
                            (is_auto,
                             improper_is_secondary_or[imsym][imp_name_orig],
                             improper2priority_or[imsym][imp_name_orig],
                             improper2priority_aa)

                        if len(improper2params) > num_impropers:
                            sys.stderr.write('DEBUG: improper['+imp_name_full+']:\n'
                                             'theta0 = ('
                                             +theta0s[0]+','+theta0s[1]+','+theta0s[2]+')\n')
                            sys.stderr.write('DEBUG: num_impropers = len(improper2params) = '
                                             +str(len(improper2params))+'\n')
                            num_impropers = len(improper2params)


                        found_at_least_one = True


            if not found_at_least_one:
                lines_warnings.append('# WARNING: Undefined rest angle (theta0) in improper: ' +
                                      #'#          the improper interaction between: ' +
                                      ' '.join(atom_names)+'\n')
                # Then we were unable to define cross terms for this interaction because
                # at least one of the equilibrium rest angles could not be determined.
                # This usually occurs because most of the .FRC files which are
                # in circulation are incomplete.  We have to handle this gracefully.
                imp_name_full = (imp_name_orig + ',X,X,X,X,X,X,X,X,X')
                reference = improper2ref_or[imsym][imp_name_orig]
                version = improper2ver_or[imsym][imp_name_orig]
                improper2ref[imsym][imp_name_full] = reference
                improper2ver[imsym][imp_name_full] = version
                improper2params[imsym][imp_name_full] = ' '.join(improper2params_or[imp_name_orig])
                CONTINUEHERE
                improper2style[imp_name_full] = 'class2'
                improper2priority[imp_name_full] = improper2priority_or[imp_name_orig]
                # substitute zeros for the cross term interactions
                improper2class2_aa[imp_name_full] = '0.0 0.0 0.0 120.0 120.0 120.0'
                improper2ref_aa[imp_name_full] = reference
                improper2ver_aa[imp_name_full] = version





        sys.stderr.write("done\n")
        sys.stderr.write("Converting to moltemplate format...\n")





        ##################### BEGIN WRITING FILE #####################





        sys.stdout.write("# This file was generated automatically using:\n")
        sys.stdout.write("# " + g_program_name + " " + " ".join(sys.argv[1:]) + "\n")
        sys.stdout.write("\n\n")
        sys.stdout.write(ffname + " {\n\n")
        
        sys.stdout.write("\n"
                         "  #        AtomType    Mass     # \"Description\" (version, reference)\n\n")
        sys.stdout.write("  write_once(\"Data Masses\") {\n")
        for atype in atom2mass:
            sys.stdout.write("    @atom:" + atype + " " + str(atom2mass[atype]))
            sys.stdout.write("  # ")
            if atype in atom2element:
                sys.stdout.write(atom2element[atype] + ", ")
            #sys.stdout.write(atom2descr[atype])
            sys.stdout.write("\"" + atom2descr[atype] + "\"")
            sys.stdout.write(" (")
            if atype in atom2numbonds:
                sys.stdout.write("nbonds="+str(atom2numbonds[atype])+", ")
            sys.stdout.write("ver=" + atom2ver[atype] +
                             ", ref=" + atom2ref[atype])
            sys.stdout.write(")\n")
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
            #ffid = atype + "_ffid" + atom2ffid[atype]
            sys.stdout.write("  replace{ @atom:" + atype +
                             " @atom:" + atom2ffid[atype] + " }\n")
        
        sys.stdout.write("\n\n\n\n")
        
        
        sys.stdout.write("  # --------------- Non-Bonded Interactions: ---------------------\n"
                         "  # Syntax:\n"
                         "  # pair_coeff    AtomType1    AtomType2   pair_style_name  parameters...\n\n")
        
        sys.stdout.write("  write_once(\"In Settings\") {\n")
                        
        for atype in pair2params:
            assert(atype in pair2style)
            if IsAutoInteraction(bond_name):
                assert(atype in atom2auto_pair)
                if include_auto_equivalences:
                    sys.stdout.write('    pair_coeff @atom:*,ap' + atom2auto_pair[atype] +
                                     ',aq*,ab*,aae*,aac*,ade*,adc*,aie*,aic*' +
                                     ' @atom:*,ap' + atom2auto_pair[atype] +
                                     ',aq*,ab*,aae*,aac*,ade*,adc*,aie*,aic*  ' +
                                     pair2style[atype] + ' ' +
                                     pair2params[atype] +
                                     '  # (ver=' + pair2ver[atype] +
                                     ', ref=' +pair2ref[atype] + ')\n')
                else:
                    continue
            else:
                assert(atype in atom2equiv_pair)
                sys.stdout.write('    pair_coeff ' +
                                 '@atom:*,p' + atom2equiv_pair[atype] + ',b*,a*,d*,i* ' + 
                                 '@atom:*,p' + atom2equiv_pair[atype] + ',b*,a*,d*,i*  ' + 
                                 pair2style[atype] + '  ' +
                                 pair2params[atype] +
                                 '  # (ver=' + pair2ver[atype] +
                                 ', ref=' +pair2ref[atype] + ')\n')
        sys.stdout.write("  } #(end of pair_coeffs)\n\n\n\n")
        






        ################# Print Charge By Bond Interactions ##################
        charge_pair_priority_high_to_low = [x[0] for x in
                                            sorted([x for x in reversed(charge_pair_priority.items())],
                                                   key=itemgetter(1),
                                                   reverse=True)]

        if len(charge_pair_priority) > 0:
            sys.stdout.write("  # ---------- Charge By Bond (a.k.a. \"bond equivalences\") ----------\n")
            # Print rules for generating (2-body) "bond" interactions:
            sys.stdout.write('\n\n\n'
                             '  write_once("Data Charge By Bond") {\n')
            for bond_name in charge_pair_priority_high_to_low:
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(bond_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(bond_name):
                    if include_auto_equivalences:
                        sys.stdout.write('    @atom:*,ap*,aq' + anames[0] +
                                         ',ab*,aae*,aac*,ade*,adc*,aie*,aic*' +
                                         ' @atom:*,ap*,aq' + anames[1] +
                                         ',ab*,aae*,aac*,ade*,adc*,aie*,aic*' +
                                         ' ' + bond2chargepair[bond_name] +
                                         "  # (ver=" + charge_pair_ver[bond_name] +
                                         ", ref=" + charge_pair_ref[bond_name] + ")\n")
                    else:
                        continue
                else:
                    sys.stdout.write('    @atom:*,p*,b' + anames[0] + ',a*,d*,i* ' +
                                     ' @atom:*,p*,b' + anames[1] + ',a*,d*,i* ' +
                                     ' ' + bond2chargepair[bond_name] +
                                     "  # (ver=" + charge_pair_ver[bond_name] +
                                     ", ref=" + charge_pair_ref[bond_name] + ")\n")
            sys.stdout.write('  } #(end of Charge by Bond (bond equivalences))\n\n'
                             '\n\n\n\n')







        ################# Print 2-body Bond Interactions ##################

        bond_names_priority_high_to_low = [x[0] for x in
                                           sorted([x for x in reversed(bond2priority.items())],
                                                  key=itemgetter(1),
                                                  reverse=True)]

        if len(bond2priority) > 0:
            sys.stdout.write("  # --------------- Bond Interactions: ---------------------\n")
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (2-body) "bond" interactions: --\n'
                             '  #  BondType  AtomType1  AtomType2\n')
            sys.stdout.write('\n'
                             '  write_once("Data Bonds By Type')
            if bond_symmetry_subgraph != '':
                sys.stdout.write(' ('+bond_symmetry_subgraph+')')
            sys.stdout.write('") {\n')
            for bond_name in bond_names_priority_high_to_low:
                if not (bond2style[bond_name] in
                        bond_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(bond_name)]
                # Did the user ask us to include "auto" interactions?
                if IsAutoInteraction(bond_name):
                    if include_auto_equivalences:
                        sys.stdout.write('    @bond:' + bond_name + ' ' +
                                         ' @atom:*,ap*,aq*,ab' + anames[0] +
                                         ',aae*,aac*,ade*,adc*,aie*,aic*' +
                                         ' @atom:*,ap*,aq*,ab' + anames[1] +
                                         ',aae*,aac*,ade*,adc*,aie*,aic*' +
                                         '\n')
                    else:
                        continue
                else:
                    sys.stdout.write('    @bond:' + bond_name + ' ' +
                                     ' @atom:*,b' + anames[0] + ',a*,d*,i* ' +
                                     ' @atom:*,b' + anames[1] + ',a*,d*,i* ' +
                                     '\n')

            sys.stdout.write('  }  # end of "Data Bonds By Type" section\n'
                             '\n')

            # Print the force-field parameters for these bond interactions:
            sys.stdout.write('\n\n'
                             '  # ------------ Bond Parameters: ----------\n')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for bond_style in bond_styles:
                if not (bond_style in bond_styles_selected):
                    continue
                sys.stdout.write('    # '+bond_style2docs[bond_style]+'\n')
            sys.stdout.write('\n'
                             '  # Syntax:  \n'
                             '  # bond_coeff BondTypeName  BondStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for bond_name in bond_names_priority_high_to_low:
                if not (bond2style[bond_name] in
                        bond_styles_selected):
                    continue
                # Did the user ask us to include "auto" interactions?
                if (IsAutoInteraction(bond_name) and
                    (not include_auto_equivalences)):
                    continue
                sys.stdout.write('    bond_coeff @bond:'+bond_name+'  '+
                                 bond2style[bond_name] + ' ' +
                                 bond2params[bond_name] +
                                 "  # (ver=" + bond2ver[bond_name] +
                                 ", ref=" +bond2ref[bond_name] + ")\n")

            sys.stdout.write('  }  # end of bond_coeff commands\n'
                             '\n\n')






        ################# Print 3-body Angle Interactions ##################

        ang_names_priority_high_to_low = [x[0] for x in
                                          sorted([x for x in reversed(angle2priority.items())],
                                                 key=itemgetter(1),
                                                 reverse=True)]

        ang_name_abbr = {}            #optional abbreviated name for each interaction
        ang_name_abbr_used = set([])  #make sure we don't reuse these abbreviated names

        if len(angle2priority) > 0:
            sys.stdout.write("  # --------------- Angle Interactions: ---------------------\n")
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (3-body) "angle" interactions: --\n'
                             '  #  AngleType AtomType1 AtomType2 AtomType3  [BondType1 BondType2]\n')
            sys.stdout.write('\n'
                             '  write_once("Data Angles By Type')
            if angle_symmetry_subgraph != '':
                sys.stdout.write(' ('+angle_symmetry_subgraph+')')
            sys.stdout.write('") {\n')
            for angle_name in ang_names_priority_high_to_low:
                if not (angle2style[angle_name] in
                        angle_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(angle_name)]

                angle_is_auto = IsAutoInteraction(angle_name)
                if angle2style[angle_name] == 'class2':
                    anm = [a for a in map(DecodeAName, anames)]
                    bnames = [[a for a in map(DecodeAName, anames[3:5])],
                              [a for a in map(DecodeAName, anames[5:7])]]
                    bond_is_auto1 = IsAutoInteraction(anames[3])
                    bond_is_auto2 = IsAutoInteraction(anames[5])

                if ((angle_is_auto or bond_is_auto1 or bond_is_auto2) and
                    (not include_auto_equivalences)):
                    continue

                # Can we ignore "auto" interactions?
                # (If so, life is much easier)
                if not (angle_is_auto or bond_is_auto1 or bond_is_auto2):
                    if angle2style[angle_name] == 'class2':
                        assert(bnames[0][1] == bnames[1][0])
                        # Optional: Shorten the angle name since some of the atom's bond names are redundant:
                        ang_name_abbr[angle_name] = EncodeInteractionName(map(EncodeAName,
                                                                              anm[0:3] +
                                                                              #[anm[3],anm[4],anm[6]],
                                                                              [bnames[0][0],bnames[0][1],bnames[1][1]]),
                                                                              angle_is_auto)
                        sys.stdout.write('    @angle:' + ang_name_abbr[angle_name] + ' ' +
                                         ' @atom:*,p*,b'+bnames[0][0]+',a'+anames[0]+',d*,i* ' +
                                         ' @atom:*,p*,b'+bnames[0][1]+',a'+anames[1]+',d*,i* ' +
                                         ' @atom:*,p*,b'+bnames[1][1]+',a'+anames[2]+',d*,i*'
                                         '\n')
                    else:
                        ang_name_abbr[angle_name] = angle_name
                        sys.stdout.write('    @angle:' + ang_name_abbr[angle_name] + ' ' +
                                         ' @atom:*,p*,b*,a'+anames[0]+',d*,i* ' +
                                         ' @atom:*,p*,b*,a'+anames[1]+',d*,i* ' +
                                         ' @atom:*,p*,b*,a'+anames[2]+',d*,i*'
                                         '\n')
                else:
                    # Consider "auto" interactions and "auto" atom equivalences
                    ang_name_abbr[angle_name] = angle_name  #(full name)
                    sys.stdout.write('    @angle:' + ang_name_abbr[angle_name] + ' ')

                    if angle2style[angle_name] == 'class2':

                        bshared = 'b*'    #(default. overidden below)
                        abshared = 'ab*'  #(default. overidden below)

                        if angle_is_auto:
                            a1 = a2 = a3 = 'a*'                #Then, dont use regular equivalences for these atoms.
                            aa1 = 'aae' + anames[0] + ',aac*'  #Instead use the corresponding "auto" equivalence names
                            aa2 = 'aae*,aac*' + anames[1]      #for these atoms. (There are different auto equivalence names depending
                            aa3 = 'aae' + anames[2] + ',aac*'  #on if the atom appears in the center (c) or end(e) of the 3-body angle)
                        else:
                            a1 = 'a' + anames[0]               #In this case, use use (regular) equivalence names
                            a2 = 'a' + anames[1]               #for these atoms
                            a3 = 'a' + anames[2]
                            aa1 = aa2 = aa3 = 'aae*,aac*'

                        if not bond_is_auto1:
                            b11 = 'b' + bnames[0][0]     #(bond atom equivalent name)
                            b12 = 'b' + bnames[0][1]     #(bond atom equivalent name)
                            bshared = 'b' + bnames[0][1] #(bond atom equivalent name)
                            ab11 = ab12 = 'ab*'
                        else:
                            b11 = b12 = 'b*'
                            ab11 = 'ab' + bnames[0][0]     #(auto bond atom name)
                            ab12 = 'ab' + bnames[0][1]     #(auto bond atom name)
                            abshared = 'ab' + bnames[0][1] #(auto bond atom name)
                        # print atom 1 information:
                        sys.stdout.write(' @atom:*,p*,'+b11+','+a1+',d*,i*,' +
                                         'ap*,aq*,'+ab11+','+aa1+
                                         ',ade*,adc*,aie*,aic*')
                        if not bond_is_auto2:
                            b21 = 'b' + bnames[1][0]  #(bond atom equivalent name)
                            b22 = 'b' + bnames[1][1]  #(bond atom equivalent name)
                            assert((bshared == 'b*') or (bshared == 'b' + bnames[1][0]))
                            bshared = 'b' + bnames[1][0]
                            ab21 = ab22 = 'ab*'
                        else:
                            b21 = b22 = 'b*'
                            ab21 = 'ab' + bnames[1][0]  #(auto bond atom name)
                            ab22 = 'ab' + bnames[1][1]  #(auto bond atom name)
                            assert((abshared == 'ab*') or (abshared == 'ab' + bnames[1][0]))
                            abshared = 'ab' + bnames[1][0]
                        # print atom 2 information:
                        sys.stdout.write(' @atom:*,p*,'+bshared+','+a2+',d*,i*,' +
                                         'ap*,aq*,'+abshared+','+aa2+
                                         ',ade*,adc*,aie*,aic*')
                        # print atom 3 information:
                        sys.stdout.write(' @atom:*,p*,'+b22+','+a3+',d*,i*,' +
                                         'ap*,aq*,'+ab22+','+aa3+
                                         ',ade*,adc*,aie*,aic*')
                        sys.stdout.write('\n')
                    else:
                        sys.stdout.write('    @angle:' + ang_name_abbr[angle_name] + ' ' +
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae'+anames[0]+'aac*,ade*,adc*,aie*,aic* '
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac'+anames[1]+',ade*,adc*,aie*,aic* '
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae'+anames[2]+'aac*,ade*,adc*,aie*,aic* '
                                         '\n')

                assert(ang_name_abbr[angle_name] not in ang_name_abbr_used)
                ang_name_abbr_used.add(ang_name_abbr[angle_name])

            sys.stdout.write('  }  # end of "Data Angles By Type" section\n'
                             '\n')

            # Print the force-field parameters for these angle interactions:
            sys.stdout.write('\n\n'
                             '  # ------- Angle Force Field Parameters: -------')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for angle_style in angle_styles:
                if not (angle_style in angle_styles_selected):
                    continue
                sys.stdout.write('    # '+angle_style2docs[angle_style]+'\n')
            sys.stdout.write('\n'
                             '  # Syntax:  \n'
                             '  # angle_coeff AngleTypeName  AngleStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for angle_name in ang_names_priority_high_to_low:
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(angle_name)]
                if not (angle2style[angle_name] in
                        angle_styles_selected):
                    continue


                # Did the user ask us to include "auto" interactions?
                #if (IsAutoInteraction(angle_name) and
                #    (not include_auto_equivalences)):
                #    continue
                # the if statement above is covered by the following:
                if angle_name not in ang_name_abbr:
                    continue

                sys.stdout.write('    angle_coeff @angle:'+ang_name_abbr[angle_name]+'  '+
                                 angle2style[angle_name] + ' ' +
                                 angle2params[angle_name] + 
                                 "  # (ver=" + angle2ver[angle_name] +
                                 ", ref=" + angle2ref[angle_name] + ")\n")
                if angle_name in angle2class2_bb:
                    sys.stdout.write('    angle_coeff @angle:'+ang_name_abbr[angle_name]+'  '+
                                     angle2style[angle_name] + ' bb ' +
                                     angle2class2_bb[angle_name] +
                                     "  # (ver=" + angle2ver_bb[angle_name] +
                                     ", ref=" + angle2ref_bb[angle_name] + ")\n")

                    assert(angle_name in angle2class2_ba)
                    sys.stdout.write('    angle_coeff @angle:'+ang_name_abbr[angle_name]+'  '+
                                     angle2style[angle_name] + ' ba ' +
                                     angle2class2_ba[angle_name] +
                                     "  # (ver=" + angle2ver_ba[angle_name] +
                                     ", ref=" + angle2ref_ba[angle_name] + ")\n")
            sys.stdout.write('  }  # end of angle_coeff commands\n'
                             '\n\n')







        ################# Print 4-body Dihedral Interactions ##################

        dih_names_priority_high_to_low = [x[0] for x in
                                          sorted([x for x in reversed(dihedral2priority.items())],
                                                 key=itemgetter(1),
                                                 reverse=True)]

        dih_name_abbr = {}            #optional abbreviated name for each interaction
        dih_name_abbr_used = set([])  #make sure we don't reuse these abbreviated names

        if len(dih_names_priority_high_to_low) > 0:
            sys.stdout.write('  # --------------- Dihedral Interactions: ---------------------\n')
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (4-body) "dihedral" interactions: --\n'
                             '  #  DihedralType AtmType1 AtmType2 AtmType3 AtmType3 [BondType1 Bnd2 Bnd3]\n')
            sys.stdout.write('\n\n'
                             '  write_once("Data Dihedrals By Type')
            if dihedral_symmetry_subgraph != '':
                sys.stdout.write(' ('+dihedral_symmetry_subgraph+')')
            sys.stdout.write('") {\n')



            for dihedral_name in dih_names_priority_high_to_low:
                if not (dihedral2style[dihedral_name] in
                        dihedral_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(dihedral_name)]

                dihedral_is_auto = IsAutoInteraction(dihedral_name)
                if dihedral2style[dihedral_name] == 'class2':
                    anm = [a for a in map(DecodeAName, anames)]
                    bnames = [[a for a in map(DecodeAName, anames[4:6])],
                              [a for a in map(DecodeAName, anames[6:8])],
                              [a for a in map(DecodeAName, anames[8:10])]]
                    bond_is_auto1 = IsAutoInteraction(anames[4])
                    bond_is_auto2 = IsAutoInteraction(anames[6])
                    bond_is_auto3 = IsAutoInteraction(anames[8])
                    ang_names = [[a for a in map(DecodeAName, anames[10:13])],
                                 [a for a in map(DecodeAName, anames[13:16])]]
                    angle_is_auto1 = IsAutoInteraction(anames[10])
                    angle_is_auto2 = IsAutoInteraction(anames[13])

                if ((dihedral_is_auto or
                     angle_is_auto1 or angle_is_auto2 or 
                     bond_is_auto1 or bond_is_auto2 or bond_is_auto3) and
                    (not include_auto_equivalences)):
                    continue

                # Can we ignore "auto" interactions?
                # (If so, life is much easier)
                if not (dihedral_is_auto or
                        angle_is_auto1 or angle_is_auto2 or 
                        bond_is_auto1 or bond_is_auto2 or bond_is_auto3):

                    if dihedral2style[dihedral_name] == 'class2':
                        assert(bnames[0][1] == bnames[1][0])
                        assert(bnames[1][1] == bnames[2][0])
                        assert(ang_names[0][1] == ang_names[1][0])
                        assert(ang_names[0][2] == ang_names[1][1])

                        # Optional: Shorten the dihedral name since some of the atom's bond names are redundant:
                        dih_name_abbr[dihedral_name] = EncodeInteractionName(map(EncodeAName,
                                                                                 anm[0:4] +
                                                                                 #[bnames[0][0], bnames[0][1],
                                                                                 # bnames[1][1], bnames[2][1]]
                                                                                 [anm[4],anm[5],anm[7],anm[9]]+
                                                                                 #[ang_names[0][0],
                                                                                 # ang_names[0][1],
                                                                                 # ang_names[0][2],
                                                                                 # ang_names[1][2]]
                                                                                 [anm[10],anm[11],anm[12],anm[15]]),
                                                                             is_auto)

                        sys.stdout.write('    @dihedral:' + dih_name_abbr[dihedral_name] + ' ' +
                                         ' @atom:*,p*,b'+bnames[0][0]+',a'+ang_names[0][0]+',d'+anames[0]+',i* ' +
                                         ' @atom:*,p*,b'+bnames[0][1]+',a'+ang_names[0][1]+',d'+anames[1]+',i* ' +
                                         ' @atom:*,p*,b'+bnames[1][1]+',a'+ang_names[0][2]+',d'+anames[2]+',i* '
                                         ' @atom:*,p*,b'+bnames[2][1]+',a'+ang_names[1][2]+',d'+anames[3]+',i*'
                                         '\n')
                    else:
                        dih_name_abbr[dihedral_name] = dihedral_name
                        sys.stdout.write('    @dihedral:' + dih_name_abbr[dihedral_name] + ' ' +
                                         ' @atom:*,p*,b*,a*,d'+anames[0]+',i* ' +
                                         ' @atom:*,p*,b*,a*,d'+anames[1]+',i* ' +
                                         ' @atom:*,p*,b*,a*,d'+anames[2]+',i* '
                                         ' @atom:*,p*,b*,a*,d'+anames[3]+',i*' +
                                         '\n')
                else:
                    # Consider "auto" interactions and "auto" atom equivalences
                    dih_name_abbr[dihedral_name] = dihedral_name  #(full name)
                    sys.stdout.write('    @dihedral:' + dih_name_abbr[dihedral_name] + ' ')

                    if dihedral2style[dihedral_name] == 'class2':

                        # equivalent names of atoms shared by more than one bond:
                        # (names ending in * mean they were unspecified for this 
                        #  dihedral interaction.  By default, this is the case.)
                        bshared1 = 'b*'       #(default. overidden below)
                        bshared2 = 'b*'       #(default. overidden below)
                        abshared1 = 'ab*'     #(default. overidden below)
                        abshared2 = 'ab*'     #(default. overidden below)

                        # equivalent names of atoms shared by more than one angle interaction:
                        # (names ending in * mean they were unspecified for this 
                        #  dihedral interaction.  By default, this is the case.)
                        ashared1 = 'a*'       #(default. overidden below)
                        ashared2 = 'a*'       #(default. overidden below)
                        aac_shared1 = 'aac*'  #(default. overidden below)
                        aae_shared1 = 'aae*'  #(default. overidden below)
                        aac_shared2 = 'aac*'  #(default. overidden below)
                        aae_shared2 = 'aae*'  #(default. overidden below)

                        if dihedral_is_auto:
                            d1 = d2 = d3 = d4 = 'd*'           #Then, dont use regular equivalences for these atoms.
                            ad1 = 'ade' + anames[0] + ',adc*'  #Instead use the corresponding "auto"
                            ad2 = 'ade*,adc*' + anames[1]      #equivalence names for these atoms.
                            ad3 = 'ade*,adc*' + anames[1]      #(There are different auto equivalence names depending upon
                            ad4 = 'ade' + anames[2] + ',adc*'  # if the atom appears in the center (c) or end(e) of the dihedral)
                        else:
                            d1 = 'd' + anames[0]               # In this case, use use (regular) equivalence names
                            d2 = 'd' + anames[1]               # for these atoms
                            d3 = 'd' + anames[2]
                            d4 = 'd' + anames[3]
                            ad1 = ad2 = ad3 = ad4 = 'ade*,adc*'

                        if not bond_is_auto1:
                            b11 = 'b' + bnames[0][0]      #(bond atom equivalent name)
                            b12 = 'b' + bnames[0][1]      #(bond atom equivalent name)
                            bshared1 = 'b' + bnames[0][1] #(bond atom equivalent name)
                            ab11 = ab12 = 'ab*'
                        else:
                            b11 = b12 = 'b*'
                            ab11 = 'ab' + bnames[0][0]      #(auto bond atom name)
                            ab12 = 'ab' + bnames[0][1]      #(auto bond atom name)
                            abshared1 = 'ab' + bnames[0][1] #(auto bond atom name)

                        if not bond_is_auto2:
                            b21 = 'b' + bnames[1][0]      #(bond atom equivalent name)
                            b22 = 'b' + bnames[1][1]      #(bond atom equivalent name)
                            assert((bshared1 == 'b*') or (bshared1 == 'b' + bnames[1][0]))
                            bshared1 = 'b' + bnames[1][0] #(bond atom equivalent name)
                            assert((bshared2 == 'b*') or (bshared2 == 'b' + bnames[1][1]))
                            bshared2 = 'b' + bnames[1][1] #(bond atom equivalent name)
                            ab21 = ab22 = 'ab*'
                        else:
                            b21 = b22 = 'b*'
                            ab21 = 'ab' + bnames[1][0]      #(auto bond atom name)
                            ab22 = 'ab' + bnames[1][1]      #(auto bond atom name)
                            assert((abshared1 == 'ab*') or (abshared1 == 'ab' + bnames[1][0]))
                            abshared1 = 'ab' + bnames[1][0] #(auto bond atom name)
                            assert((abshared2 == 'ab*') or (abshared2 == 'ab' + bnames[1][1]))
                            abshared2 = 'ab' + bnames[1][1] #(auto bond atom name)

                        if not bond_is_auto3:
                            b31 = 'b' + bnames[2][0]      #(bond atom equivalent name)
                            b32 = 'b' + bnames[2][1]      #(bond atom equivalent name)
                            assert((bshared2 == 'b*') or (bshared2 == 'b' + bnames[2][0]))
                            bshared2 = 'b' + bnames[2][0] #(bond atom equivalent name)
                            ab31 = ab32 = 'ab*'
                        else:
                            b31 = b32 = 'b*'
                            ab31 = 'ab' + bnames[2][0]      #(auto bond atom name)
                            ab32 = 'ab' + bnames[2][1]      #(auto bond atom name)
                            assert((abshared2 == 'ab*') or (abshared2 == 'ab' + bnames[2][0]))
                            abshared2 = 'ab' + bnames[2][0] #(auto bond atom name)

                        if not angle_is_auto1:
                            a11 = 'a' + ang_names[0][0]      #(angle atom equivalent name)
                            a12 = 'a' + ang_names[0][1]      #(angle atom equivalent name)
                            a13 = 'a' + ang_names[0][2]      #(angle atom equivalent name)
                            ashared1 = 'a' + ang_names[0][1] #(angle atom equivalent name)
                            ashared2 = 'a' + ang_names[0][2] #(angle atom equivalent name)
                            aa11 = 'aae*,aac*'
                            aa12 = 'aae*,aac*'
                            aa13 = 'aae*,aac*'
                        else:
                            a11 = a12 = a13 = 'a*'
                            aa11 = 'aae'+ang_names[0][0]+'aac*'  #(auto angle atom name)
                            aa12 = 'aae*,aac'+ang_names[0][1]    #(auto angle atom name)
                            aa13 = 'aae'+ang_names[0][2]+'aac*'  #(auto angle atom name)
                            aac_shared1 = 'aac'+ang_names[0][1]  #(auto angle atom name)
                            aae_shared2 = 'aae'+ang_names[0][2]  #(auto angle atom name)

                        if not angle_is_auto2:
                            a21 = 'a' + ang_names[1][0]      #(angle atom equivalent name)
                            a22 = 'a' + ang_names[1][1]      #(angle atom equivalent name)
                            a23 = 'a' + ang_names[1][2]      #(angle atom equivalent name)
                            assert((ashared1 == 'a*') or (ashared1 == 'a' + ang_names[1][0]))
                            ashared1 = 'a' + ang_names[1][0] #(angle atom equivalent name)
                            assert((ashared2 == 'a*') or (ashared2 == 'a' + ang_names[1][1]))
                            ashared2 = 'a' + ang_names[1][1] #(angle atom equivalent name)
                            aa21 = 'aae*,aac*'
                            aa22 = 'aae*,aac*'
                            aa23 = 'aae*,aac*'
                        else:
                            a21 = a22 = a23 = 'a*'
                            aa21 = 'aae'+ang_names[1][0]+',aac*'  #(auto angle atom name)
                            aa22 = 'aae*,aac'+ang_names[1][1]     #(auto angle atom name)
                            aa23 = 'aae'+ang_names[1][2]+',aac*'  #(auto angle atom name)
                            aae_shared1 = 'aae'+ang_names[1][0]   #(auto angle atom name)
                            aac_shared2 = 'aac'+ang_names[1][1]   #(auto angle atom name)


                        # print atom 1 information:
                        sys.stdout.write(' @atom:*,p*,'+b11+','+a11+','+d1+',i*,' +
                                         'ap*,aq*,'+ab11+','+aa11+',' +
                                         ad1+',aie*,aic*')
                        # print atom 2 information:
                        sys.stdout.write(' @atom:*,p*,'+bshared1+','+ashared1+','+d2+',i*,' +
                                         'ap*,aq*,'+abshared1+','+aae_shared1+','+aac_shared1+',' +
                                         ad2+',aie*,aic*')
                        # print atom 3 information:
                        sys.stdout.write(' @atom:*,p*,'+bshared2+','+ashared2+','+d3+',i*,' +
                                         'ap*,aq*,'+abshared2+','+aae_shared2+','+aac_shared2+',' +
                                         ad3+',aie*,aic*')
                        # print atom 4 information:
                        sys.stdout.write(' @atom:*,p*,'+b32+','+a23+','+d4+',i*,' +
                                         'ap*,aq*,'+ab32+','+aa23+',' +
                                         ad4+',aie*,aic*')
                        sys.stdout.write('\n')
                    else:
                        assert(dihedral_is_auto)  #(so we should use "auto" equivalence names for these atoms)
                        sys.stdout.write('    @dihedral:' + dih_name_abbr[dihedral_name] + ' ' +
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade'+anames[0]+',adc*,aie*,aic* '
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade*,adc'+anames[1]+',aie*,aic* '
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade*,adc'+anames[2]+',aie*,aic* '
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade'+anames[3]+',adc*,aie*,aic* '
                                         '\n')




                assert(dih_name_abbr[dihedral_name] not in dih_name_abbr_used)
                dih_name_abbr_used.add(dih_name_abbr[dihedral_name])

            sys.stdout.write('  }  # end of "Data Dihedrals By Type" section\n'
                             '\n')

            # Print the force-field parameters for these dihedral interactions:
            sys.stdout.write('\n\n'
                             '  # ------- Dihedral Force Field Parameters: -------\n')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for dihedral_style in dihedral_styles:
                if not (dihedral_style in dihedral_styles_selected):
                    continue
                sys.stdout.write('    # '+dihedral_style2docs[dihedral_style]+'\n')
            sys.stdout.write('\n'
                             '  # Syntax:  \n'
                             '  # dihedral_coeff DihedralTypeName  DihedralStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for dihedral_name in dih_names_priority_high_to_low:
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(dihedral_name)]
                #if (len(anames) == 4) and dihedral2style[dihedral_name] == 'class2':
                #    continue
                if not (dihedral2style[dihedral_name] in
                        dihedral_styles_selected):
                    continue

                # Did the user ask us to include "auto" interactions?
                #if (IsAutoInteraction(dihedral_name) and
                #    (not include_auto_equivalences)):
                #    continue
                # the if statement above is covered by the following:
                if dihedral_name not in dih_name_abbr:
                    continue

                sys.stdout.write('    dihedral_coeff @dihedral:'+dih_name_abbr[dihedral_name]+'  '+
                                 dihedral2style[dihedral_name] + ' ' +
                                 dihedral2params[dihedral_name] +
                                 "  # (ver=" + dihedral2ver[dihedral_name] +
                                 ", ref=" + dihedral2ref[dihedral_name] + ")\n")
                if dihedral_name in dihedral2class2_mbt:
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dih_name_abbr[dihedral_name]+'  '+
                                     dihedral2style[dihedral_name] + ' mbt ' +
                                     dihedral2class2_mbt[dihedral_name] +
                                     "  # (ver=" + dihedral2ver_mbt[dihedral_name] +
                                     ", ref=" + dihedral2ref_mbt[dihedral_name] + ")\n")

                    assert(dihedral_name in dihedral2class2_ebt)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dih_name_abbr[dihedral_name]+'  '+
                                     dihedral2style[dihedral_name] + ' ebt ' +
                                     dihedral2class2_ebt[dihedral_name] +
                                     "  # (ver=" + dihedral2ver_ebt[dihedral_name] +
                                     ", ref=" + dihedral2ref_ebt[dihedral_name] + ")\n")

                    assert(dihedral_name in dihedral2class2_at)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dih_name_abbr[dihedral_name]+'  '+
                                     dihedral2style[dihedral_name] + ' at ' +
                                     dihedral2class2_at[dihedral_name] +
                                     "  # (ver=" + dihedral2ver_at[dihedral_name] +
                                     ", ref=" + dihedral2ref_at[dihedral_name] + ")\n")

                    assert(dihedral_name in dihedral2class2_aat)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dih_name_abbr[dihedral_name]+'  '+
                                     dihedral2style[dihedral_name] + ' aat ' +
                                     dihedral2class2_aat[dihedral_name] +
                                     "  # (ver=" + dihedral2ver_aat[dihedral_name] +
                                     ", ref=" + dihedral2ref_aat[dihedral_name] + ")\n")
                    assert(dihedral_name in dihedral2class2_bb13)
                    sys.stdout.write('    dihedral_coeff @dihedral:'+dih_name_abbr[dihedral_name]+'  '+
                                     dihedral2style[dihedral_name] + ' bb13 ' +
                                     dihedral2class2_bb13[dihedral_name] +
                                     "  # (ver=" + dihedral2ver_bb13[dihedral_name] +
                                     ", ref=" + dihedral2ref_bb13[dihedral_name] + ")\n")
            sys.stdout.write('  }  # end of dihedral_coeff commands\n'
                             '\n\n')





        ################# Print 4-body Improper Interactions ##################

        imp_names_priority_high_to_low = [x[0] for x in
                                          sorted([x for x in reversed(improper2priority.items())],
                                                 key=itemgetter(1),
                                                 reverse=True)]

        imp_name_abbr = {}            #optional abbreviated name for each interaction
        imp_name_abbr_used = set([])  #make sure we don't reuse these abbreviated names

        if len(imp_names_priority_high_to_low) > 0:
            sys.stdout.write("  # --------------- Improper Interactions: ---------------------\n")
            sys.stdout.write('\n'
                             '\n'
                             '  # -- Rules for generating (4-body) "improper" interactions: --\n'
                             '  #  ImproperType AtmType1 AtmType2 AtmType3 AtmType3 [BondType1 Bnd2 Bnd3]\n')
            sys.stdout.write('\n'
                             '  write_once("Data Impropers By Type')
            if improper_symmetry_subgraph != '':
                sys.stdout.write(' ('+improper_symmetry_subgraph+')')
            sys.stdout.write('") {\n')
            for improper_name in imp_names_priority_high_to_low:
                if not (improper2style[improper_name] in
                        improper_styles_selected):
                    continue
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(improper_name)]
                #if (len(anames) == 4) and improper2style[improper_name] == 'class2':
                #    continue
                ang_names = [[a for a in map(DecodeAName, anames[4:7])],
                             [a for a in map(DecodeAName, anames[7:10])],
                             [a for a in map(DecodeAName, anames[10:13])]]
                anm = [a for a in map(DecodeAName, anames)]

                improper_is_auto = IsAutoInteraction(improper_name)
                if improper2style[improper_name] == 'class2':
                    angle_is_auto1 = IsAutoInteraction(anames[4])
                    angle_is_auto2 = IsAutoInteraction(anames[7])
                    angle_is_auto3 = IsAutoInteraction(anames[10])

                if ((improper_is_auto or
                     angle_is_auto1 or
                     angle_is_auto2 or
                     angle_is_auto3) and
                    (not include_auto_equivalences)):
                    continue

                # Can we ignore "auto" interactions?
                # (If so, life is much easier)
                if not (improper_is_auto or
                        angle_is_auto1 or
                        angle_is_auto2 or
                        angle_is_auto3):
                    if improper2style[improper_name] == 'class2':
                        # NOTE: atom orderings here are LAMMPS implementation specific.
                        # http://lammps.sandia.gov/doc/improper_class2.html
                        #ang_names[0] <==> (a1, a2, a3) <==>  (i, j, k)
                        #ang_names[1] <==> (a1, a2, a4) <==>  (i, j, l)
                        #ang_names[2] <==> (a3, a2, a4) <==>  (k, j, l)
                        assert(ang_names[0][1] == ang_names[1][1] == ang_names[2][1])
                        assert(ang_names[0][0] == ang_names[1][0])
                        assert(ang_names[1][2] == ang_names[2][2])
                        assert(ang_names[2][0] == ang_names[0][2])

                        # Optional: Shorten the improper name since some of the atom's bond names are redundant:
                        imp_name_abbr[improper_name] = EncodeInteractionName(map(EncodeAName, anm[0:4] +
                                                                                 [ang_names[0][0],
                                                                                  ang_names[0][1],
                                                                                  ang_names[0][2],
                                                                                  ang_names[1][2]]),
                                                                                  #[anm[4],anm[5],anm[6],
                                                                                  #anm[9]],
                                                                             improper_is_auto)
                        sys.stdout.write('    @improper:' + imp_name_abbr[improper_name] + ' ' +
                                         ' @atom:*,p*,b*,a'+ang_names[0][0]+',d*,i' + anames[0] +
                                         ' @atom:*,p*,b*,a'+ang_names[0][1]+',d*,i' + anames[1] +
                                         ' @atom:*,p*,b*,a'+ang_names[0][2]+',d*,i' + anames[2] +
                                         ' @atom:*,p*,b*,a'+ang_names[1][2]+',d*,i' + anames[3] +
                                         '\n')
                    else:
                        imp_name_abbr[improper_name] = improper_name
                        sys.stdout.write('    @improper:' + imp_name_abbr[improper_name] + ' ' +
                                         ' @atom:*,p*,b*,a*,d*,i' + anames[0] +
                                         ' @atom:*,p*,b*,a*,d*,i' + anames[1] +
                                         ' @atom:*,p*,b*,a*,d*,i' + anames[2] +
                                         ' @atom:*,p*,b*,a*,d*,i' + anames[3] +
                                         '\n')
                else:
                    # Consider "auto" interactions and "auto" atom equivalences
                    imp_name_abbr[improper_name] = improper_name  #(full name)
                    sys.stdout.write('    @improper:' + imp_name_abbr[improper_name] + ' ')

                    if improper2style[improper_name] == 'class2':

                        #ang_names[0] <==> (a1, a2, a3) <==>  (i, j, k)
                        #ang_names[1] <==> (a1, a2, a4) <==>  (i, j, l)
                        #ang_names[2] <==> (a3, a2, a4) <==>  (k, j, l)

                        # default angle atom equivalence names:
                        ashared1 = 'a*'          #(default for a1 <-> ang_names[0][0], ang_names[1][0])
                        ashared2 = 'a*'          #(default for a2 <-> ang_names[0][1], ang_names[1][1], ang_names[2][1])
                        ashared3 = 'a*'          #(default for a3 <-> ang_names[2][0], ang_names[0][2])
                        ashared4 = 'a*'          #(default for a4 <-> ang_names[1][2], ang_names[2][2])

                        # default auto angle atom equivalence names:
                        aashared1 = 'aae*,aac*'  #(default for a1 <-> ang_names[0][0], ang_names[1][0])
                        aashared2 = 'aae*,aac*'  #(default for a2 <-> ang_names[0][1], ang_names[1][1], ang_names[2][1])
                        aashared3 = 'aae*,aac*'  #(default for a3 <-> ang_names[2][0], ang_names[0][2])
                        aashared4 = 'aae*,aac*'  #(default for a4 <-> ang_names[1][2], ang_names[2][2])

                        if improper_is_auto:
                            i1 = i2 = i3 = i4 = 'i*'           #Then, dont use regular equivalences for these atoms.
                            ai1 = 'aie' + anames[0] + ',aic*'  #Instead use the corresponding "auto" equivalence names
                            ai2 = 'aie*,aic*' + anames[1]      #for these atoms. (There are different auto equivalence names depending
                            ai3 = 'aie' + anames[2] + ',aic*'  #on if the atom appears in the center (c) or end(e)
                            ai4 = 'aie' + anames[3] + ',aic*'
                        else:
                            i1 = 'i' + anames[0]               #In this case, use use (regular) equivalence names
                            i2 = 'i' + anames[1]               #for these atoms
                            i3 = 'i' + anames[2]
                            i4 = 'i' + anames[3]
                            ai1 = ai2 = ai3 = 'aie*,aic*'

                        #For reference, LAMMPS-specific atom ordering:
                        #ang_names[0] <==> (a1, a2, a3) <==>  (i, j, k)
                        #ang_names[1] <==> (a1, a2, a4) <==>  (i, j, l)
                        #ang_names[2] <==> (a3, a2, a4) <==>  (k, j, l)
                        if not angle_is_auto1:
                            ashared1 = 'a' + ang_names[0][0]
                            ashared2 = 'a' + ang_names[0][1]
                            ashared3 = 'a' + ang_names[0][2]
                        else:
                            aashared1 = 'aae' + ang_names[0][0] + ',aac*'
                            aashared2 = 'aae*,aac' + ang_names[0][1]
                            aashared3 = 'aae' + ang_names[0][2] + ',aac*'

                        #For reference, LAMMPS-specific atom ordering:
                        #ang_names[0] <==> (a1, a2, a3) <==>  (i, j, k)
                        #ang_names[1] <==> (a1, a2, a4) <==>  (i, j, l)
                        #ang_names[2] <==> (a3, a2, a4) <==>  (k, j, l)
                        if not angle_is_auto2:
                            assert((ashared1 == 'a*') or (ashared1 == 'a' + ang_names[1][0]))
                            ashared1 = 'a' + ang_names[1][0]
                            assert((ashared2 == 'a*') or (ashared2 == 'a' + ang_names[1][1]))
                            ashared2 = 'a' + ang_names[1][1]
                            ashared4 = 'a' + ang_names[1][2]
                        else:
                            assert((aashared1 == 'aae*,aac*') or (aashared1 == 'aae' + ang_names[1][0] + ',aac*'))
                            aashared1 = 'aae' + ang_names[1][0] + ',aac*'
                            assert((aashared2 == 'aae*,aac*') or (aashared2 == 'aae*,aac' + ang_names[1][1]))
                            aashared2 = 'aae*,aac' + ang_names[1][1]
                            aashared4 = 'aae' + ang_names[1][2] + ',aac*'

                        #For reference, LAMMPS-specific atom ordering:
                        #ang_names[0] <==> (a1, a2, a3) <==>  (i, j, k)
                        #ang_names[1] <==> (a1, a2, a4) <==>  (i, j, l)
                        #ang_names[2] <==> (a3, a2, a4) <==>  (k, j, l)
                        if not angle_is_auto3:
                            assert((ashared3 == 'a*') or (ashared3 == 'a' + ang_names[2][0]))
                            ashared3 = 'a' + ang_names[2][0]
                            assert((ashared2 == 'a*') or (ashared2 == 'a' + ang_names[2][1]))
                            ashared2 = 'a' + ang_names[2][1]
                            assert((ashared4 == 'a*') or (ashared4 == 'a' + ang_names[2][2]))
                            ashared4 = 'a' + ang_names[2][2]
                        else:
                            assert((aashared3 == 'aae*,aac*') or (aashared3 == 'aae' + ang_names[2][0] + ',aac*'))
                            aashared3 = 'aae' + ang_names[2][0] + ',aac*'
                            assert((aashared2 == 'aae*,aac*') or (aashared2 == 'aae*,aac' + ang_names[2][1]))
                            aashared2 = 'aae*,aac' + ang_names[2][1]
                            assert((aashared4 == 'aae*,aac*') or (aashared4 == 'aae' + ang_names[2][2] + ',aac*'))
                            aashared4 = 'aae' + ang_names[2][2] + ',aac*'

                        # print atom 1 information:
                        sys.stdout.write(' @atom:*,p*,b*,'+ashared1+',d*,'+i1+','+
                                         'ap*,aq*,ab*,'+aashared1+',ad*,'+ai1)
                        # print atom 2 information:
                        sys.stdout.write(' @atom:*,p*,b*,'+ashared2+',d*,'+i2+','+
                                         'ap*,aq*,ab*,'+aashared2+',ad*,'+ai2)
                        # print atom 3 information:
                        sys.stdout.write(' @atom:*,p*,b*,'+ashared3+',d*,'+i3+','+
                                         'ap*,aq*,ab*,'+aashared3+',ad*,'+ai3)
                        # print atom 4 information:
                        sys.stdout.write(' @atom:*,p*,b*,'+ashared4+',d*,'+i4+','+
                                         'ap*,aq*,ab*,'+aashared4+',ad*,'+ai4)
                        sys.stdout.write('\n')
                    else:
                        sys.stdout.write('    @improper:' + imp_name_abbr[improper_name] + ' ' +
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade*,adc*,aie*,aie'+anames[0]+',aic*'
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade*,adc*,aie*,aie*,aic'+anames[1]+
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade*,adc*,aie*,aie'+anames[2]+',aic*'
                                         ' @atom:*,p*,b*,d*,i*,' +
                                         'ap*,aq*,ab*,aae*,aac*,ade*,adc*,aie*,aie'+anames[3]+',aic*'
                                         '\n')

                assert(imp_name_abbr[improper_name] not in imp_name_abbr_used)
                imp_name_abbr_used.add(imp_name_abbr[improper_name])




            sys.stdout.write('  }  # end of "Data Impropers By Type" section\n'
                             '\n')

            # Print the force-field parameters for these improper interactions:
            sys.stdout.write('\n\n'
                             '  # ------- Improper Force Field Parameters: -------\n')
            sys.stdout.write('  # For an explanation of these parameters, visit:\n')
            for improper_style in improper_styles:
                if not (improper_style in improper_styles_selected):
                    continue
                sys.stdout.write('    # '+improper_style2docs[improper_style]+'\n')
            sys.stdout.write('\n'
                             '# Syntax:  \n'
                             '  # improper_coeff ImproperTypeName  ImproperStyle  parameters...\n\n')
            sys.stdout.write('\n'
                             '  write_once("In Settings") {\n')
            for improper_name in imp_names_priority_high_to_low:
                anames = ['*' if x=='X' else x
                          for x in ExtractANames(improper_name)]
                #if (len(anames) == 4) and improper2style[improper_name] == 'class2':
                #    continue
                # Optional: Shorten the angle name since some of the bnames are redundant:

                is_auto = IsAutoInteraction(improper_name)

                if not (improper2style[improper_name] in
                        improper_styles_selected):
                    continue

                # Did the user ask us to include "auto" interactions?
                #if (IsAutoInteraction(improper_name) and
                #    (not include_auto_equivalences)):
                #    continue
                # the if statement above is covered by the following:
                if improper_name not in imp_name_abbr:
                    continue

                sys.stdout.write('    improper_coeff @improper:'+imp_name_abbr[improper_name]+'  '+
                                 improper2style[improper_name] + '  ' +
                                 improper2params[improper_name] +
                                 "  # (ver=" + improper2ver[improper_name] +
                                 ", ref=" + improper2ref[improper_name] + ")\n")
                if improper_name in improper2class2_aa:
                    sys.stdout.write('    improper_coeff @improper:'+imp_name_abbr[improper_name]+'  '+
                                     improper2style[improper_name] + ' aa ' +
                                     improper2class2_aa[improper_name] +
                                     "  # (ver=" + improper2ver_aa[improper_name] +
                                     ", ref=" + improper2ref[improper_name] + ")\n")
            sys.stdout.write('  }  # end of improper_coeff commands\n'
                             '\n\n')



        sys.stdout.write('\n\n\n\n'
                         '  # -------------------- Select LAMMPS style(s) ------------------\n'
                         '\n')

        
        sys.stdout.write('\n'
                         '  # LAMMPS supports many different kinds of bonded and non-bonded\n'
                         '  # interactions which can be selected at run time.  Eventually\n'
                         '  # we must inform LAMMPS which of them we will need.  We specify\n'
                         '  # this in the "In Init" section: \n\n')
        
        sys.stdout.write('  write_once("In Init") {\n')
        sys.stdout.write('    units real\n')
        sys.stdout.write('    atom_style full\n')

        if len(bond_styles) > 0:
            sys.stdout.write('    bond_style hybrid')
            for bond_style in bond_styles:
                if not (bond_style in bond_styles_selected):
                    continue
                sys.stdout.write(' ' + bond_style)
            sys.stdout.write('\n')
            for bond_style in bond_styles:
                if not (bond_style in bond_styles_selected):
                    continue
                sys.stdout.write('    # '+bond_style2docs[bond_style]+'\n')
            sys.stdout.write('\n')

        if len(angle_styles) > 0:
            sys.stdout.write('    angle_style hybrid')
            for angle_style in angle_styles:
                if not (angle_style in angle_styles_selected):
                    continue
                sys.stdout.write(' ' + angle_style)
            sys.stdout.write('\n')
            for angle_style in angle_styles:
                if not (angle_style in angle_styles_selected):
                    continue
                sys.stdout.write('    # '+angle_style2docs[angle_style]+'\n')
            sys.stdout.write('\n')

        if len(dihedral_styles) > 0:
            sys.stdout.write('    dihedral_style hybrid')
            for dihedral_style in dihedral_styles:
                if not (dihedral_style in dihedral_styles_selected):
                    continue
                sys.stdout.write(' ' + dihedral_style)
            sys.stdout.write('\n')
            for dihedral_style in dihedral_styles:
                if not (dihedral_style in dihedral_styles_selected):
                    continue
                sys.stdout.write('    # '+dihedral_style2docs[dihedral_style]+'\n')
            sys.stdout.write('\n')

        if len(improper_styles) > 0:
            sys.stdout.write('    improper_style hybrid')
            for improper_style in improper_styles:
                if not (improper_style in improper_styles_selected):
                    continue
                sys.stdout.write(' ' + improper_style)
            sys.stdout.write('\n')
            for improper_style in improper_styles:
                if not (improper_style in improper_styles_selected):
                    continue
                sys.stdout.write('    # '+improper_style2docs[improper_style]+'\n')
            sys.stdout.write('\n')

        if len(pair_styles) > 0:
            sys.stdout.write('    pair_style hybrid')
            for pair_style in pair_styles:
                if not (pair_style in pair_styles_selected):
                    continue
                sys.stdout.write(' ' + pair_style +
                                 ' ' + pair_style_args[pair_style])
            sys.stdout.write('\n')
            for pair_style in pair_styles:
                sys.stdout.write('    # '+pair_style2docs[pair_style]+'\n')
            sys.stdout.write('\n')

        sys.stdout.write('    pair_modify mix ' + pair_mixing_style + '\n')
        sys.stdout.write('    ' + special_bonds_command + '\n')
        sys.stdout.write('    ' + kspace_style + '\n')
        sys.stdout.write('  } #end of init parameters\n\n')
        sys.stdout.write('}  # ' + ffname + '\n\n')
        
        
        sys.stdout.write("#\n"
                         "# WARNING: The following 1-2, 1-3, and 1-4 weighting parameters were ASSUMED:\n")
        sys.stdout.write("#          " + special_bonds_command + "\n")
        sys.stdout.write("#          (See http://lammps.sandia.gov/doc/special_bonds.html for details)\n")

        #sys.stderr.write(' done.\n')


        if len(lines_templates) > 0:
            sys.stdout.write('\n\n\n\n'
                             '# ---- templates from the original .frc file used for atom type selection: ---\n')
            for line in lines_templates:
                sys.stdout.write('# '+line)

        if len(lines_references) > 0:
            sys.stdout.write('\n\n\n\n'
                             '# ---- references from the original .frc file: ----\n\n')
            for ref_number,lines in sorted(lines_references.items()):
                sys.stdout.write('# reference '+str(ref_number)+'\n')
                for line in lines:
                    sys.stdout.write('# '+line)
                sys.stdout.write('\n')


        if len(lines_warnings) > 0:
            sys.stdout.write('\n\n\n\n'
                             '# ---- additional warnings: ----\n')
            for line in lines_warnings:
                sys.stdout.write(line)


        if filename_in != '':
            file_in.close()




    except InputError as err:
        sys.stderr.write('\n\n' + str(err) + '\n')
        sys.exit(1)



if __name__ == '__main__':
    main()
