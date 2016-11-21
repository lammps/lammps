#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
#         http://www.chem.ucsb.edu/~sheagroup
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2012, Regents of the University of California
# All rights reserved.


import sys
from nbody_graph_search import *
#from collections import namedtuple
if sys.version < '2.7':
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
from collections import defaultdict
from ttree_lex import MatchesPattern, MatchesAll, InputError
#import gc


def GenInteractions_int(G_system, 
                        g_bond_pattern,
                        typepattern_to_coefftypes,
                        canonical_order, #function to sort atoms and bonds
                        atomtypes_int2str,
                        bondtypes_int2str,
                        report_progress = False): #print messages to sys.stderr?
    """
    GenInteractions() automatically determines a list of interactions 
    present in a system of bonded atoms (argument "G_system"),
    which satisfy the bond topology present in "g_bond_pattern", and 
    satisfy the atom and bond type requirements in "typepattern_to_coefftypes".

    Whenever a set of atoms in "G_system" are bonded together in a way which
    matches "g_bond_pattern", and when the atom and bond types is consistent 
    with one of the entries in "typepattern_to_coefftypes", the corresponding 
    list of atoms from G_system is appended to the list of results.

    These results (the list of lists of atoms participating in an interaction)
    are organized according their corresponding "coefftype", a string 
    which identifies the type of interaction they obey as explained above.
    results are returned as a dictionary using "coefftype" as the lookup key.

    Arguments:
    
     -- typepattern_to_coefftypes is a list of 2-tuples --
    The first element of the 2-tuple is the "typepattern".
    It contains a string describing a list of atom types and bond types.
    The typepattern is associated with a "coefftype",
    which is the second element of the 2-tuple.  This is a string 
    which identifies the type of interaction between the atoms.  
    Later on, this string can be used to lookup the force field 
    parameters for this interaction elsewhere.)

     -- Arguments: G_system, g_bond_pattern, atomtypes_int2str, bondtypes_int2str --

    G_system stores a list of atoms and bonds, and their attributes in 
    "Ugraph" format.  In this format:
    Atom ID numbers are represented by indices into the G_system.verts[] list.
    Bond ID numbers are represented by indices into the G_system.edges[] list.
    Atom types are represented as integers in the G_system.verts[i].attr list.
    Bond types are represented as integers in the G_system.edges[i].attr list.
    They are converted into strings using
    atomtypes_int2str, and bondtypes_int2str.

    g_bond_pattern is a graph which specifies the type of bonding between
    the atoms required for a match. It is in Ugraph format (however the 
    atom and bond types are left blank.)

    Atom and bond types are supplied by the user in string format. (These 
    strings typically encode integers, but could be any string in principle.)
    The string-version of the ith atom type is stored in 
       atomtypes_int2str[ G_system.verts[i].attr ]
    The string-version of the ith bond type is stored in 
       bondtypes_int2str[ G_system.edges[i].attr ]

     -- The "canonical_order" argument: --

    The search for atoms with a given bond pattern often yields 
    redundant matches.  There is no difference for example 
    between the angle formed between three consecutively 
    bonded atoms (named, 1, 2, 3, for example), and the
    angle between the same atoms in reverse order (3, 2, 1). 
    However both triplets of atoms will be returned by the subgraph-
    matching algorithm when searching for ALL 3-body interactions.)

    To eliminate this redundancy, the caller must supply a "canonical_order" 
    argument.  This is a function which sorts the atoms and bonds in a way
    which is consistent with the type of N-body interaction being considered.
    The atoms (and bonds) in a candidate match are rearranged by the 
    canonical_order().  Then the re-ordered list of atom and bond ids is 
    tested against the list of atom/bond ids in the matches-found-so-far,
    before it is added.

    """

    if report_progress:
        startatomid = 0
        sys.stderr.write('  searching for matching bond patterns:\n')
        sys.stderr.write('    0%')

    # Figure out which atoms from "G_system" bond together in a way which 
    # matches the "g_bond_pattern" argument.  Organize these matches by 
    # atom and bond types and store all of the non-redundant ones in
    # the "interactions_by_type" variable.

    gm = GraphMatcher(G_system, g_bond_pattern)

    interactions_by_type = defaultdict(list)

    for atombondids in gm.Matches():
        # "atombondids" is a tuple.
        #  atombondids[0] has atomIDs from G_system corresponding to g_bond_pattern
        #     (These atomID numbers are indices into the G_system.verts[] list.)
        #  atombondids[1] has bondIDs from G_system corresponding to g_bond_pattern
        #     (These bondID numbers are indices into the G_system.edges[] list.)

        # It's convenient to organize the list of interactions-between-
        # atoms in a dictionary indexed by atomtypes and bondtypes.
        # (Because many atoms and bonds typically share the same type, 
        #  organizing the results this way makes it faster to check 
        #  whether a given interaction matches a "typepattern" defined
        #  by the user.  We only have to check once for the whole group.)

        atombondtypes = \
            (tuple([G_system.GetVert(Iv).attr for Iv in atombondids[0]]),
             tuple([G_system.GetEdge(Ie).attr for Ie in atombondids[1]]))

        interactions_by_type[atombondtypes].append(atombondids)

        if report_progress:
            # GraphMatcher.Matches() searches for matches in an order
            # that selects a different atomid number from G_system, 
            # starting at 0, and continuing up to the number of atoms (-1)
            # in the system (G_system.nv-1), and using this as the first
            # atom in the match (ie match[0][0]). This number can be used 
            # to guess much progress has been made so far.
            oldatomid = startatomid
            startatomid = atombondids[0][0]
            percent_complete = (100 * startatomid) // G_system.GetNumVerts()
            # report less often as more progress made
            if percent_complete <= 4:
                old_pc = (100 * oldatomid) // G_system.GetNumVerts()
                if percent_complete > old_pc:
                    sys.stderr.write('  '+str(percent_complete)+'%')
            elif percent_complete <= 8:
                pc_d2    = (100 * startatomid) // (2*G_system.GetNumVerts())
                oldpc_d2 = (100 * oldatomid)   // (2*G_system.GetNumVerts())
                if pc_d2 > oldpc_d2:
                    sys.stderr.write('  '+str(percent_complete)+'%')
            elif percent_complete <= 20:
                pc_d4    = (100 * startatomid) // (4*G_system.GetNumVerts())
                oldpc_d4 = (100 * oldatomid)   // (4*G_system.GetNumVerts())
                if pc_d4 > oldpc_d4:
                    sys.stderr.write('  '+str(percent_complete)+'%')
            else:
                pc_d10    = (100 * startatomid) // (10*G_system.GetNumVerts())
                oldpc_d10 = (100 * oldatomid)   // (10*G_system.GetNumVerts())
                if pc_d10 > oldpc_d10:
                    sys.stderr.write('  '+str(percent_complete)+'%')

    if report_progress:
        sys.stderr.write('  100%\n')
        #sys.stderr.write('    ...done\n')
        #sys.stderr.write('    Looking up available atom and bond types...')

    #coefftype_to_atomids = defaultdict(list)
    #abids_to_coefftypes = defaultdict(list)
    coefftype_to_atomids = OrderedDict()
    abids_to_coefftypes = OrderedDict()




    # -------------------- reporting progress -----------------------
    if report_progress:
        # The next interval of code is not technically necessary, but it makes 
        # the printed output easier to read by excluding irrelevant interactions
        # Now, test each match to see if the atoms and bonds involved match
        # any of the type-patterns in the "typepattern_to_coefftypes" argument.

        types_atoms_all_str = set([])
        types_bonds_all_str = set([])
        for typepattern, coefftype in typepattern_to_coefftypes:
            for atombondtypes, abidslist in interactions_by_type.items():
                for Iv in atombondtypes[0]:
                    types_atoms_all_str.add(atomtypes_int2str[Iv])
                for Ie in atombondtypes[1]:
                    types_bonds_all_str.add(bondtypes_int2str[Ie])
    # ------------------ reporting progress (end) -------------------



    count = 0
    for typepattern, coefftype in typepattern_to_coefftypes:


        # ------------------ reporting progress -----------------------
        # The next interval of code is not technically necessary, but it makes 
        # the printed output easier to read by excluding irrelevant interactions

        if report_progress:

            # Check to see if the atoms or bonds referred to in typepattern
            # are (potentially) satisfied by any of the atoms present in the system.
            # If any of the required atoms for this typepattern are not present
            # in this system, then skip to the next typepattern.
            atoms_available_Iv = [False for Iv in range(0, g_bond_pattern.GetNumVerts())]
            for Iv in range(0, g_bond_pattern.GetNumVerts()):
                for type_atom_str in types_atoms_all_str:
                    if MatchesPattern(type_atom_str, typepattern[Iv]):
                        atoms_available_Iv[Iv] = True
            atoms_available = True
            for Iv in range(0, g_bond_pattern.GetNumVerts()):
                if not atoms_available_Iv[Iv]:
                    atoms_available = False

            bonds_available_Ie = [False for Ie in range(0, g_bond_pattern.GetNumEdges())]
            for Ie in range(0, g_bond_pattern.GetNumEdges()):
                for type_bond_str in types_bonds_all_str:
                    if MatchesPattern(type_bond_str,
                                      typepattern[g_bond_pattern.GetNumVerts()+Ie]):
                        bonds_available_Ie[Ie] = True
            bonds_available = True
            for Ie in range(0, g_bond_pattern.GetNumEdges()):
                if not bonds_available_Ie[Ie]:
                    bonds_available = False

            if atoms_available and bonds_available:

                # Explanation:
                # (Again) only if ALL of the atoms and bond requirements for
                # this typepattern are satisfied by at least SOME of the atoms
                # present in the this system, ...THEN print a status message.
                # (Because for complex all-atom force-fields, the number of
                # possible atom types, and typepatterns far exceeds the number
                # of atom types typically present in the system.  Otherwise
                # hundreds of kB of irrelevant information can be printed.)

                sys.stderr.write('    checking '+coefftype+' type requirements:'
                                 #' (atom-types,bond-types) '
                                 '\n     '+str(typepattern)+'\n')

        # ------------------ reporting progress (end) -------------------





        for atombondtypes, abidslist in interactions_by_type.items():
            # express atom & bond types in a tuple of the original string format
            types_atoms  = [atomtypes_int2str[Iv] for Iv in atombondtypes[0]]
            types_bonds  = [bondtypes_int2str[Ie] for Ie in atombondtypes[1]]
            type_strings = types_atoms + types_bonds
            # use string comparisons to check for a match with typepattern
            if MatchesAll(type_strings, typepattern): #<-see "ttree_lex.py"
                for abids in abidslist:

                    # Re-order the atoms (and bonds) in a "canonical" way. 
                    # Only add new interactions to the list after re-ordering 
                    # them and checking that they have not been added earlier.
                    # (...well not when using the same coefftype at least.
                    #  This prevents the same triplet of atoms from 
                    #  being used to calculate the bond-angle twice: 
                    #  once for 1-2-3 and 3-2-1, for example.)
                    abids = canonical_order(abids)
                    redundant = False
                    if abids in abids_to_coefftypes:
                        coefftypes = abids_to_coefftypes[abids]
                        if coefftype in coefftypes:
                            redundant = True

                    if not redundant:
                        # (It's too bad python does not
                        #  have an Ordered defaultdict)
                        if coefftype in coefftype_to_atomids:
                            coefftype_to_atomids[coefftype].append(abids[0])
                        else:
                            coefftype_to_atomids[coefftype]=[abids[0]]
                        if abids in abids_to_coefftypes:
                            abids_to_coefftypes[abids].append(coefftype)
                        else:
                            abids_to_coefftypes[abids] = [coefftype]
                        count += 1

    if report_progress:
        sys.stderr.write('  (found '+
                         str(count)+' non-redundant matches)\n')

    return coefftype_to_atomids








def GenInteractions_str(bond_pairs,
                        g_bond_pattern,
                        typepattern_to_coefftypes,
                        canonical_order, #function to sort atoms and bonds
                        atomids_str,
                        atomtypes_str,
                        bondids_str,
                        bondtypes_str,
                        report_progress = False): #print messages to sys.stderr?


    assert(len(atomids_str) == len(atomtypes_str))
    assert(len(bondids_str) == len(bondtypes_str))
    # The atomids and atomtypes and bondtypes are strings.
    # First we assign a unique integer id to each string.

    atomids_str2int = {}
    atomtypes_str2int = {}
    atomtypes_int2str = []
    atomtype_int = 0
    for i in range(0, len(atomids_str)):
        if atomids_str[i] in atomids_str2int:
            raise InputError('Error: multiple atoms have the same id ('+
                             str(atomids_str[i])+')')
        atomids_str2int[atomids_str[i]] = i
        #atomtypes_int = len(atomtypes_int)+1
        if (not (atomtypes_str[i] in atomtypes_str2int)):
            atomtypes_str2int[atomtypes_str[i]] = atomtype_int
            atomtypes_int2str.append(atomtypes_str[i])
            atomtype_int += 1
        #atomtypes_int.append(atomtype_int)

    bondids_str2int = {}
    bondtypes_str2int = {}
    bondtypes_int2str = []
    bondtype_int = 0
    for i in range(0, len(bondids_str)):
        if bondids_str[i] in bondids_str2int:
            raise InputError('Error: multiple bonds have the same id ('+
                             str(bondids_str[i])+')')
        bondids_str2int[bondids_str[i]] = i
        #bondtype_int = len(bondtypes_int)+1
        if (not (bondtypes_str[i] in bondtypes_str2int)):
            bondtypes_str2int[bondtypes_str[i]] = bondtype_int
            bondtypes_int2str.append(bondtypes_str[i])
            bondtype_int += 1

    # Now convert "bond_pairs" into the UGraph format
    G_system = Ugraph()
    for iv in range(0, len(atomtypes_str)):
        G_system.AddVertex(iv, atomtypes_str2int[atomtypes_str[iv]])

    for ie in range(0, len(bond_pairs)):
        atomid1_str  = bond_pairs[ie][0]
        atomid2_str  = bond_pairs[ie][1]
        if (atomid1_str not in atomids_str2int):
            raise InputError('Error in Bonds Section:\n'
                             '  '+atomid1_str+' is not defined in Atoms section\n')
        if (atomid2_str not in atomids_str2int):
            raise InputError('Error in Bonds Section:\n'
                             '  '+atomid2_str+' is not defined in Atoms section\n')
        G_system.AddEdge(atomids_str2int[atomid1_str],
                         atomids_str2int[atomid2_str],
                         bondtypes_str2int[bondtypes_str[ie]])

    coefftype_to_atomids_int = GenInteractions_int(G_system,
                                                   g_bond_pattern,
                                                   typepattern_to_coefftypes,
                                                   canonical_order,
                                                   atomtypes_int2str,
                                                   bondtypes_int2str,
                                                   report_progress)
    coefftype_to_atomids_str = OrderedDict()
    for coefftype, atomidss_int in coefftype_to_atomids_int.items():
        if report_progress:
            sys.stderr.write('    processing coefftype: '+str(coefftype)+'\n')
        for atomids_int in atomidss_int:
            if coefftype in coefftype_to_atomids_str:
                coefftype_to_atomids_str[coefftype].append(
                [atomids_str[iv] for iv in atomids_int])
            else:
                coefftype_to_atomids_str[coefftype] = \
                [[atomids_str[iv] for iv in atomids_int]]
        #gc.collect()

    return coefftype_to_atomids_str


