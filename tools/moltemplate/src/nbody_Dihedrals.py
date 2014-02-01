from nbody_graph_search import Ugraph

#    To find 4-body "dihedral" interactions, we would use this subgraph:
#
#                              1st bond connects atoms 0 and 1
#       *---*---*---*      =>  2nd bond connects atoms 1 and 2
#       0   1   2   3          3rd bond connects atoms 2 and 3
#                               

bond_pattern = Ugraph([(0,1), (1,2), (2,3)])
# (Ugraph atom indices begin at 0, not 1)


def canonical_order(match):
    """
    When searching for atoms with matching bond patterns GraphMatcher 
    often returns redundant results. We must define a "canonical_order"
    function which sorts the atoms and bonds in a way which is consistent 
    with the type of N-body interaction being considered.
    The atoms (and bonds) in a candidate match are rearranged by the 
    canonical_order().  Then the re-ordered list of atom and bond ids is 
    tested against the list of atom/bond ids in the matches-found-so-far,
    before it is added to the list of interactions found so far.
    (For example, it does not make sense to define a separate 4-body dihedral-
    angle interaction between atoms 1, 2, 3, 4  AND 4, 3, 2, 1.  The dihedral-
    angle is not altered when the order of atoms is reversed, so we arbitrarily
    choose whichever order causes the first atom to have the lower atomID.)

    """

    # match[0][0:3] contains the ID numbers of the 4 atoms in the match
    atom0 = match[0][0]  
    atom1 = match[0][1]
    atom2 = match[0][2]
    atom3 = match[0][3]
    # match[1][0:2] contains the ID numbers of the the 3 bonds
    bond0 = match[1][0]  
    bond1 = match[1][1]
    bond2 = match[1][2]
    if atom0 < atom3:
        #return ((atom0, atom1, atom2, atom3), (bond0, bond1, bond2))  same as:
        return match
    else:
        return ((atom3, atom2, atom1, atom0), (bond2, bond1, bond0))
