from nbody_graph_search import Ugraph

#    To find 4-body "improper" interactions, we would use this subgraph:
#           3
#           *                  1st bond connects atoms 0 and 1
#           |              =>  2nd bond connects atoms 0 and 2
#         _.*._                3rd bond connects atoms 0 and 3
#       *'  0  `*              
#      1         2
#

bond_pattern = Ugraph([(0,1), (0,2), (0,3)])
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
    (For example, it does not make sense to define a separate 4-body improper-
    angle interaction between atoms 1, 2, 3, 4  AND 1, 3, 2, 4.  The improper-
    angle is defined as the angle between planes formed by atoms 1,2,3 & 2,3,4.
    This is not effected by swapping the middle pair of atoms so we arbitrarily
    sort them so that the second atom has a lower atomID than the third atom.)

    """
    atom0 = match[0][0]
    atom1 = match[0][1]
    atom2 = match[0][2]
    atom3 = match[0][3]
    if atom1 <= atom2:
        #return ((atom0,atom1,atom2,atom3), match[1])  same thing as:
        return match
    else:
        return ((atom0,atom2,atom1,atom3), match[1])
