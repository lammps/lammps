from nbody_graph_search import Ugraph

#    To find 4-body "improper" interactions, we would use this subgraph:
#           3
#           *                  1st bond connects atoms 1 and 0
#           |              =>  2nd bond connects atoms 1 and 2
#         _.*._                3rd bond connects atoms 1 and 3
#       *'  1  `*              
#      0         2
#

bond_pattern = Ugraph([(1,0), (1,2), (1,3)])
# (Note: Ugraph atom-index counters begin at 0, not 1)


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

    In this version, we assume that the central atom is the second one,
    and that the interaction is invariant with respect to swapping the first
    and third atoms (I and K).
    (This would be the case if the improper angle was defined as the angle
     between planes JIK and JKL  {assuming the potential function is even}.)

    """
    atom0 = match[0][0]
    atom1 = match[0][1]
    atom2 = match[0][2]
    atom3 = match[0][3]
    # match[1][0:2] contains the ID numbers for the 3 bonds
    bond0 = match[1][0]
    bond1 = match[1][1]
    bond2 = match[1][2]
    if atom0 <= atom2:
        #return ((atom0,atom1,atom2,atom3), (bond0, bond1, bond2))
        # But this is the same thing as:
        return match
    else:
        return ((atom2,atom1,atom0,atom3), (bond1, bond0, bond2))
