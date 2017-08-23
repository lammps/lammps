try:
    from ..nbody_graph_search import Ugraph
except:
    # not installed as a module
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
    before it is added to the list of interactions found so far.  In this
    case we assume the second atom is the central atom (the "hub"), and the
    energy is invariant with respect to permutations of the other 3 atoms.
    So we arbitrarily sort these other 3 atoms in increasing order
    (as well as the bonds which connect the central atom to them).

    """
    atom0 = match[0][0]
    atom1 = match[0][1]
    atom2 = match[0][2]
    atom3 = match[0][3]
    bonds = match[1]
    ab=[(atom1,0), (atom2,1), (atom3,2)]
    ab.sort()
    return ((atom0, ab[0][0], ab[1][0], ab[2][0]),
            (bonds[ab[0][1]], bonds[ab[1][1]], bonds[ab[2][1]]))
