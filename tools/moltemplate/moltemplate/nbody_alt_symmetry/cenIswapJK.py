try:
    from ..nbody_graph_search import Ugraph
except:
    # not installed as a module
    from nbody_graph_search import Ugraph

#    To find 4-body "improper" interactions,
#    (by default, most of the time), we would use this subgraph:
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
    angle interaction between atoms 0, 1, 2, 3  AND 0, 2, 1, 3.
    The "improper angle" is often defined as the angle between planes formed
    by atoms 0,1,2 & 1,2,3.  Alternately, it may instead be defined as the
    angle between the 0,1,2 plane and atom 3.  Either way, this angle does
    not change when swapping the middle pair of atoms (1 and 2)
    (except for a change of sign, which does not matter since the energy functions
     used are typically sign invariant.  Furthermore, neither of OUTER pair of atoms
     are the central atom. There are 3!=6 ways of ordering the remaining 3 atoms.)
    Consequently it does not make sense to define a separate 4-body improper-
    interaction between atoms 0,1,2,3   AS WELL AS between  0,2,1,3.
    So we sort the atoms and bonds so that the first atom has a always has
    a lower atomID than the last atom.  (Later we will check to see if we
    have already defined an interaction between these 4 atoms.  If not then
    we create a new one.)

    """
    atom0 = match[0][0]
    atom1 = match[0][1]
    atom2 = match[0][2]
    atom3 = match[0][3]
    # match[1][0:2] contains the ID numbers for the 3 bonds
    bond0 = match[1][0]
    bond1 = match[1][1]
    bond2 = match[1][2]
    if atom1 <= atom2:
        #return ((atom0,atom1,atom2,atom3), (bond0, bond1, bond2))
        # But this is the same thing as:
        return match
    else:
        return ((atom0,atom2,atom1,atom3), (bond1, bond0, bond2))
