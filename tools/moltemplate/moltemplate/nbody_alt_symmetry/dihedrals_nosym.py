try:
    from ..nbody_graph_search import Ugraph
except:
    # not installed as a module
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
    However, some dihedral_styles (such as dihedral_style class2)
    have no symmetry (at least not for arbitrary choices of parameters).
    These force-field styles, the different permulations of atom-order
    are not equivalent.  So we do not want to rearrange the order of
    the atoms (and bonds) in the match, because the formula for the
    interaction between atoms 1,2,3,4 is not the same as the formula
    for the interaction between atoms 4,3,2,1.
    In this case, this function returns
    the original "match" argument unmodified.

    """

    return match
