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
    However, some improper_styles (such as improper_style class2)
    have no symmetry (at least not for arbitrary choices of parameters).
    These force-field styles, the different permulations of atom-order
    are not equivalent.  So we do not want to rearrange the order of
    the atoms (and bonds) in the match, because the resulting interaction
    is not equivalent.  In this case, this function returns
    the original "match" argument unmodified.

    """

    return match
