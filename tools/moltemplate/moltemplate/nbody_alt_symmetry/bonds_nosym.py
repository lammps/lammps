try:
    from .nbody_graph_search import Ugraph
except (SystemError, ValueError):
    # not installed as a package
    from nbody_graph_search import Ugraph


#    To find 2-body "bond" interactions, we would use this subgraph:
#
#
#       *---*           =>  one bond connects atoms 0 and 1
#       0   1
#

bond_pattern = Ugraph([(0, 1)])
# (Ugraph atom indices begin at 0, not 1)


#    The next function eliminates the redundancy between 0-1-2 and 2-1-0:
def canonical_order(match):
    """
    When searching for atoms with matching bond patterns GraphMatcher
    often returns redundant results. We must define a "canonical_order"
    function which sorts the atoms and bonds in a way which is consistent
    with the type of N-body interaction being considered.
    However, occasionally we DON'T want to modify the atom order.
    In this case, this function returns
    the original "match" argument unmodified.

    """

    return match
