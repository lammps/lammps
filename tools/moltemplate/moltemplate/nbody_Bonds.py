try:
    from .nbody_graph_search import Ugraph
except (ImportError, SystemError, ValueError):
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


#    The next function eliminates the redundancy between 0-1 and 1-0:
def canonical_order(match):
    """
    It does not make sense to define a separate bond between atoms 1 and 2,
    and between atoms 2 and 1.  This function will swap the atoms in the bond
    if the first atom > second atom.

    """
    # match[0][0:2] contains the ID numbers for the 2 atoms in the match
    atom0 = match[0][0]
    atom1 = match[0][1]
    # match[1][0:1] contains the ID numbers for the 1 bond
    bond0 = match[1][0]
    if atom0 < atom1:
        # return ((atom0, atom1), (bond0))  same thing as:
        return match
    else:
        return ((atom1, atom0), (bond0))
