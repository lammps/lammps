from nbody_graph_search import Ugraph

#    To find 3-body "angle" interactions, we would use this subgraph:
#
#                               
#       *---*---*           =>  1st bond connects atoms 0 and 1
#       0   1   2               2nd bond connects atoms 1 and 2
#

bond_pattern = Ugraph([(0,1), (1,2)])
# (Ugraph atom indices begin at 0, not 1)


#    The next function eliminates the redundancy between 0-1-2 and 2-1-0:
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
    (For example, it does not make sense to define a separate 3-body angle 
    interaction between atoms 1, 2, 3  AND  3, 2, 1.  This is the same triplet
    of atoms, and the angle between them is the same.)

    """
    # match[0][0:2] contains the ID numbers for the 3 atoms in the match
    atom0 = match[0][0]  
    atom1 = match[0][1]
    atom2 = match[0][2]
    # match[1][0:1] contains the ID numbers for the 2 bonds
    bond0 = match[1][0]
    bond1 = match[1][1]
    if atom0 < atom2:
        #return ((atom0, atom1, atom2), (bond0, bond1))  same thing as:
        return match
    else:
        return ((atom2, atom1, atom0), (bond1, bond0))
