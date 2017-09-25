try:
    from ..nbody_graph_search import Ugraph
except:
    # not installed as a module
    from nbody_graph_search import Ugraph

# This file defines how improper interactions are generated in AMBER (GAFF).
# To use it, add "(gaff_imp.py)" to the name of the "Data Impropers By Type"
# section, and make sure this file is located in the "common" directory.
# For example:
# write_once("Data Impropers By Type (gaff_imp.py)") {
#   ...
# }


#    To find 4-body "improper" interactions,
#    (by default, most of the time), we would use this subgraph:
#           0
#           *                  1st bond connects atoms 2 and 0
#           |              =>  2nd bond connects atoms 2 and 1
#         _.*._                3rd bond connects atoms 2 and 3
#       *'  2  `*
#      1         3
#
# In AMBER/GAFF, the central atom is the third atom ("2").
# http://archive.ambermd.org/201307/0519.html
# This differs from other force-fields.
# We take this detail into account in the line below:

bond_pattern = Ugraph([(2,0), (2,1), (2,3)])

# As with other force-fields, the improper-angle is the angle between the planes
# defined by the first three atoms (0,1,2) and last three atoms (1,2,3).
# (This is implemented in LAMMPS using an improper_style which requires
#  that the atoms in the interaction will be listed in this order: 0,1,2,3.)

def canonical_order(match):
    """
    Before defining a new interaction, we must check to see if an
    interaction between these same 4 atoms has already been created
     (perhaps listed in a different, but equivalent order).
    If we don't check for this this, we will create many unnecessary redundant
    interactions (which can slow down he simulation).
    To avoid this, I define a "canonical_order" function which sorts the atoms
    and bonds in a way which is consistent with the symmetry of the interaction
    being generated...  Later the re-ordered list of atom and bond ids will be
    tested against the list of atom/bond ids in the matches-found-so-far,
    before it is added to the list of interactions found so far.  Note that
    the energy of an improper interactions is a function of the improper angle.
    The "improper angle" is often defined as the angle between planes formed
    by atoms 0,1,2 & 1,2,3.  (Alternately, it is sometimes defined as the
    angle between the 0,1,2 plane and atom 3.)
    This angle does not change when swapping the OUTER pair of atoms (0 and 3)
    (except for a change of sign, which does not matter since the energy functions
     used are typically sign invariant.  Furthermore, neither of OUTER pair of atoms
     are the central atom. There are 3!=6 ways of ordering the remaining 3 atoms.)
    Consequently it does not make sense to define a separate 4-body improper-
    interaction between atoms 0,1,2,3   AS WELL AS between  3,1,2,0.
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
    if atom0 <= atom3:
        #return ((atom0,atom1,atom2,atom3), (bond0, bond1, bond2))
        # But this is the same thing as:
        return match
    else:
        return ((atom3,atom1,atom2,atom0), (bond2,bond1,bond0))
