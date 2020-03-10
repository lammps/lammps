.. index:: delete_atoms

delete_atoms command
====================

Syntax
""""""


.. code-block:: LAMMPS

   delete_atoms style args keyword value ...

* style = *group* or *region* or *overlap* or *porosity*

  .. parsed-literal::

       *group* args = group-ID
       *region* args = region-ID
       *overlap* args = cutoff group1-ID group2-ID
         cutoff = delete one atom from pairs of atoms within the cutoff (distance units)
         group1-ID = one atom in pair must be in this group
         group2-ID = other atom in pair must be in this group
       *porosity* args = region-ID fraction seed
         region-ID = region within which to perform deletions
         fraction = delete this fraction of atoms
         seed = random number seed (positive integer)

* zero or more keyword/value pairs may be appended
* keyword = *compress* or *bond* or *mol*

  .. parsed-literal::

       *compress* value = *no* or *yes*
       *bond* value = *no* or *yes*
       *mol* value = *no* or *yes*



Examples
""""""""


.. code-block:: LAMMPS

   delete_atoms group edge
   delete_atoms region sphere compress no
   delete_atoms overlap 0.3 all all
   delete_atoms overlap 0.5 solvent colloid
   delete_atoms porosity cube 0.1 482793 bond yes

Description
"""""""""""

Delete the specified atoms.  This command can be used to carve out
voids from a block of material or to delete created atoms that are too
close to each other (e.g. at a grain boundary).

For style *group*\ , all atoms belonging to the group are deleted.

For style *region*\ , all atoms in the region volume are deleted.
Additional atoms can be deleted if they are in a molecule for which
one or more atoms were deleted within the region; see the *mol*
keyword discussion below.

For style *overlap* pairs of atoms whose distance of separation is
within the specified cutoff distance are searched for, and one of the
2 atoms is deleted.  Only pairs where one of the two atoms is in the
first group specified and the other atom is in the second group are
considered.  The atom that is in the first group is the one that is
deleted.

Note that it is OK for the two group IDs to be the same (e.g. group
*all*\ ), or for some atoms to be members of both groups.  In these
cases, either atom in the pair may be deleted.  Also note that if
there are atoms which are members of both groups, the only guarantee
is that at the end of the deletion operation, enough deletions will
have occurred that no atom pairs within the cutoff will remain
(subject to the group restriction).  There is no guarantee that the
minimum number of atoms will be deleted, or that the same atoms will
be deleted when running on different numbers of processors.

For style *porosity* a specified *fraction* of atoms are deleted
within the specified region.  For example, if fraction is 0.1, then
10% of the atoms will be deleted.  The atoms to delete are chosen
randomly.  There is no guarantee that the exact fraction of atoms will
be deleted, or that the same atoms will be deleted when running on
different numbers of processors.

If the *compress* keyword is set to *yes*\ , then after atoms are
deleted, then atom IDs are re-assigned so that they run from 1 to the
number of atoms in the system.  Note that this is not done for
molecular systems (see the :doc:`atom_style <atom_style>` command),
regardless of the *compress* setting, since it would foul up the bond
connectivity that has already been assigned.  However, the
:doc:`reset_ids <reset_ids>` command can be used after this command to
accomplish the same thing.

Note that the re-assignment of IDs is not really a compression, where
gaps in atom IDs are removed by decrementing atom IDs that are larger.
Instead the IDs for all atoms are erased, and new IDs are assigned so
that the atoms owned by individual processors have consecutive IDs, as
the :doc:`create_atoms <create_atoms>` command explains.

A molecular system with fixed bonds, angles, dihedrals, or improper
interactions, is one where the topology of the interactions is
typically defined in the data file read by the
:doc:`read_data <read_data>` command, and where the interactions
themselves are defined with the :doc:`bond_style <bond_style>`,
:doc:`angle_style <angle_style>`, etc commands.  If you delete atoms
from such a system, you must be careful not to end up with bonded
interactions that are stored by remaining atoms but which include
deleted atoms.  This will cause LAMMPS to generate a "missing atoms"
error when the bonded interaction is computed.  The *bond* and *mol*
keywords offer two ways to do that.

It the *bond* keyword is set to *yes* then any bond or angle or
dihedral or improper interaction that includes a deleted atom is also
removed from the lists of such interactions stored by non-deleted
atoms.  Note that simply deleting interactions due to dangling bonds
(e.g. at a surface) may result in a inaccurate or invalid model for
the remaining atoms.

It the *mol* keyword is set to *yes*\ , then for every atom that is
deleted, all other atoms in the same molecule (with the same molecule
ID) will also be deleted.  This is not done for atoms with molecule ID
= 0, since such an ID is assumed to flag isolated atoms that are not
part of molecules.

.. note::

   The molecule deletion operation in invoked after all individual
   atoms have been deleted using the rules described above for each
   style.  This means additional atoms may be deleted that are not in the
   group or region, that are not required by the overlap cutoff
   criterion, or that will create a higher fraction of porosity than was
   requested.

Restrictions
""""""""""""


The *overlap* styles requires inter-processor communication to acquire
ghost atoms and build a neighbor list.  This means that your system
must be ready to perform a simulation before using this command (force
fields setup, atom masses set, etc).  Since a neighbor list is used to
find overlapping atom pairs, it also means that you must define a
:doc:`pair style <pair_style>` with the minimum force cutoff distance
between any pair of atoms types (plus the :doc:`neighbor <neighbor>`
skin) >= the specified overlap cutoff.

If the :doc:`special_bonds <special_bonds>` command is used with a
setting of 0, then a pair of bonded atoms (1-2, 1-3, or 1-4) will not
appear in the neighbor list, and thus will not be considered for
deletion by the *overlap* styles.  You probably don't want to be
deleting one atom in a bonded pair anyway.

The *bond yes* option cannot be used with molecular systems defined
using molecule template files via the :doc:`molecule <molecule>` and
:doc:`atom_style template <atom_style>` commands.

Related commands
""""""""""""""""

:doc:`create_atoms <create_atoms>`, :doc:`reset_ids <reset_ids>`

Default
"""""""

The option defaults are compress = yes, bond = no, mol = no.
