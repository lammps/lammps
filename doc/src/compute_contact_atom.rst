.. index:: compute contact/atom

compute contact/atom command
============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID contact/atom group2-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* contact/atom = style name of this compute command
* group2-ID = optional argument to restrict which atoms to consider for contacts (see below)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all contact/atom
   compute 1 all contact/atom mygroup

Description
"""""""""""

Define a computation that calculates the number of contacts
for each atom in a group.

The contact number is defined for finite-size spherical particles as
the number of neighbor atoms which overlap the central particle,
meaning that their distance of separation is less than or equal to the
sum of the radii of the two particles.

The value of the contact number will be 0.0 for atoms not in the
specified compute group.

The optional *group2-ID* argument allows to specify from which group atoms
contribute to the coordination number. Default setting is group 'all'.

Output info
"""""""""""

This compute calculates a per-atom vector, whose values can be
accessed by any command that uses per-atom values from a compute as
input.  See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

The per-atom vector values will be a number :math:`\ge 0.0`, as explained
above.

Restrictions
""""""""""""

This compute is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

This compute requires that atoms store a radius as defined by the
:doc:`atom_style sphere <atom_style>` command.

Related commands
""""""""""""""""

:doc:`compute coord/atom <compute_coord_atom>`

Default
"""""""

*group2-ID* = all
