.. index:: compute contact/atom

compute contact/atom command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID contact/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* contact/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all contact/atom

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

**Output info:**

This compute calculates a per-atom vector, whose values can be
accessed by any command that uses per-atom values from a compute as
input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

The per-atom vector values will be a number >= 0.0, as explained
above.

Restrictions
""""""""""""

This compute requires that atoms store a radius as defined by the
:doc:`atom_style sphere <atom_style>` command.

Related commands
""""""""""""""""

:doc:`compute coord/atom <compute_coord_atom>`

**Default:** none
