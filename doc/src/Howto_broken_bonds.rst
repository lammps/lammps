Broken Bonds
============

Typically, molecular bond interactions persist for the duration of a
simulation in LAMMPS.  However, some commands break bonds dynamically,
including the following:

* :doc:`bond_style quartic <bond_quartic>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`fix bond/react <fix_bond_react>`
* :doc:`BPM package <Howto_bpm>` bond styles

A bond can break if it is stretched beyond a user-defined threshold or
more generally if other criteria are met.

For the quartic bond style, when a bond is broken its bond type is set
to 0 to effectively break it and pairwise forces between the two atoms
in the broken bond are "turned on".  Angles, dihedrals, etc cannot be
defined for a system when :doc:`bond_style quartic <bond_quartic>` is
used.

Similarly, bond styles in the BPM package are also incompatible with
angles, dihedrals, etc. and when a bond breaks its type is set to zero.
However, in the BPM package one can either turn off all pair interactions
between bonded particles or leave them on, overlaying pair forces on
top of bond forces. To remove pair forces, the special bond list is
dynamically updated.  More details can be found on the :doc:`Howto BPM
<Howto_bpm>` page.

The :doc:`fix bond/break <fix_bond_break>` and :doc:`fix bond/react
<fix_bond_react>` commands allow breaking of bonds within a molecular
topology with may also define angles, dihedrals, etc.  These commands
update internal topology data structures to remove broken bonds, as
well as the appropriate angle, dihedral, etc interactions which
include the bond.  They also trigger a rebuild of the neighbor list
when this occurs, to turn on the appropriate pairwise forces.

Note that when bonds are dumped to a file via the :doc:`dump local
<dump>` command, bonds with type 0 are not included.

The :doc:`delete_bonds <delete_bonds>` command can be used to query
the status of broken bonds with type = 0 or permanently delete them,
e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove

The compute :doc:`count/type <compute_count_type>` command tallies the
current number of bonds (or angles, etc) for each bond (angle, etc)
type.  It also tallies broken bonds with type = 0.

The compute :doc:`nbond/atom <compute_nbond_atom>` command tallies the
current number of bonds each atom is part of, excluding broken bonds
with type = 0.
