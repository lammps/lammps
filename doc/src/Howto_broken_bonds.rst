Broken Bonds
============

Typically, bond interactions persist for the duration of a simulation
in LAMMPS.  However, some commands allow bonds to break, including the
following:

* :doc:`bond_style quartic <bond_quartic>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`fix bond/react <fix_bond_react>`
* :doc:`BPM package <Howto_bpm>` bond styles

A bond can break if it is stretched beyond a user-defined threshold or
more generally if other criteria are met.

For the quartic bond style, when a bond is broken its bond type is set
to 0 and pairwise forces between the two atoms in the broken bond are
"turned on".  Angles, dihedrals, etc cannot be defined for the system
when :doc:`bond_style quartic <bond_quartic>` is used.

The :doc:`fix bond/break <fix_bond_break>` and :doc:`fix bond/react
<fix_bond_react>` commands allow breaking of bonds within a molecular
topology with also defines angles, dihedrals, etc. These fixes will
update internal topology data structures when bonds are broken, so
that the appropriate angle, dihederal, etc interactions are also
turned off.  They will also trigger a rebuild of the neighbor list
when this occurs, to turn on the appropriate pairwise forces.

In the BPM package, one can either turn off all pair interactions
between bonded particles or leave them on, overlaying pair forces on
top of bond forces. To remove pair forces, the special bond list is
dynamically updated. More details can be found on the :doc:`Howto BPM
<Howto_bpm>` page.

Note that when bonds are dumped to a file via the :doc:`dump local
<dump>` command, bonds with type 0 are not included.

The :doc:`delete_bonds <delete_bonds>` command can also be used to
query the status of broken bonds (type = 0) or permanently delete
them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove

The compute :doc:`count/type bond <compute_count_type>` command
tallies the current number of bonds for each bond type.  It also
tallies broken bonds with type = 0.

The compute :doc:`nbond/atom <compute_nbond_atom>` command tallies the
current number of bonds each atom is part of, excluding broken bonds
with type = 0.
