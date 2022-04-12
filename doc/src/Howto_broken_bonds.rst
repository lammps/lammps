Broken Bonds
============

Typically, bond interactions persist for the duration of a simulation
in LAMMPS. However, there are some exceptions that allow for bonds to
break including the :doc:`quartic bond style <bond_quartic>` and the
bond styles in the :doc:`BPM package <Howto_bpm>` which contains the
:doc:`bpm/spring <bond_bpm_spring>` and
:doc:`bpm/rotational <bond_bpm_rotational>` bond styles. In these cases,
a bond can be broken if it is stretched beyond a user-defined threshold.
LAMMPS accomplishes this by setting the bond type to zero such that the
bond force is no longer computed.

Users are normally able to weight the contribution of pair forces to atoms
that are bonded using the :doc:`special_bonds command <special_bonds>`.
When bonds break, this is not always the case. For the quartic bond style,
pair forces are always turned off between bonded particles. LAMMPS does
this via a computational sleight-of-hand. It subtracts the pairwise
interaction as part of the bond computation. When the bond breaks, the
subtraction stops.  For this to work, the pairwise interaction must always
be computed by the :doc:`pair_style <pair_style>` command, whether the bond
is broken or not.  This means that :doc:`special_bonds <special_bonds>` must
be set to 1,1,1. After the bond breaks, the pairwise interaction between the
two atoms is turned on, since they are no longer bonded.

In the BPM package, one can either turn off all pair interactions between
bonded particles or leave them on, overlaying pair forces on top of bond
forces. To remove pair forces, the special bond list is dynamically
updated. More details can be found on the :doc:`Howto BPM <Howto_bpm>`
page.

Bonds can also be broken by fixes which change bond topology, including
:doc:`fix bond/break <fix_bond_break>` and
:doc:`fix bond/react <fix_bond_react>`. These fixes will automatically
trigger a rebuild of the neighbor list and update special bond data structures
when bonds are broken.

Note that when bonds are dumped to a file via the :doc:`dump local <dump>` command, bonds with type 0 are not included.  The
:doc:`delete_bonds <delete_bonds>` command can also be used to query the
status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove

The compute :doc:`nbond/atom <compute_nbond_atom>` can also be used
to tally the current number of bonds per atom, excluding broken bonds.
