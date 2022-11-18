.. index:: reset_image_flags

reset image_flags sub-command
=============================

Syntax
""""""

.. parsed-literal::

   reset image_flags group-ID

* group-ID = ID of group of atoms whose image flags will be reset

Examples
""""""""

.. code-block:: LAMMPS

   reset image_flags all
   reset image_flags mobile

Description
"""""""""""

.. versionadded:: TBD

Reset the image flags of atoms so that molecular fragments with atoms in
different periodic images retain their geometry.  This avoids
inconsistent image flags resulting from resetting them all to zero.

Only images flags for atoms in the specified group are reset; all others
remain unchanged.  No check is made whether the group covers complete
molecule fragments and thus whether the command will result in
inconsistent image flags.

Molecular fragments are identified by the algorithm used :doc:`compute
fragment/atom <compute_cluster_atom>`.  For each fragment the average of
the largest and the smallest image flag in each direction across all
atoms in the fragment is computed and subtracted from the current image
flag in the same direction.

This can be a useful operation to perform after running longer
equilibration runs of mobile systems where molecules would pass through
the system multiple times and thus produce non-zero image flags.

.. note::

   Same as explained for the :doc:`compute fragment/atom
   <compute_cluster_atom>` command, molecules are identified using the
   current bond topology.  This will **not** account for bonds broken by
   the :doc:`bond_style quartic <bond_quartic>` command, because this
   bond style does not perform a full update of the bond topology data
   structures within LAMMPS.  In that case, using the :doc:`delete_bonds
   all bond 0 remove <delete_bonds>` will permanently delete such
   broken bonds and should thus be used first.

Restrictions
""""""""""""

The command can only be used when the atom style supports bonds.

Related commands
""""""""""""""""

:doc:`reset mol_ids <reset_mol_ids>`,
:doc:`compute fragment/atom <compute_cluster_atom>`

Defaults
""""""""

none
