.. index:: reset_atom_ids

reset_atom_ids command
======================

Syntax
""""""

.. code-block:: LAMMPS

   reset_atom_ids keyword values ...

   * zero or more keyword/value pairs may be appended
   * keyword = *sort*

.. parsed-literal::

   *sort* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   reset_atom_ids
   reset_atom_ids sort yes

Description
"""""""""""

Reset atom IDs for the system, including all the global IDs stored
for bond, angle, dihedral, improper topology data.  This will
create a set of IDs that are numbered contiguously from 1 to N
for a N atoms system.

This can be useful to do after performing a "delete_atoms" command for
a molecular system.  The delete_atoms compress yes option will not
perform this operation due to the existence of bond topology.  It can
also be useful to do after any simulation which has lost atoms,
e.g. due to atoms moving outside a simulation box with fixed
boundaries (see the "boundary command"), or due to evaporation (see
the "fix evaporate" command).

If the *sort* keyword is used with a setting of *yes*, then the
assignment of new atom IDs will be the same no matter how many
processors LAMMPS is running on.  This is done by first doing a
spatial sort of all the atoms into bins and sorting them within each
bin.  Because the set of bins is independent of the number of
processors, this enables a consistent assignment of new IDs to each
atom.

This can be useful to do after using the "create_atoms" command and/or
"replicate" command.  In general those commands do not guarantee
assignment of the same atom ID to the same physical atom when LAMMPS
is run on different numbers of processors.  Enforcing consistent IDs
can be useful for debugging or comparing output from two different
runs.

Note that the spatial sort requires communication of atom IDs and
coordinates between processors in an all-to-all manner.  This is done
efficiently in LAMMPS, but it is more expensive than how atom IDs are
reset without sorting.

Note that whether sorting or not, the resetting of IDs is not a
compression, where gaps in atom IDs are removed by decrementing atom
IDs that are larger.  Instead the IDs for all atoms are erased, and
new IDs are assigned so that the atoms owned by an individual
processor have consecutive IDs, as the :doc:`create_atoms
<create_atoms>` command explains.

.. note::

   If this command is used before a :doc:`pair style <pair_style>` is
   defined, an error about bond topology atom IDs not being found may
   result.  This is because the cutoff distance for ghost atom
   communication was not sufficient to find atoms in bonds, angles, etc
   that are owned by other processors.  The :doc:`comm_modify cutoff <comm_modify>` command can be used to correct this issue.
   Or you can define a pair style before using this command.  If you do
   the former, you should unset the comm_modify cutoff after using
   reset_atom_ids so that subsequent communication is not inefficient.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`delete_atoms <delete_atoms>`

Default
"""""""

By default, *sort* is no.
