.. index:: reset_ids

reset_ids command
=================

Syntax
""""""

.. code-block:: LAMMPS

   reset_ids

Examples
""""""""

.. code-block:: LAMMPS

   reset_ids

Description
"""""""""""

Reset atom IDs for the system, including all the global IDs stored
for bond, angle, dihedral, improper topology data.  This will
create a set of IDs that are numbered contiguously from 1 to N
for a N atoms system.

This can be useful to do after performing a "delete\_atoms" command for
a molecular system.  The delete\_atoms compress yes option will not
perform this operation due to the existence of bond topology.  It can
also be useful to do after any simulation which has lost atoms,
e.g. due to atoms moving outside a simulation box with fixed
boundaries (see the "boundary command"), or due to evaporation (see
the "fix evaporate" command).

Note that the resetting of IDs is not really a compression, where gaps
in atom IDs are removed by decrementing atom IDs that are larger.
Instead the IDs for all atoms are erased, and new IDs are assigned so
that the atoms owned by an individual processor have consecutive IDs,
as the :doc:`create_atoms <create_atoms>` command explains.

.. note::

   If this command is used before a :doc:`pair style <pair_style>` is
   defined, an error about bond topology atom IDs not being found may
   result.  This is because the cutoff distance for ghost atom
   communication was not sufficient to find atoms in bonds, angles, etc
   that are owned by other processors.  The :doc:`comm_modify cutoff <comm_modify>` command can be used to correct this issue.
   Or you can define a pair style before using this command.  If you do
   the former, you should unset the comm\_modify cutoff after using
   reset\_ids so that subsequent communication is not inefficient.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`delete_atoms <delete_atoms>`

**Default:** none
