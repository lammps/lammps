.. index:: reset_mol_ids

reset_mol_ids command
=====================

Syntax
""""""

Syntax
""""""

.. parsed-literal::

   reset_mol_ids group-ID

* group-ID = ID of group of atoms whose molecule IDs will be reset

Examples
""""""""

.. code-block:: LAMMPS

   reset_mol_ids all
   reset_mol_ids solvent

Description
"""""""""""

Reset molecule IDs for a group of atoms.  This will create a new set
of molecule IDs that are numbered contiguously from 1 to N, if there
are N different molecule IDs used by the group.  Only molecule IDs for
atoms in the specified group are reset.

This can be useful to invoke after performing a reactive molecular
dynamics run with :doc:`fix bond/react <fix_bond_react>`, :doc:`fix
bond/create <fix_bond_create>`, or :doc:`fix bond/break
<fix_bond_break>`. It can also be useful after molecules have been
deleted with :doc:`delete_atoms <delete_atoms>` or after a simulation
which has lost molecules, e.g. via the :doc:`fix evaporate
<fix_evaporate>` command.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`reset_ids <reset_ids>`, :doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/create <fix_bond_create>`,
:doc:`fix bond/break <fix_bond_break>`,
:doc:`fix evaporate <fix_evaporate>`,
:doc:`delete_atoms <delete_atoms>`

**Default:** none
