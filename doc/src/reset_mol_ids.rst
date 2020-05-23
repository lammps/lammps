.. index:: reset_mol_ids

reset_mol_ids command
=================

Syntax
""""""

.. code-block:: LAMMPS

   reset_mol_ids

Examples
""""""""

.. code-block:: LAMMPS

   reset_mol_ids

Description
"""""""""""

Reset molecule IDs for the system.  This will create a set of molecule
IDs that are numbered contiguously from 1 to N for a N molecules
system.

This can be useful to do after performing a reactive molecular
dynamics run with :doc:`fix bond/react <fix_bond_react>`.  It can also
be useful to do after any simulation which has lost atoms, e.g. due to
atoms moving outside a simulation box with fixed boundaries (see the
"boundary command"), or due to evaporation (see the "fix evaporate"
command).

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`reset_ids <reset_ids>`, :doc:`fix bond/react <fix_bond_react>`

**Default:** none
