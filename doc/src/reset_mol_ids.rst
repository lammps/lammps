.. index:: reset_mol_ids

reset_mol_ids command
=====================

Syntax
""""""

Syntax
""""""

.. parsed-literal::

   reset_mol_ids group-ID keyword value ...

* group-ID = ID of group of atoms whose molecule IDs will be reset
* zero or more keyword/value pairs may be appended
* keyword = *offset*

  .. parsed-literal::

       *offset* value = *Noffset*
       
Examples
""""""""

.. code-block:: LAMMPS

   reset_mol_ids all
   reset_mol_ids solvent offset 1000

Description
"""""""""""

Reset molecule IDs for a group of atoms.  This will create a new set
of molecule IDs that are numbered contiguously from 1 to N, if there
are N different molecules in the group.  Only molecule IDs for atoms
in the specified group are reset.

For purposes of this operation, molecules are identified by the bond
topology of the system, not by the current molecule IDs.  A molecule
is a set of atoms, each is which is bonded to one or more atoms in the
set.  If an atom is not bonded to any other atom, it becomes its own
1-atom molecule.  Once new molecules are identified, this command will
overwrite the current molecule ID for each atom with a new ID.

This can be a useful operation to perform after running reactive
molecular dynamics run with :doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/create <fix_bond_create>`, or :doc:`fix bond/break
<fix_bond_break>`, all of which can change molecule topologies. It can
also be useful after molecules have been deleted with the
:doc:`delete_atoms <delete_atoms>` command or after a simulation which
has lost molecules, e.g. via the :doc:`fix evaporate <fix_evaporate>`
command.

.. note::

   The same as explained for the :doc:`compute fragment/atom
   <compute_cluster_atom>` command, molecules are identified using the
   current bond topology within each fragment.  This will not account
   for bonds broken by the :doc:`bond_style quartic <bond_quartic>`
   command because it does not perform a full update of the bond
   topology data structures within LAMMPS.

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
