.. index:: reset_mol_ids

reset_mol_ids command
=====================

Syntax
""""""

.. parsed-literal::

   reset_mol_ids group-ID keyword value ...

* group-ID = ID of group of atoms whose molecule IDs will be reset
* zero or more keyword/value pairs may be appended
* keyword = *compress* or *offset* or *single*

  .. parsed-literal::

       *compress* value = *yes* or *no*
       *offset* value = *Noffset* >= -1
       *single* value = *yes* or *no* to treat single atoms (no bonds) as molecules

Examples
""""""""

.. code-block:: LAMMPS

   reset_mol_ids all
   reset_mol_ids all offset 10 single yes
   reset_mol_ids solvent compress yes offset 100
   reset_mol_ids solvent compress no

Description
"""""""""""

Reset molecule IDs for a group of atoms based on current bond
connectivity.  This will typically create a new set of molecule IDs
for atoms in the group.  Only molecule IDs for atoms in the specified
group are reset; molecule IDs for atoms not in the group are not
changed.

For purposes of this operation, molecules are identified by the current
bond connectivity in the system, which may or may not be consistent with
the current molecule IDs.  A molecule in this context is a set of atoms
connected to each other with explicit bonds.  The specific algorithm
used is the one of :doc:`compute fragment/atom <compute_cluster_atom>`
Once the molecules are identified and a new molecule ID computed for
each, this command will update the current molecule ID for all atoms in
the group with the new molecule ID.  Note that if the group excludes
atoms within molecules, one (physical) molecule may become two or more
(logical) molecules.  For example if the group excludes atoms in the
middle of a linear chain, then each end of the chain is considered an
independent molecule and will be assigned a different molecule ID.

This can be a useful operation to perform after running reactive
molecular dynamics run with :doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/create <fix_bond_create>`, or :doc:`fix bond/break
<fix_bond_break>`, all of which can change molecule topologies. It can
also be useful after molecules have been deleted with the
:doc:`delete_atoms <delete_atoms>` command or after a simulation which
has lost molecules, e.g. via the :doc:`fix evaporate <fix_evaporate>`
command.

The *compress* keyword determines how new molecule IDs are computed.  If
the setting is *yes* (the default) and there are N molecules in the
group, the new molecule IDs will be a set of N contiguous values.  See
the *offset* keyword for details on selecting the range of these values.
If the setting is *no*, the molecule ID of every atom in the molecule
will be set to the smallest atom ID of any atom in the molecule.

The *single* keyword determines whether single atoms (not bonded to
another atom) are treated as one-atom molecules or not, based on the
*yes* or *no* setting.  If the setting is *no* (the default), their
molecule IDs are set to 0.  This setting can be important if the new
molecule IDs will be used as input to other commands such as
:doc:`compute chunk/atom molecule <compute_chunk_atom>` or :doc:`fix
rigid molecule <fix_rigid>`.

The *offset* keyword is only used if the *compress* setting is *yes*.
Its default value is *Noffset* = -1.  In that case, if the specified
group is *all*, then the new compressed molecule IDs will range from 1
to N.  If the specified group is not *all* and the largest molecule ID
of atoms outside that group is M, then the new compressed molecule IDs will
range from M+1 to M+N, to avoid collision with existing molecule
IDs.  If an *Noffset* >= 0 is specified, then the new compressed
molecule IDs will range from *Noffset*\ +1 to *Noffset*\ +N.  If the group
is not *all* there may be collisions with the molecule IDs of other atoms.

.. note::

   The same as explained for the :doc:`compute fragment/atom
   <compute_cluster_atom>` command, molecules are identified using the
   current bond topology.  This will not account for bonds broken by
   the :doc:`bond_style quartic <bond_quartic>` command because it
   does not perform a full update of the bond topology data structures
   within LAMMPS.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`reset_atom_ids <reset_atom_ids>`, :doc:`fix bond/react <fix_bond_react>`,
:doc:`fix bond/create <fix_bond_create>`,
:doc:`fix bond/break <fix_bond_break>`,
:doc:`fix evaporate <fix_evaporate>`,
:doc:`delete_atoms <delete_atoms>`,
:doc:`compute fragment/atom <compute_cluster_atom>`

Default
"""""""

The default keyword settings are compress = yes, single = no, and
offset = -1.
