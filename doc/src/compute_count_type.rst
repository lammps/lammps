.. index:: compute count/type

compute count/type command
====================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID count/type mode

* ID, group-ID are documented in :doc:`compute <compute>` command
* count/type = style name of this compute command
* mode = {atom} or {bond}
  
Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all count/type atom
   compute 1 flowmols count/type bond

Description
"""""""""""

Define a computation that counts the current number of atoms by atom
type or the number of bonds by bond type.  The latter can be useful in
reactive simulations where bonds are broken or created.

Note that for this command, bonds are the topological ones enumerated
in a data file, initially read by the :doc:`read_data <read_data>`
command.  They do not refer to bonds defined on-the-fly by bond-order
or reactive pair styles.

These commands can create and break toplogical bonds:

* :doc:`fix bond/react <fix_bond_react>`
* :doc:`fix bond/create <fix_bond_create>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`bond_style quartic <bond_quartic>`
* :doc:`BPM package <Howto_bpm>` bond styles

If the {mode} setting is {atom} then the count of atoms for each atom
type is tallied.  Only atoms in the specified group are counted.

If the {mode} setting is {bond} then the count of bonds for each bond
type is tallied.  Only bonds with both atoms in the specified group
are counted.

For {mode} = {bond}, broken bonds with a bond type of zero are also
counted.  See the :doc:`Howto broken bonds <Howto_broken_bonds>` doc
page for details.  Note that the group setting is ignored for broken
bonds; all broken bonds in the system are counted.

----------

Output info
"""""""""""

This compute calculates a global vector of counts.  If the mode is
{atom}, the vector length is the number of atom types.  If the mode is
{bond}, the vector length is the number of bond types.

If the mode is {bond} this compute also calculates a global scalar
which is the number of broken bonds.

These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output
<Howto_output>` page for an overview of LAMMPS output options.

The scalar and vector values calculated by this compute are "extensive".

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

none

Default
"""""""

none
