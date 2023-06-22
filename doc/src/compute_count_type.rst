.. index:: compute count/type

compute count/type command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID count/type mode

* ID, group-ID are documented in :doc:`compute <compute>` command
* count/type = style name of this compute command
* mode = {atom} or {bond} or {angle} or {dihedral} or {improper}

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all count/type atom
   compute 1 flowmols count/type bond

Description
"""""""""""

.. versionadded:: 15Jun2023

Define a computation that counts the current number of atoms for each
atom type.  Or the number of bonds (angles, dihedrals, impropers) for
each bond (angle, dihedral, improper) type.

The former can be useful if atoms are added to or deleted from the
system in random ways, e.g. via the :doc:`fix deposit <fix_deposit>`,
:doc:`fix pour <fix_pour>`, or :doc:`fix evaporate <fix_evaporate>`
commands.  The latter can be useful in reactive simulations where
molecular bonds are broken or created, as well as angles, dihedrals,
impropers.

Note that for this command, bonds (angles, etc) are the topological
kind enumerated in a data file, initially read by the :doc:`read_data
<read_data>` command or defined by the :doc:`molecule <molecule>`
command.  They do not refer to implicit bonds defined on-the-fly by
bond-order or reactive pair styles based on the current conformation
of small clusters of atoms.

These commands can turn off topological bonds (angles, etc) by setting
their bond (angle, etc) types to negative values.  This command
includes the turned-off bonds (angles, etc) in the count for each
type:

* :doc:`fix shake <fix_shake>`
* :doc:`delete_bonds <delete_bonds>`

These commands can create and/or break topological bonds (angles,
etc).  In the case of breaking, they remove the bond (angle, etc) from
the system, so that they no longer exist (:doc:`bond_style quartic
<bond_quartic>` and :doc:`BPM bond styles <Howto_bpm>` are exceptions,
see the discussion below).  Thus they are not included in the counts
for each type:

* :doc:`delete_bonds remove <delete_bonds>`
* :doc:`bond_style quartic <bond_quartic>`
* :doc:`fix bond/react <fix_bond_react>`
* :doc:`fix bond/create <fix_bond_create>`
* :doc:`fix bond/break <fix_bond_break>`
* :doc:`BPM package <Howto_bpm>` bond styles

----------

If the {mode} setting is {atom} then the count of atoms for each atom
type is tallied.  Only atoms in the specified group are counted.

If the {mode} setting is {bond} then the count of bonds for each bond
type is tallied.  Only bonds with both atoms in the specified group
are counted.

For {mode} = {bond}, broken bonds with a bond type of zero are also
counted.  The :doc:`bond_style quartic <bond_quartic>` and :doc:`BPM
bond styles <Howto_bpm>` break bonds by doing this.  See the :doc:`
Howto broken bonds <Howto_broken_bonds>` doc page for more details.
Note that the group setting is ignored for broken bonds; all broken
bonds in the system are counted.

If the {mode} setting is {angle} then the count of angles for each
angle type is tallied.  Only angles with all 3 atoms in the specified
group are counted.

If the {mode} setting is {dihedral} then the count of dihedrals for
each dihedral type is tallied.  Only dihedrals with all 4 atoms in the
specified group are counted.

If the {mode} setting is {improper} then the count of impropers for
each improper type is tallied.  Only impropers with all 4 atoms in the
specified group are counted.

----------

Output info
"""""""""""

This compute calculates a global vector of counts.  If the mode is
{atom} or {bond} or {angle} or {dihedral} or {improper}, then the
vector length is the number of atom types or bond types or angle types
or dihedral types or improper types, respectively.

If the mode is {bond} this compute also calculates a global scalar
which is the number of broken bonds with type = 0, as explained above.

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
