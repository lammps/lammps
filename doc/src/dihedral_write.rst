.. index:: dihedral_write

dihedral_write command
======================

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_write dtype N file keyword

* dtype = dihedral type
* N = # of values
* file = name of file to write values to
* keyword = section name in file for this set of tabulated values

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_write 1  500 table.txt Harmonic_1
   dihedral_write 3 1000 table.txt Harmonic_3

Description
"""""""""""

.. versionadded:: 8Feb2023

Write energy and force values to a file as a function of the dihedral
angle for the currently defined dihedral potential.  Force in this
context means the force with respect to the dihedral angle, not the
force on individual atoms.  This is useful for plotting the potential
function or otherwise debugging its values.  The resulting file can also
be used as input for use with :doc:`dihedral style table
<dihedral_table>`.

If the file already exists, the table of values is appended to the end
of the file to allow multiple tables of energy and force to be included
in one file.  The individual sections may be identified by the *keyword*.

The energy and force values are computed for dihedrals ranging from 0
degrees to 360 degrees for 4 interacting atoms forming an dihedral type
dtype, using the appropriate :doc:`dihedral_coeff <dihedral_coeff>`
coefficients. N evenly spaced dihedrals are used.  Since 0 and 360
degrees are the same dihedral angle, the latter entry is skipped.

For example, for N = 6, values would be computed at
:math:`\phi = 0, 60, 120, 180, 240, 300`.

The file is written in the format used as input for the
:doc:`dihedral_style table <dihedral_table>` option with *keyword* as
the section name.  Each line written to the file lists an index number
(1-N), an dihedral angle (in degrees), an energy (in energy units), and
a force (in force units per radians^2).  In case a new file is created,
the first line will be a comment with a "DATE:" and "UNITS:" tag with
the current date and :doc:`units <units>` settings.  For subsequent
invocations of the *dihedral_write* command for the same file, data will
be appended and the current units settings will be compared to the data
from the header, if present. The *dihedral_write* will refuse to add a
table to an existing file if the units are not the same.

Restrictions
""""""""""""

All force field coefficients for dihedrals and other kinds of interactions
must be set before this command can be invoked.

The table of the dihedral energy and force data data is created by using a
separate, internally created, new LAMMPS instance with a dummy system of
4 atoms for which the dihedral potential energy is computed after
transferring the dihedral style and coefficients and arranging the 4 atoms
into the corresponding geometries.  The dihedral force is then determined
from the potential energies through numerical differentiation.  As a
consequence of this approach, not all dihedral styles are compatible. The
following conditions must be met:

- The dihedral style must be able to write its coefficients to a data file.
  This condition excludes for example :doc:`dihedral style hybrid <dihedral_hybrid>` and
  :doc:`dihedral style table <dihedral_table>`.
- The potential function must not have any terms that depend on geometry
  properties other than the dihedral.  This condition excludes for
  example :doc:`dihedral style class2 <dihedral_class2>`.  Please note
  that the *write_dihedral* command has no way of checking for this
  condition.  It will check the style name against an internal list of
  known to be incompatible styles.  The resulting tables may be bogus
  for unlisted dihedral styles if the requirement is not met.  It is
  thus recommended to make careful tests for any created tables.

Related commands
""""""""""""""""

:doc:`dihedral_style table <dihedral_table>`, :doc:`bond_write <bond_write>`,
:doc:`angle_write <angle_write>`, :doc:`dihedral_style <dihedral_style>`,
:doc:`dihedral_coeff <dihedral_coeff>`

Default
"""""""

none
