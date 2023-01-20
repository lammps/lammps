.. index:: angle_write

angle_write command
===================

Syntax
""""""

.. code-block:: LAMMPS

   angle_write atype N file keyword

* atype = angle type
* N = # of values
* file = name of file to write values to
* keyword = section name in file for this set of tabulated values

Examples
""""""""

.. code-block:: LAMMPS

   angle_write 1  500 table.txt Harmonic_1
   angle_write 3 1000 table.txt Harmonic_3

Description
"""""""""""

.. versionadded:: TBD

Write energy and force values to a file as a function of angle for the
currently defined angle potential.  Force in this context means the
force with respect to the angle, not the force on individual atoms.
This is useful for plotting the potential function or otherwise
debugging its values.  The resulting file can also be used as input for
use with :doc:`angle style table <angle_table>`.

If the file already exists, the table of values is appended to the end
of the file to allow multiple tables of energy and force to be included
in one file.  The individual sections may be identified by the *keyword*.

The energy and force values are computed for angles ranging from 0
degrees to 180 degrees for 3 interacting atoms forming an angle type
atype, using the appropriate :doc:`angle_coeff <angle_coeff>`
coefficients. N evenly spaced angles are used.

For example, for N = 6, values are computed at :math:`\theta = 0, 36,
72, 108, 144, 180`.

The file is written in the format used as input for the
:doc:`angle_style table <angle_table>` option with *keyword* as the
section name.  Each line written to the file lists an index number
(1-N), an angle (in degrees), an energy (in energy units), and a force
(in force units per radians^2).  In case a new file is created, the
first line will be a comment with a "DATE:" and "UNITS:" tag with the
current date and :doc:`units <units>` settings.  For subsequent
invocations of the *angle_write* command for the same file, data will be
appended and the current units settings will be compared to the data
from the header, if present. The *angle_write* will refuse to add a
table to an existing file if the units are not the same.

Restrictions
""""""""""""

All force field coefficients for angle and other kinds of interactions
must be set before this command can be invoked.

The table of the angle energy and force data data is created by using a
separate, internally created, new LAMMPS instance with a dummy system of
3 atoms for which the angle potential energy is computed after
transferring the angle style and coefficients and arranging the 3 atoms
into the corresponding geometries.  The angle force is then determined
from the potential energies through numerical differentiation.  As a
consequence of this approach, not all angle styles are compatible. The
following conditions must be met:

- The angle style must be able to write its coefficients to a data file.
  This condition excludes for example :doc:`angle style hybrid <angle_hybrid>` and
  :doc:`angle style table <angle_table>`.
- The potential function must not have any terms that depend on geometry
  properties other than the angle. This condition excludes for example
  :doc:`angle style class2 <angle_class2>` all angle types for
  :doc:`angle style charmm <angle_charmm>` that have non-zero
  Urey-Bradley terms.  Please note that the *write_angle* command has no
  way of checking for this condition, so the resulting tables may be
  bogus if the requirement is not met.  It is thus recommended to make
  careful tests for any created tables.

Related commands
""""""""""""""""

:doc:`angle_style table <angle_table>`, :doc:`bond_write <bond_write>`,
:doc:`dihedral_write <dihedral_write>`, :doc:`angle_style <angle_style>`,
:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none
