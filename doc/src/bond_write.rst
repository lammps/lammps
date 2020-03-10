.. index:: bond_write

bond_write command
==================

Syntax
""""""

.. code-block:: LAMMPS

   bond_write btype N inner outer file keyword itype jtype

* btype = bond types
* N = # of values
* inner,outer = inner and outer bond length (distance units)
* file = name of file to write values to
* keyword = section name in file for this set of tabulated values
* itype,jtype = 2 atom types (optional)

Examples
""""""""

.. code-block:: LAMMPS

   bond_write 1 500 0.5 3.5 table.txt Harmonic_1
   bond_write 3 1000 0.1 6.0 table.txt Morse

Description
"""""""""""

Write energy and force values to a file as a function of distance for
the currently defined bond potential.  This is useful for plotting the
potential function or otherwise debugging its values.  If the file
already exists, the table of values is appended to the end of the file
to allow multiple tables of energy and force to be included in one
file.

The energy and force values are computed at distances from inner to
outer for 2 interacting atoms forming a bond of type btype, using the
appropriate :doc:`bond_coeff <bond_coeff>` coefficients. N evenly spaced
distances are used.

For example, for N = 7, inner = 1.0, and outer = 4.0,
values are computed at r = 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0.

The file is written in the format used as input for the
:doc:`bond_style <bond_style>` *table* option with *keyword* as the
section name.  Each line written to the file lists an index number
(1-N), a distance (in distance units), an energy (in energy units),
and a force (in force units).

Restrictions
""""""""""""

All force field coefficients for bond and other kinds of interactions
must be set before this command can be invoked.

Due to how the bond force is computed, an inner value > 0.0 must
be specified even if the potential has a finite value at r = 0.0.

Related commands
""""""""""""""""

:doc:`bond_style table <bond_table>`,
:doc:`bond_style <bond_style>`, :doc:`bond_coeff <bond_coeff>`

**Default:** none
