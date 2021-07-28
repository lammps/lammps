.. index:: pair_write

pair_write command
==================

Syntax
""""""

.. code-block:: LAMMPS

   pair_write itype jtype N style inner outer file keyword Qi Qj

* itype,jtype = 2 atom types
* N = # of values
* style = *r* or *rsq* or *bitmap*
* inner,outer = inner and outer cutoff (distance units)
* file = name of file to write values to
* keyword = section name in file for this set of tabulated values
* Qi,Qj = 2 atom charges (charge units) (optional)

Examples
""""""""

.. code-block:: LAMMPS

   pair_write 1 3 500 r 1.0 10.0 table.txt LJ
   pair_write 1 1 1000 rsq 2.0 8.0 table.txt Yukawa_1_1 -0.5 0.5

Description
"""""""""""

Write energy and force values to a file as a function of distance for
the currently defined pair potential.  This is useful for plotting the
potential function or otherwise debugging its values.  If the file
already exists, the table of values is appended to the end of the file
to allow multiple tables of energy and force to be included in one
file.  In case a new file is created, the first line will be a comment
containing a "DATE:" and "UNITS:" tag with the current date and the
current :doc:`units <units>` setting as argument.  For subsequent
invocations of the pair_write command, the current units setting is
compared against the entry in the file, if present, and pair_write
will refuse to add a table if the units are not the same.

The energy and force values are computed at distances from inner to
outer for 2 interacting atoms of type itype and jtype, using the
appropriate :doc:`pair_coeff <pair_coeff>` coefficients.  If the style
is *r*, then N distances are used, evenly spaced in r; if the style is
*rsq*, N distances are used, evenly spaced in r\^2.

For example, for N = 7, style = *r*, inner = 1.0, and outer = 4.0,
values are computed at r = 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0.

If the style is *bitmap*, then 2\^N values are written to the file in a
format and order consistent with how they are read in by the
:doc:`pair_coeff <pair_coeff>` command for pair style *table*\ .  For
reasonable accuracy in a bitmapped table, choose N >= 12, an *inner*
value that is smaller than the distance of closest approach of 2
atoms, and an *outer* value <= cutoff of the potential.

If the pair potential is computed between charged atoms, the charges
of the pair of interacting atoms can optionally be specified.  If not
specified, values of Qi = Qj = 1.0 are used.

The file is written in the format used as input for the
:doc:`pair_style <pair_style>` *table* option with *keyword* as the
section name.  Each line written to the file lists an index number
(1-N), a distance (in distance units), an energy (in energy units),
and a force (in force units).

Restrictions
""""""""""""

All force field coefficients for pair and other kinds of interactions
must be set before this command can be invoked.

Due to how the pairwise force is computed, an inner value > 0.0 must
be specified even if the potential has a finite value at r = 0.0.

For EAM potentials, the pair_write command only tabulates the
pairwise portion of the potential, not the embedding portion.

Related commands
""""""""""""""""

:doc:`pair_style table <pair_table>`,
:doc:`pair_style <pair_style>`, :doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
