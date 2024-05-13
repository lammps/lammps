.. index:: bond_style table
.. index:: bond_style table/omp

bond_style table command
========================

Accelerator Variants: *table/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style table style N

* style = *linear* or *spline* = method of interpolation
* N = use N values in table

Examples
""""""""

.. code-block:: LAMMPS

   bond_style table linear 1000
   bond_coeff 1 file.table ENTRY1

Description
"""""""""""

Style *table* creates interpolation tables of length *N* from bond
potential and force values listed in a file(s) as a function of bond
length.  The files are read by the :doc:`bond_coeff <bond_coeff>`
command.

The interpolation tables are created by fitting cubic splines to the
file values and interpolating energy and force values at each of *N*
distances.  During a simulation, these tables are used to interpolate
energy and force values as needed.  The interpolation is done in one
of 2 styles: *linear* or *spline*.

For the *linear* style, the bond length is used to find 2 surrounding
table values from which an energy or force is computed by linear
interpolation.

For the *spline* style, a cubic spline coefficients are computed and
stored at each of the *N* values in the table.  The bond length is
used to find the appropriate set of coefficients which are used to
evaluate a cubic polynomial which computes the energy or force.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above.

* filename
* keyword

The filename specifies a file containing tabulated energy and force
values.  The keyword specifies a section of the file.  The format of
this file is described below.

----------

Suitable tables for use with this bond style can be created by LAMMPS
itself from existing bond styles using the :doc:`bond_write
<bond_write>` command.  This can be useful to have a template file for
testing the bond style settings and to build a compatible custom file.
Another option to generate tables is the Python code in the
``tools/tabulate`` folder of the LAMMPS source code distribution.

The format of a tabulated file is as follows (without the
parenthesized comments):

.. parsed-literal::

   # Bond potential for harmonic (one or more comment or blank lines)

   HAM                           (keyword is the first text on line)
   N 101 FP 0 0 EQ 0.5           (N, FP, EQ  parameters)
                                 (blank line)
   1 0.00 338.0000 1352.0000     (index, bond-length, energy, force)
   2 0.01 324.6152 1324.9600
   ...
   101 1.00 338.0000 -1352.0000

A section begins with a non-blank line whose first character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the
:doc:`bond_coeff <bond_coeff>` command.  The next line lists (in any
order) one or more parameters for the table.  Each parameter is a
keyword followed by one or more numeric values.

The parameter "N" is required and its value is the number of table
entries that follow.  Note that this may be different than the *N*
specified in the :doc:`bond_style table <bond_style>` command.  Let
Ntable = *N* in the bond_style command, and Nfile = "N" in the
tabulated file.  What LAMMPS does is a preliminary interpolation by
creating splines using the Nfile tabulated values as nodal points.  It
uses these to interpolate as needed to generate energy and force
values at Ntable different points.  The resulting tables of length
Ntable are then used as described above, when computing energy and
force for individual bond lengths.  This means that if you want the
interpolation tables of length Ntable to match exactly what is in the
tabulated file (with effectively no preliminary interpolation), you
should set Ntable = Nfile.

The "FP" parameter is optional.  If used, it is followed by two values
fplo and fphi, which are the derivatives of the force at the innermost
and outermost bond lengths.  These values are needed by the spline
construction routines.  If not specified by the "FP" parameter, they
are estimated (less accurately) by the first two and last two force
values in the table.

The "EQ" parameter is also optional.  If used, it is followed by a the
equilibrium bond length, which is used, for example, by the :doc:`fix
shake <fix_shake>` command.  If not used, the equilibrium bond length is
to the distance in the table with the lowest potential energy.

Following a blank line, the next N lines list the tabulated values.
On each line, the first value is the index from 1 to N, the second value is
the bond length r (in distance units), the third value is the energy (in
energy units), and the fourth is the force (in force units).  The bond
lengths must range from a LO value to a HI value, and increase from
one line to the next.  If the actual bond length is ever smaller than
the LO value or larger than the HI value, then the calculation is
aborted with an error, so it is advisable to cover the whole range
of possible bond lengths.

Note that one file can contain many sections, each with a tabulated
potential.  LAMMPS reads the file section by section until it finds
one that matches the specified keyword.

----------

.. include:: accel_styles.rst

----------

Restart info
""""""""""""

This bond style writes the settings for the "bond_style table" command
to :doc:`binary restart files <restart>`, so a bond_style command does
not need to specified in an input script that reads a restart file.
However, the coefficient information is not stored in the restart file,
since it is tabulated in the potential files.  Thus, bond_coeff commands
do need to be specified in the restart input script.

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`,
:doc:`bond_write <bond_write>`

Default
"""""""

none
