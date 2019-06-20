.. index:: bond\_style table

bond\_style table command
=========================

bond\_style table/omp command
=============================

Syntax
""""""


.. parsed-literal::

   bond_style table style N

* style = *linear* or *spline* = method of interpolation
* N = use N values in table

Examples
""""""""


.. parsed-literal::

   bond_style table linear 1000
   bond_coeff 1 file.table ENTRY1

Description
"""""""""""

Style *table* creates interpolation tables of length *N* from bond
potential and force values listed in a file(s) as a function of bond
length.  The files are read by the :doc:`bond\_coeff <bond_coeff>`
command.

The interpolation tables are created by fitting cubic splines to the
file values and interpolating energy and force values at each of *N*
distances.  During a simulation, these tables are used to interpolate
energy and force values as needed.  The interpolation is done in one
of 2 styles: *linear* or *spline*\ .

For the *linear* style, the bond length is used to find 2 surrounding
table values from which an energy or force is computed by linear
interpolation.

For the *spline* style, a cubic spline coefficients are computed and
stored at each of the *N* values in the table.  The bond length is
used to find the appropriate set of coefficients which are used to
evaluate a cubic polynomial which computes the energy or force.

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above.

* filename
* keyword

The filename specifies a file containing tabulated energy and force
values.  The keyword specifies a section of the file.  The format of
this file is described below.


----------


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

A section begins with a non-blank line whose 1st character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the
:doc:`bond\_coeff <bond_coeff>` command.  The next line lists (in any
order) one or more parameters for the table.  Each parameter is a
keyword followed by one or more numeric values.

The parameter "N" is required and its value is the number of table
entries that follow.  Note that this may be different than the *N*
specified in the :doc:`bond\_style table <bond_style>` command.  Let
Ntable = *N* in the bond\_style command, and Nfile = "N" in the
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
equilibrium bond length, which is used, for example, by the :doc:`fix shake <fix_shake>` command.  If not used, the equilibrium bond
length is to the distance in the table with the lowest potential energy.

Following a blank line, the next N lines list the tabulated values.
On each line, the 1st value is the index from 1 to N, the 2nd value is
the bond length r (in distance units), the 3rd value is the energy (in
energy units), and the 4th is the force (in force units).  The bond
lengths must range from a LO value to a HI value, and increase from
one line to the next.  If the actual bond length is ever smaller than
the LO value or larger than the HI value, then the bond energy and
force is evaluated as if the bond were the LO or HI length.

Note that one file can contain many sections, each with a tabulated
potential.  LAMMPS reads the file section by section until it finds
one that matches the specified keyword.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`, :doc:`delete\_bonds <delete_bonds>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
