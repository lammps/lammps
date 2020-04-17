.. index:: angle_style table

angle_style table command
=========================

angle_style table/omp command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style table style N

* style = *linear* or *spline* = method of interpolation
* N = use N values in table

Examples
""""""""

.. code-block:: LAMMPS

   angle_style table linear 1000
   angle_coeff 3 file.table ENTRY1

Description
"""""""""""

Style *table* creates interpolation tables of length *N* from angle
potential and derivative values listed in a file(s) as a function of
angle The files are read by the :doc:`angle_coeff <angle_coeff>`
command.

The interpolation tables are created by fitting cubic splines to the
file values and interpolating energy and derivative values at each of
*N* angles.  During a simulation, these tables are used to interpolate
energy and force values on individual atoms as needed.  The
interpolation is done in one of 2 styles: *linear* or *spline*\ .

For the *linear* style, the angle is used to find 2 surrounding table
values from which an energy or its derivative is computed by linear
interpolation.

For the *spline* style, a cubic spline coefficients are computed and
stored at each of the *N* values in the table.  The angle is used to
find the appropriate set of coefficients which are used to evaluate a
cubic polynomial which computes the energy or derivative.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above.

* filename
* keyword

The filename specifies a file containing tabulated energy and
derivative values.  The keyword specifies a section of the file.  The
format of this file is described below.

----------

The format of a tabulated file is as follows (without the
parenthesized comments):

.. parsed-literal::

   # Angle potential for harmonic (one or more comment or blank lines)

   HAM                           (keyword is the first text on line)
   N 181 FP 0 0 EQ 90.0          (N, FP, EQ parameters)
                                 (blank line)
   N 181 FP 0 0                  (N, FP parameters)
   1 0.0 200.5 2.5               (index, angle, energy, derivative)
   2 1.0 198.0 2.5
   ...
   181 180.0 0.0 0.0

A section begins with a non-blank line whose first character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the
:doc:`angle_coeff <angle_coeff>` command.  The next line lists (in any
order) one or more parameters for the table.  Each parameter is a
keyword followed by one or more numeric values.

The parameter "N" is required and its value is the number of table
entries that follow.  Note that this may be different than the *N*
specified in the :doc:`angle_style table <angle_style>` command.  Let
Ntable = *N* in the angle_style command, and Nfile = "N" in the
tabulated file.  What LAMMPS does is a preliminary interpolation by
creating splines using the Nfile tabulated values as nodal points.  It
uses these to interpolate as needed to generate energy and derivative
values at Ntable different points.  The resulting tables of length
Ntable are then used as described above, when computing energy and
force for individual angles and their atoms.  This means that if you
want the interpolation tables of length Ntable to match exactly what
is in the tabulated file (with effectively no preliminary
interpolation), you should set Ntable = Nfile.

The "FP" parameter is optional.  If used, it is followed by two values
fplo and fphi, which are the second derivatives at the innermost and
outermost angle settings.  These values are needed by the spline
construction routines.  If not specified by the "FP" parameter, they
are estimated (less accurately) by the first two and last two
derivative values in the table.

The "EQ" parameter is also optional.  If used, it is followed by a the
equilibrium angle value, which is used, for example, by the :doc:`fix shake <fix_shake>` command.  If not used, the equilibrium angle is
set to 180.0.

Following a blank line, the next N lines list the tabulated values.
On each line, the first value is the index from 1 to N, the second value is
the angle value (in degrees), the third value is the energy (in energy
units), and the fourth is -dE/d(theta) (also in energy units).  The third
term is the energy of the 3-atom configuration for the specified
angle.  The last term is the derivative of the energy with respect to
the angle (in degrees, not radians).  Thus the units of the last term
are still energy, not force.  The angle values must increase from one
line to the next.  The angle values must also begin with 0.0 and end
with 180.0, i.e. span the full range of possible angles.

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

**Restart info:**

This angle style writes the settings for the "angle_style table"
command to :doc:`binary restart files <restart>`, so a angle_style
command does not need to specified in an input script that reads a
restart file.  However, the coefficient information is not stored in
the restart file, since it is tabulated in the potential files.  Thus,
angle_coeff commands do need to be specified in the restart input
script.

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

**Default:** none
