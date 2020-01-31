.. index:: pair\_style table/rx

pair\_style table/rx command
============================

pair\_style table/rx/kk command
===============================

Syntax
""""""


.. parsed-literal::

   pair_style table style N ...

* style = *lookup* or *linear* or *spline* or *bitmap* = method of interpolation
* N = use N values in *lookup*\ , *linear*\ , *spline* tables
* weighting = fractional or molecular (optional)

Examples
""""""""


.. parsed-literal::

   pair_style table/rx linear 1000
   pair_style table/rx linear 1000 fractional
   pair_style table/rx linear 1000 molecular
   pair_coeff \* \* rxn.table ENTRY1 h2o h2o 10.0
   pair_coeff \* \* rxn.table ENTRY1 1fluid 1fluid 10.0
   pair_coeff \* 3 rxn.table ENTRY1 h2o no2 10.0

Description
"""""""""""

Style *table/rx* is used in reaction DPD simulations,where the
coarse-grained (CG) particles are composed of *m* species whose
reaction rate kinetics are determined from a set of *n* reaction rate
equations through the :doc:`fix rx <fix_rx>` command.  The species of
one CG particle can interact with a species in a neighboring CG
particle through a site-site interaction potential model. Style
*table/rx* creates interpolation tables of length *N* from pair
potential and force values listed in a file(s) as a function of
distance.  The files are read by the :doc:`pair_coeff <pair_coeff>`
command.

The interpolation tables are created by fitting cubic splines to the
file values and interpolating energy and force values at each of *N*
distances.  During a simulation, these tables are used to interpolate
energy and force values as needed.  The interpolation is done in one
of 4 styles: *lookup*\ , *linear*\ , *spline*\ , or *bitmap*\ .

For the *lookup* style, the distance between 2 atoms is used to find
the nearest table entry, which is the energy or force.

For the *linear* style, the pair distance is used to find 2
surrounding table values from which an energy or force is computed by
linear interpolation.

For the *spline* style, a cubic spline coefficients are computed and
stored at each of the *N* values in the table.  The pair distance is
used to find the appropriate set of coefficients which are used to
evaluate a cubic polynomial which computes the energy or force.

For the *bitmap* style, the N means to create interpolation tables
that are 2\^N in length.  The pair distance is used to index into the
table via a fast bit-mapping technique :ref:`(Wolff) <Wolff>` and a linear
interpolation is performed between adjacent table values.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above.

* filename
* keyword
* species1
* species2
* cutoff (distance units)

The filename specifies a file containing tabulated energy and force
values.  The keyword specifies a section of the file.  The cutoff is
an optional coefficient.  If not specified, the outer cutoff in the
table itself (see below) will be used to build an interpolation table
that extend to the largest tabulated distance.  If specified, only
file values up to the cutoff are used to create the interpolation
table.  The format of this file is described below.

The species tags define the site-site interaction potential between
two species contained within two different particles.  The species
tags must either correspond to the species defined in the reaction
kinetics files specified with the :doc:`fix rx <fix_rx>` command or they
must correspond to the tag "1fluid", signifying interaction with a
product species mixture determined through a one-fluid approximation.
The interaction potential is weighted by the geometric average of
either the mole fraction concentrations or the number of molecules
associated with the interacting coarse-grained particles (see the
*fractional* or *molecular* weighting pair style options). The coarse-grained potential is
stored before and after the reaction kinetics solver is applied, where
the difference is defined to be the internal chemical energy (uChem).


----------


Here are some guidelines for using the pair\_style table/rx command to
best effect:

* Vary the number of table points; you may need to use more than you think
  to get good resolution.
* Always use the :doc:`pair_write <pair_write>` command to produce a plot
  of what the final interpolated potential looks like.  This can show up
  interpolation "features" you may not like.
* Start with the linear style; it's the style least likely to have problems.
* Use *N* in the pair\_style command equal to the "N" in the tabulation
  file, and use the "RSQ" or "BITMAP" parameter, so additional interpolation
  is not needed.  See discussion below.
* Make sure that your tabulated forces and tabulated energies are consistent
  (dE/dr = -F) along the entire range of r values.
* Use as large an inner cutoff as possible.  This avoids fitting splines
  to very steep parts of the potential.



----------


The format of a tabulated file is a series of one or more sections,
defined as follows (without the parenthesized comments):


.. parsed-literal::

   # Morse potential for Fe   (one or more comment or blank lines)

   MORSE_FE                   (keyword is first text on line)
   N 500 R 1.0 10.0           (N, R, RSQ, BITMAP, FPRIME parameters)
                              (blank)
   1 1.0 25.5 102.34          (index, r, energy, force)
   2 1.02 23.4 98.5
   ...
   500 10.0 0.001 0.003

A section begins with a non-blank line whose 1st character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the pair\_coeff
command.  The next line lists (in any order) one or more parameters
for the table.  Each parameter is a keyword followed by one or more
numeric values.

The parameter "N" is required and its value is the number of table
entries that follow.  Note that this may be different than the *N*
specified in the :doc:`pair_style table/rx <pair_style>` command.  Let
Ntable = *N* in the pair\_style command, and Nfile = "N" in the
tabulated file.  What LAMMPS does is a preliminary interpolation by
creating splines using the Nfile tabulated values as nodal points.  It
uses these to interpolate as needed to generate energy and force
values at Ntable different points.  The resulting tables of length
Ntable are then used as described above, when computing energy and
force for individual pair distances.  This means that if you want the
interpolation tables of length Ntable to match exactly what is in the
tabulated file (with effectively no preliminary interpolation), you
should set Ntable = Nfile, and use the "RSQ" or "BITMAP" parameter.
The internal table abscissa is RSQ (separation distance squared).

All other parameters are optional.  If "R" or "RSQ" or "BITMAP" does
not appear, then the distances in each line of the table are used
as-is to perform spline interpolation.  In this case, the table values
can be spaced in *r* uniformly or however you wish to position table
values in regions of large gradients.

If used, the parameters "R" or "RSQ" are followed by 2 values *rlo*
and *rhi*\ .  If specified, the distance associated with each energy and
force value is computed from these 2 values (at high accuracy), rather
than using the (low-accuracy) value listed in each line of the table.
The distance values in the table file are ignored in this case.
For "R", distances uniformly spaced between *rlo* and *rhi* are
computed; for "RSQ", squared distances uniformly spaced between
*rlo\*rlo* and *rhi\*rhi* are computed.

If used, the parameter "BITMAP" is also followed by 2 values *rlo* and
*rhi*\ .  These values, along with the "N" value determine the ordering
of the N lines that follow and what distance is associated with each.
This ordering is complex, so it is not documented here, since this
file is typically produced by the :doc:`pair_write <pair_write>` command
with its *bitmap* option.  When the table is in BITMAP format, the "N"
parameter in the file must be equal to 2\^M where M is the value
specified in the pair\_style command.  Also, a cutoff parameter cannot
be used as an optional 3rd argument in the pair\_coeff command; the
entire table extent as specified in the file must be used.

If used, the parameter "FPRIME" is followed by 2 values *fplo* and
*fphi* which are the derivative of the force at the innermost and
outermost distances listed in the table.  These values are needed by
the spline construction routines.  If not specified by the "FPRIME"
parameter, they are estimated (less accurately) by the first 2 and
last 2 force values in the table.  This parameter is not used by
BITMAP tables.

Following a blank line, the next N lines list the tabulated values.
On each line, the 1st value is the index from 1 to N, the 2nd value is
r (in distance units), the 3rd value is the energy (in energy units),
and the 4th is the force (in force units).  The r values must increase
from one line to the next (unless the BITMAP parameter is specified).

Note that one file can contain many sections, each with a tabulated
potential.  LAMMPS reads the file section by section until it finds
one that matches the specified keyword.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This pair style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

The :doc:`pair_modify <pair_modify>` shift, table, and tail options are
not relevant for this pair style.

This pair style writes the settings for the "pair\_style table/rx" command
to :doc:`binary restart files <restart>`, so a pair\_style command does
not need to specified in an input script that reads a restart file.
However, the coefficient information is not stored in the restart
file, since it is tabulated in the potential files.  Thus, pair\_coeff
commands do need to be specified in the restart input script.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


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


This command is part of the USER-DPD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** fractional weighting


----------


.. _Wolff:



**(Wolff)** Wolff and Rudd, Comp Phys Comm, 120, 200-32 (1999).
