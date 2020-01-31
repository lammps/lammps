.. index:: pair\_style tri/lj

pair\_style tri/lj command
==========================

Syntax
""""""


.. parsed-literal::

   pair_style tri/lj cutoff

cutoff = global cutoff for interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style tri/lj 3.0
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 1.0 1.5 2.5

Description
"""""""""""

Style *tri/lj* treats particles which are triangles as a set of small
spherical particles that tile the triangle surface as explained below.
Interactions between two triangles, each with N1 and N2 spherical
particles, are calculated as the pairwise sum of N1\*N2 Lennard-Jones
interactions.  Interactions between a triangle with N spherical
particles and a point particle are treated as the pairwise sum of N
Lennard-Jones interactions.  See the :doc:`pair_style lj/cut <pair_lj>`
doc page for the definition of Lennard-Jones interactions.

The cutoff distance for an interaction between 2 triangles, or between
a triangle and a point particle, is calculated from the position of
the triangle (its centroid), not between pairs of individual spheres
comprising the triangle.  Thus an interaction is either calculated in
its entirety or not at all.

The set of non-overlapping spherical particles that represent a
triangle, for purposes of this pair style, are generated in the
following manner.  Assume the triangle is of type I, and sigma\_II has
been specified.  We want a set of spheres with centers in the plane of
the triangle, none of them larger in diameter than sigma\_II, which
completely cover the triangle's area, but with minimal overlap and a
minimal total number of spheres.  This is done in a recursive manner.
Place a sphere at the centroid of the original triangle.  Calculate
what diameter it must have to just cover all 3 corner points of the
triangle.  If that diameter is equal to or smaller than sigma\_II, then
include a sphere of the calculated diameter in the set of covering
spheres.  It the diameter is larger than sigma\_II, then split the
triangle into 2 triangles by bisecting its longest side.  Repeat the
process on each sub-triangle, recursing as far as needed to generate a
set of covering spheres.  When finished, the original criteria are
met, and the set of covering spheres should be near minimal in number
and overlap, at least for input triangles with a reasonable
aspect-ratio.

The LJ interaction between 2 spheres on different triangles of types
I,J is computed with an arithmetic mixing of the sigma values of the 2
spheres and using the specified epsilon value for I,J atom types.
Note that because the sigma values for triangles spheres is computed
using only sigma\_II values, specific to the triangles's type, this
means that any specified sigma\_IJ values (for I != J) are effectively
ignored.

For style *tri/lj*\ , the following coefficients must be defined for
each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command
as in the examples above, or in the data file or restart files read by
the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* epsilon (energy units)
* sigma (distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global cutoff
is used.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of this pair style can be mixed.  The
default mix value is *geometric*\ .  See the "pair\_modify" command for
details.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


This style is part of the ASPHERE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Defining particles to be triangles so they participate in tri/tri or
tri/particle interactions requires the use the :doc:`atom_style tri <atom_style>` command.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style line/lj <pair_line_lj>`

**Default:** none


