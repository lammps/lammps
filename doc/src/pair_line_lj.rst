.. index:: pair_style line/lj

pair_style line/lj command
==========================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style line/lj cutoff

cutoff = global cutoff for interactions (distance units)

Examples
""""""""


.. code-block:: LAMMPS

   pair_style line/lj 3.0
   pair_coeff * * 1.0 1.0 1.0 0.8 1.12
   pair_coeff 1 2 1.0 2.0 1.0 1.5 1.12 5.0
   pair_coeff 1 2 1.0 0.0 1.0 1.0 2.5

Description
"""""""""""

Style *line/lj* treats particles which are line segments as a set of
small spherical particles that tile the line segment length as
explained below.  Interactions between two line segments, each with N1
and N2 spherical particles, are calculated as the pairwise sum of
N1\*N2 Lennard-Jones interactions.  Interactions between a line segment
with N spherical particles and a point particle are treated as the
pairwise sum of N Lennard-Jones interactions.  See the :doc:`pair_style lj/cut <pair_lj>` doc page for the definition of Lennard-Jones
interactions.

The set of non-overlapping spherical sub-particles that represent a
line segment are generated in the following manner.  Their size is a
function of the line segment length and the specified sub-particle
size for that particle type.  If a line segment has a length L and is
of type I, then the number of spheres N that represent the segment is
calculated as N = L/sizeI, rounded up to an integer value.  Thus if L
is not evenly divisible by sizeI, N is incremented to include one
extra sphere.  The centers of the spheres are spaced equally along the
line segment.  Imagine N+1 equally-space points, which include the 2
end points of the segment.  The sphere centers are halfway between
each pair of points.

The LJ interaction between 2 spheres on different line segments (or a
sphere on a line segment and a point particles) is computed with
sub-particle epsilon, sigma, and cutoff values that are set by the
pair\_coeff command, as described below.  If the distance between the 2
spheres is greater than the sub-particle cutoff, there is no
interaction.  This means that some pairs of sub-particles on 2 line
segments may interact, but others may not.

For purposes of creating the neighbor list for pairs of interacting
line segments or lines/point particles, a regular particle-particle
cutoff is used, as defined by the *cutoff* setting above in the
pair\_style command or overridden with an optional argument in the
pair\_coeff command for a type pair as discussed below.  The distance
between the centers of 2 line segments, or the center of a line
segment and a point particle, must be less than this distance (plus
the neighbor skin; see the :doc:`neighbor <neighbor>` command), for
the pair of particles to be included in the neighbor list.

.. note::

   This means that a too-short value for the *cutoff* setting can
   exclude a pair of particles from the neighbor list even if pairs of
   their sub-particle spheres would interact, based on the sub-particle
   cutoff specified in the pair\_coeff command.  E.g. sub-particles at the
   ends of the line segments that are close to each other.  Which may not
   be what you want, since it means the ends of 2 line segments could
   pass through each other.  It is up to you to specify a *cutoff*
   setting that is consistent with the length of the line segments you
   are using and the sub-particle cutoff settings.

For style *line/lj*\ , the following coefficients must be defined for
each pair of atom types via the :doc:`pair_coeff <pair_coeff>` command
as in the examples above, or in the data file or restart files read by
the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* sizeI (distance units)
* sizeJ (distance units)
* epsilon (energy units)
* sigma (distance units)
* subcutoff (distance units)
* cutoff (distance units)

The *sizeI* and *sizeJ* coefficients are the sub-particle sizes for
line particles of type I and type J.  They are used to define the N
sub-particles per segment as described above.  These coefficients are
actually stored on a per-type basis.  Thus if there are multiple
pair\_coeff commands that involve type I, as either the first or
second atom type, you should use consistent values for sizeI or sizeJ
in all of them.  If you do not do this, the last value specified for
sizeI will apply to all segments of type I.  If typeI or typeJ refers
to point particles, the corresponding sizeI or sizeJ is ignored; it
can be set to 0.0.

The *epsilon*\ , *sigma*\ , and *subcutoff* coefficients are used to
compute an LJ interactions between a pair of sub-particles on 2 line
segments (of type I and J), or between a sub particle/point particle
pair.  As discussed above, the *subcutoff* and *cutoff* params are
different.  The latter is only used for building the neighbor list
when the distance between centers of two line segments or one segment
and a point particle is calculated.

The *cutoff* coefficient is optional.  If not specified, the global
cutoff is used.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, coefficients must be specified.
No default mixing rules are used.

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

Defining particles to be line segments so they participate in
line/line or line/particle interactions requires the use the
:doc:`atom_style line <atom_style>` command.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style tri/lj <pair_tri_lj>`

**Default:** none
