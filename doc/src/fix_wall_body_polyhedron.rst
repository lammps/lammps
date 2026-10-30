.. index:: fix wall/body/polyhedron

fix wall/body/polyhedron command
================================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID wall/body/polyhedron k_n c_n c_t wallstyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/body/polyhedron = style name of this fix command
* k_n = normal repulsion strength (force/distance units or pressure units - see discussion below)
* c_n = normal damping coefficient (force/distance units or pressure units - see discussion below)
* c_t = tangential damping coefficient (force/distance units or pressure units - see discussion below)
* wallstyle = *xplane* or *yplane* or *zplane*
* args = list of arguments for a particular style

  .. parsed-literal::

       *xplane* or *yplane* or *zplane* args = lo hi
         lo,hi = position of lower and upper plane (distance units), either can be NULL)


* zero or more keyword/value pairs may be appended to args
* keyword = *wiggle* or *history*

  .. parsed-literal::

       *wiggle* values = dim amplitude period
         dim = *x* or *y* or *z*
         amplitude = size of oscillation (distance units)
         period = time of oscillation (time units)
       *history* values = mu k_t
         mu = friction coefficient between the wall and the particles
         k_t = tangential stiffness (same units as k_n), or NULL for 2/7 k_n

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all wall/body/polyhedron 1000.0 20.0 5.0 xplane -10.0 10.0
   fix 1 all wall/body/polyhedron 1000.0 20.0 5.0 yplane 0.0 NULL history 0.5 NULL

Description
"""""""""""

This fix is for use with 3d models of body particles of style
*rounded/polyhedron*\ .  It bounds the simulation domain with wall(s).
All particles in the group interact with the wall when they are close
enough to touch it.  The nature of the interaction between the wall
and the polygon particles is the same as that between the polygon
particles themselves, which is similar to a Hookean potential.  See
the :doc:`Howto body <Howto_body>` page for more details on using
body particles.

The parameters *k_n*, *c_n*, *c_t* have the same meaning and units as
those specified with the :doc:`pair_style body/rounded/polyhedron <pair_body_rounded_polyhedron>` command.

.. versionchanged:: 30Sep2026

Each vertex of a particle, or the center of a sphere, is repelled
by the wall with the force :math:`k_n (r_v - d)`, when its signed distance
:math:`d` from the wall along the wall normal is smaller than its rounded
radius :math:`r_v`, also when the vertex has moved past the wall.  The
damping forces act at the contact point between the rounded surface and
the wall, using the velocity of the particle at that point relative to the
wall, so that they also exert torques.  There is no cohesion with the
wall, and no friction force unless the *history* keyword is used (see
below), and the contact forces are not scaled by
the size of the contact region.  Previously, a vertex that had moved past
the wall was no longer repelled, or even pushed further out.  A
particle interacts with both walls of a pair, e.g. in a narrow channel.
Previously, only the wall nearer to the center of the particle was
checked, so that an elongated particle could overlap the other wall
without being repelled.

The *wallstyle* is planar.  The 3 options specify a pair of walls in a
dimension.  Wall positions are given by
*lo* and *hi*\ .  Either of the values can be specified as NULL if a
single wall is desired.

Optionally, the wall can be moving, if the *wiggle* keyword is appended.

For the *wiggle* keyword, the wall oscillates sinusoidally, similar to
the oscillations of particles which can be specified by the :doc:`fix move <fix_move>` command.  This is useful in packing simulations of
particles.  The arguments to the *wiggle* keyword specify a dimension
for the motion, as well as its *amplitude* and *period*\ .  Note that
if the dimension is in the plane of the wall, this is effectively a
shearing motion.  If the dimension is perpendicular to the wall, it is
more of a shaking motion.

Each timestep, the position of a wiggled wall in the appropriate *dim*
is set according to this equation:

.. parsed-literal::

   position = coord + A - A cos (omega \* delta)

where *coord* is the specified initial position of the wall, *A* is
the *amplitude*, *omega* is 2 PI / *period*, and *delta* is the time
elapsed since the fix was specified.  The velocity of the wall is set
to the derivative of this expression.

.. versionadded:: 30Sep2026

With the *history* keyword, the wall also exerts a friction force from a
tangential spring on the particles, similar to the *history* keyword of
:doc:`pair_style body/rounded/polyhedron <pair_body_rounded_polyhedron>` between
two particles.  The tangential deformation :math:`\xi` of a particle
touching the wall is accumulated from its tangential velocity relative
to the wall, :math:`\xi \leftarrow \xi + v_t \Delta t`, and the
friction force is :math:`-k_t \xi`, with a magnitude of at most
:math:`\mu F_n`, where :math:`F_n` is the sum of the elastic normal
forces of all vertices of the particle touching the wall.  Once the
spring force exceeds this limit, the particle slides and :math:`\xi` is
reduced accordingly.  The friction force acts at the average of the
contact points of these vertices, weighted by their elastic normal
forces, and :math:`v_t` is the velocity of the particle at that point,
so that a particle resting with a face on the wall is not subject to a
spurious torque.  The tangential deformation is reset to zero when the
particle no longer touches the wall.  A particle touching both walls of
a pair has a separate tangential deformation and friction force for
each wall.  If *k_t* is specified as NULL, it
is set to :math:`\frac{2}{7} k_n` as for the pair style.  The damping
forces are unchanged.  This allows particles to stay at rest on the
wall, e.g. under gravity that is inclined with respect to the wall.

-----------------

Dump image info
"""""""""""""""

.. versionadded:: 11Feb2026

This fix supports the *fix* keyword of :doc:`dump image <dump_image>`.
The fix will pass geometry information about *xplane*\, *yplane*\, and
*zplane* style walls to *dump image* so that the walls will be included
in the rendered image.  Please note, that for :doc:`2d systems
<dimension>`, a wall rendered as a plane would be invisible and it is
thus rendered as a cylinder.

The color of the wall is by default that of the first atom type when
using color styles "type" or "element".  With color style "const" the
default value of "white" can be changed using :doc:`dump_modify fcolor
<dump_image>`.  The transparency is by default fully opaque and can be
changed globally with *dump\_modify ftrans*\ .

For 2d systems, the *fflag1* setting determines whether the cylinder
representing the wall is capped with a sphere at the ends: 0 means no caps, 1
means the lower end is capped, 2 means the upper end is capped, and 3
means both ends are capped.  The *fflag2* setting allows to set the
radius of the rendered cylinders.

For 3d systems, both *fflag1* and *fflag2* are ignored.

------------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With the *history* keyword, this fix writes the tangential deformations
of the particles touching the wall to :doc:`binary restart files
<restart>`, so that a simulation can continue correctly.  See the
:doc:`read_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file, so that the
operation of the fix continues in an uninterrupted fashion.  Otherwise
no information about this fix is written to binary restart files.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BODY package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Any dimension (xyz) that has a wall must be non-periodic.

Related commands
""""""""""""""""

:doc:`atom_style body <atom_style>`, :doc:`pair_style body/rounded/polyhedron <pair_body_rounded_polyhedron>`

Default
"""""""

The *history* keyword is not used.
