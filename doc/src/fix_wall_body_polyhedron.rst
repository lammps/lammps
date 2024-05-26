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
* keyword = *wiggle*

  .. parsed-literal::

       *wiggle* values = dim amplitude period
         dim = *x* or *y* or *z*
         amplitude = size of oscillation (distance units)
         period = time of oscillation (time units)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all wall/body/polyhedron 1000.0 20.0 5.0 xplane -10.0 10.0

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

The *wallstyle* can be planar or cylindrical.  The 3 planar options
specify a pair of walls in a dimension.  Wall positions are given by
*lo* and *hi*\ .  Either of the values can be specified as NULL if a
single wall is desired.

Optionally, the wall can be moving, if the *wiggle* keyword is appended.

For the *wiggle* keyword, the wall oscillates sinusoidally, similar to
the oscillations of particles which can be specified by the :doc:`fix move <fix_move>` command.  This is useful in packing simulations of
particles.  The arguments to the *wiggle* keyword specify a dimension
for the motion, as well as it's *amplitude* and *period*\ .  Note that
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

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

none
