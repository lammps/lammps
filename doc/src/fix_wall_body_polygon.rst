.. index:: fix wall/body/polygon

fix wall/body/polygon command
=============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID wall/body/polygon k_n c_n c_t wallstyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/body/polygon = style name of this fix command
* k\_n = normal repulsion strength (force/distance or pressure units)
* c\_n = normal damping coefficient (force/distance or pressure units)
* c\_t = tangential damping coefficient (force/distance or pressure units)
* wallstyle = *xplane* or *yplane* or *zplane* or *zcylinder*
* args = list of arguments for a particular style
  
  .. parsed-literal::
  
       *xplane* or *yplane* args = lo hi
         lo,hi = position of lower and upper plane (distance units), either can be NULL)
       *zcylinder* args = radius
         radius = cylinder radius (distance units)

* zero or more keyword/value pairs may be appended to args
* keyword = *wiggle*
  
  .. parsed-literal::
  
       *wiggle* values = dim amplitude period
         dim = *x* or *y* or *z*
         amplitude = size of oscillation (distance units)
         period = time of oscillation (time units)



Examples
""""""""

fix 1 all wall/body/polygon 1000.0 20.0 5.0 xplane -10.0 10.0

Description
"""""""""""

This fix is for use with 2d models of body particles of style
*rounded/polygon*\ .  It bounds the simulation domain with wall(s).  All
particles in the group interact with the wall when they are close
enough to touch it.  The nature of the interaction between the wall
and the polygon particles is the same as that between the polygon
particles themselves, which is similar to a Hookean potential.  See
the :doc:`Howto body <Howto_body>` doc page for more details on using
body particles.

The parameters *k\_n*, *c\_n*, *c\_t* have the same meaning and units as
those specified with the :doc:`pair\_style body/rounded/polygon <pair_body_rounded_polygon>` command.

The *wallstyle* can be planar or cylindrical.  The 2 planar options
specify a pair of walls in a dimension.  Wall positions are given by
*lo* and *hi*\ .  Either of the values can be specified as NULL if a
single wall is desired.  For a *zcylinder* wallstyle, the cylinder's
axis is at x = y = 0.0, and the radius of the cylinder is specified.

Optionally, the wall can be moving, if the *wiggle* keyword is
appended.

For the *wiggle* keyword, the wall oscillates sinusoidally, similar to
the oscillations of particles which can be specified by the :doc:`fix move <fix_move>` command.  This is useful in packing simulations of
particles.  The arguments to the *wiggle* keyword specify a dimension
for the motion, as well as it's *amplitude* and *period*\ .  Note that
if the dimension is in the plane of the wall, this is effectively a
shearing motion.  If the dimension is perpendicular to the wall, it is
more of a shaking motion.  A *zcylinder* wall can only be wiggled in
the z dimension.

Each timestep, the position of a wiggled wall in the appropriate *dim*
is set according to this equation:


.. parsed-literal::

   position = coord + A - A cos (omega \* delta)

where *coord* is the specified initial position of the wall, *A* is
the *amplitude*\ , *omega* is 2 PI / *period*\ , and *delta* is the time
elapsed since the fix was specified.  The velocity of the wall is set
to the derivative of this expression.

**Restart, fix\_modify, output, run start/stop, minimize info:**

None of the :doc:`fix\_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the BODY package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Any dimension (xy) that has a wall must be non-periodic.

Related commands
""""""""""""""""

:doc:`atom\_style body <atom_style>`, :doc:`pair\_style body/rounded/polygon <pair_body_rounded_polygon>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
