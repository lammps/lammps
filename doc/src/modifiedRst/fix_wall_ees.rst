.. index:: fix wall/ees

fix wall/ees command
====================

fix wall/region/ees command
===========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID style args

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *wall/ees* or *wall/region/ees*
  
  .. parsed-literal::
  
       args for style *wall/ees*\ : one or more *face parameters* groups may be appended
       face = *xlo* or *xhi* or *ylo* or *yhi* or *zlo* or *zhi*
       parameters = coord epsilon sigma cutoff
         coord = position of wall = EDGE or constant or variable
           EDGE = current lo or hi edge of simulation box
           constant = number like 0.0 or -30.0 (distance units)
           variable = :doc:`equal-style variable <variable>` like v_x or v_wiggle
         epsilon = strength factor for wall-particle interaction (energy or energy/distance\^2 units)
           epsilon can be a variable (see below)
         sigma = size factor for wall-particle interaction (distance units)
           sigma can be a variable (see below)
         cutoff = distance from wall at which wall-particle interaction is cut off (distance units)

  
  .. parsed-literal::
  
       args for style *wall/region/ees*\ : *region-ID* *epsilon* *sigma* *cutoff*
         region-ID = region whose boundary will act as wall
         epsilon = strength factor for wall-particle interaction (energy or energy/distance\^2 units)
         sigma = size factor for wall-particle interaction (distance units)
         cutoff = distance from wall at which wall-particle interaction is cut off (distance units)



Examples
""""""""


.. parsed-literal::

   fix wallhi all wall/ees xlo -1.0 1.0 1.0 2.5 units box
   fix wallhi all wall/ees xhi EDGE 1.0 1.0 2.5
   fix wallhi all wall/ees v_wiggle 23.2 1.0 1.0 2.5
   fix zwalls all wall/ees zlo 0.0 1.0 1.0 0.858 zhi 40.0 1.0 1.0 0.858

   fix ees_cube all wall/region/ees myCube 1.0 1.0 2.5

Description
"""""""""""

Fix *wall/ees* bounds the simulation domain on one or more of its
faces with a flat wall that interacts with the ellipsoidal atoms in the
group by generating a force on the atom in a direction perpendicular to
the wall and a torque parallel with the wall.  The energy of
wall-particle interactions E is given by:

.. math::

E = \epsilon \left[ \frac{2  \sigma_{LJ}^{12} \left(7 r^5+14 r^3 \sigma_{n}^2+3 r \sigma_{n}^4\right) }{945 \left(r^2-\sigma_{n}^2\right)^7} -\frac{ \sigma_{LJ}^6 \left(2 r \sigma_{n}^3+\sigma_{n}^2 \left(r^2-\sigma_{n}^2\right)\log{ \left[\frac{r-\sigma_{n}}{r+\sigma_{n}}\right]}\right) }{12 \sigma_{n}^5 \left(r^2-\sigma_{n}^2\right)} \right]\qquad \sigma_n < r < r_c


Introduced by Babadi and Ejtehadi in :ref:`(Babadi) <BabadiEjtehadi>`. Here,
*r* is the distance from the particle to the wall at position *coord*\ ,
and Rc is the *cutoff* distance at which the particle and wall no
longer interact. Also, sigma\_n is the distance between center of
ellipsoid and the nearest point of its surface to the wall. The energy
of the wall is:

.. image:: JPG/fix_wall_ees_image.jpg
   :align: center

Details of using this command and specifications are the same as
fix/wall command. You can also find an example in USER/ees/ under
examples/ directory.

The prefactor *epsilon* can be thought of as an
effective Hamaker constant with energy units for the strength of the
ellipsoid-wall interaction.  More specifically, the *epsilon* pre-factor
= 8 \* pi\^2 \* rho\_wall \* rho\_ellipsoid \* epsilon
\* sigma\_a \* sigma\_b \* sigma\_c, where epsilon is the LJ parameters for
the constituent LJ particles and sigma\_a, sigma\_b, and sigma\_c are radii
of ellipsoidal particles. Rho\_wall and rho\_ellipsoid are the number
density of the constituent particles, in the wall and ellipsoid
respectively, in units of 1/volume.

.. note::

   You must insure that r is always bigger than sigma\_n for
   all particles in the group, or LAMMPS will generate an error.  This
   means you cannot start your simulation with particles touching the wall
   position *coord* (r = sigma\_n) or with particles penetrating the wall
   (0 =< r < sigma\_n) or with particles on the wrong side of the
   wall (r < 0).

Fix *wall/region/ees* treats the surface of the geometric region defined
by the *region-ID* as a bounding wall which interacts with nearby
ellipsoidal particles according to the EES potential introduced above.

Other details of this command are the same as for the :doc:`fix wall/region <fix_wall_region>` command.  One may also find an example
of using this fix in the examples/USER/misc/ees/ directory.

Restrictions
""""""""""""


This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix requires that atoms be ellipsoids as defined by the
:doc:`atom\_style ellipsoid <atom_style>` command.

Related commands
""""""""""""""""

:doc:`fix wall <fix_wall>`,
:doc:`pair resquared <pair_resquared>`

Default
"""""""

none


----------


.. _BabadiEjtehadi:



**(Babadi)** Babadi and Ejtehadi, EPL, 77 (2007) 23002.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
