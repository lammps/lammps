.. index:: fix wall/reflect/stochastic

fix wall/reflect/stochastic command
===================================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID wall/reflect/stochastic rstyle seed face args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/reflect/stochastic = style name of this fix command
* rstyle = diffusive or maxwell or ccl
* seed = random seed for stochasticity (positive integer)
* one or more face/args pairs may be appended
* face = *xlo* or *xhi* or *ylo* or *yhi* or *zlo* or *zhi*
  
  .. parsed-literal::
  
       args = pos temp velx vely velz accomx accomy accomz
         pos = EDGE or constant
           EDGE = current lo or hi edge of simulation box
           constant = number like 0.0 or 30.0 (distance units)
         temp = wall temperature (temperature units)
         velx,vely,velz = wall velocity in x,y,z directions (velocity units)
         accomx,accomy,accomz = accommodation coeffs in x,y,z directions (unitless)
           not specified for rstyle = diffusive
           single accom coeff specified for rstyle maxwell
           all 3 coeffs specified for rstyle cll

* zero or more keyword/value pairs may be appended
* keyword = *units*
  
  .. parsed-literal::
  
       *units* value = *lattice* or *box*
         *lattice* = the wall position is defined in lattice units
         *box* = the wall position is defined in simulation box units



Examples
""""""""


.. parsed-literal::

   fix zwalls all wall/reflect/stochastic diffusive 23424 zlo EDGE 300 0.1 0.1 0 zhi EDGE 200 0.1 0.1 0
   fix ywalls all wall/reflect/stochastic maxwell 345533 ylo 5.0 300 0.1 0.0 0.0 0.8 yhi 10.0 300 0.1 0.0 0.0 0.8
   fix xwalls all wall/reflect/stochastic cercignanilampis 2308 xlo 0.0 300 0.0 0.1 0.9 0.8 0.7 xhi EDGE 300 0.0 0.1 0 0.9 0.8 0.7 units box

Description
"""""""""""

Bound the simulation with one or more walls which reflect particles
in the specified group when they attempt to move through them.

Reflection means that if an atom moves outside the wall on a timestep
(e.g. due to the :doc:`fix nve <fix_nve>` command), then it is put back
inside the wall with a changed velocity.

This fix models treats the wall as a moving solid boundary with a
finite temperature, which can exchange energy with particles that
collide with it.  This is different than the simpler :doc:`fix wall/reflect <fix_wall_reflect>` command which models mirror
reflection.  For this fix, the post collision velocity of each
particle is treated stochastically.  The randomness can come from many
sources: thermal motion of the wall atoms, surface roughness, etc.
Three stochastic reflection models are currently implemented.

For rstyle *diffusive*\ , particles are reflected diffusively. Their
velocity distribution corresponds to an equilibrium distribution of
particles at the wall temperature.  No accommodation coefficients
are specified.

For rstyle *maxwell*\ , particle reflection is Maxwellian which means
partially diffusive and partially specular (:ref:`Maxwell <Maxwell>`).  A
single accommodation coeff is specified which must be between 0.0 and
1.0 inclusive.  It determines the fraction of the collision which is
diffusive versus specular.  An accommodation coefficient of 1.0 is fully
diffusive; a coefficient of 0.0 is fully specular.

For rstyle *cll*\ , particle collisions are computed by the
Cercignani/Lampis model.  See :ref:`CL <CL>` and :ref:`To <To>` for details.
Three accommodations coefficient are specified.  Each must be between
0.0 and 1.0 inclusive.  Two are velocity accommodation coefficients;
one is a normal kinetic energy accommodation.  The normal coeff is the
one corresponding to the normal of the wall itself.  For example if
the wall is *ylo* or *yhi*\ , *accomx* and *accomz* are the tangential
velocity accommodation coefficients, and *accomy* is the normal
kinetic energy accommodation coefficient.

The optional *units* keyword determines the distance units used to
define a wall position.  A *box* value selects standard distance units
as defined by the :doc:`units <units>` command, e.g. Angstroms for units
= real or metal.  A *lattice* value means the distance units are in
lattice spacings. The :doc:`lattice <lattice>` command must have been
previously used to define the lattice spacings.


----------


Restrictions
""""""""""""


This fix has the same limitations as the :doc:`fix wall/reflect <fix_wall_reflect>` command.  Any dimension (xyz) that
has a wall must be non-periodic.  It should not be used with rigid
bodies such as those defined by the :doc:`fix rigid <fix_rigid>`
command.  The wall velocity must lie on the same plane as the wall
itself.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix wall/reflect <fix_wall_reflect>`

Default
"""""""

The default for the units keyword is lattice.


----------


.. _Maxwell:



**(Maxwell)** J.C. Maxwell, Philos. Tans. Royal Soc. London, 157: 49-88
(1867).

.. _CL:



**(Cercignani)** C. Cercignani and M. Lampis. Trans. Theory
Stat. Phys. 1, 2, 101 (1971).

.. _To:



**(To)** Q.D. To, V.H. Vu, G. Lauriat, and
C. Leonard. J. Math. Phys. 56, 103101 (2015).


