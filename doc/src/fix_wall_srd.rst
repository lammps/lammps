.. index:: fix wall/srd

fix wall/srd command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID wall/srd face arg ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/srd = style name of this fix command
* one or more face/arg pairs may be appended
* face = *xlo* or *xhi* or *ylo* or *yhi* or *zlo* or *zhi*
  
  .. parsed-literal::
  
       *xlo*\ ,\ *ylo*\ ,\ *zlo* arg = EDGE or constant or variable
         EDGE = current lo edge of simulation box
         constant = number like 0.0 or -30.0 (distance units)
         variable = :doc:`equal-style variable <variable>` like v_x or v_wiggle
       *xhi*\ ,\ *yhi*\ ,\ *zhi* arg = EDGE or constant or variable
         EDGE = current hi edge of simulation box
         constant = number like 50.0 or 100.3 (distance units)
         variable = :doc:`equal-style variable <variable>` like v_x or v_wiggle

* zero or more keyword/value pairs may be appended
* keyword = *units*
  
  .. parsed-literal::
  
       *units* value = *lattice* or *box*
         *lattice* = the wall position is defined in lattice units
         *box* = the wall position is defined in simulation box units



Examples
""""""""


.. parsed-literal::

   fix xwalls all wall/srd xlo EDGE xhi EDGE
   fix walls all wall/srd xlo 0.0 ylo 10.0 units box
   fix top all wall/srd zhi v_pressdown

Description
"""""""""""

Bound the simulation with one or more walls which interact with
stochastic reaction dynamics (SRD) particles as slip (smooth) or
no-slip (rough) flat surfaces.  The wall interaction is actually
invoked via the :doc:`fix srd <fix_srd>` command, only on the group of
SRD particles it defines, so the group setting for the fix wall/srd
command is ignored.

A particle/wall collision occurs if an SRD particle moves outside the
wall on a timestep.  This alters the position and velocity of the SRD
particle and imparts a force to the wall.

The *collision* and *Tsrd* settings specified via the :doc:`fix srd <fix_srd>` command affect the SRD/wall collisions.  A *slip*
setting for the *collision* keyword means that the tangential
component of the SRD particle momentum is preserved.  Thus only a
normal force is imparted to the wall.  The normal component of the new
SRD velocity is sampled from a Gaussian distribution at temperature
*Tsrd*\ .

For a *noslip* setting of the *collision* keyword, both the normal and
tangential components of the new SRD velocity are sampled from a
Gaussian distribution at temperature *Tsrd*\ .  Additionally, a new
tangential direction for the SRD velocity is chosen randomly.  This
collision style imparts both a normal and tangential force to the
wall.

Up to 6 walls or faces can be specified in a single command: *xlo*\ ,
*xhi*\ , *ylo*\ , *yhi*\ , *zlo*\ , *zhi*\ .  A *lo* face reflects particles
that move to a coordinate less than the wall position, back in the
*hi* direction.  A *hi* face reflects particles that move to a
coordinate higher than the wall position, back in the *lo* direction.

The position of each wall can be specified in one of 3 ways: as the
EDGE of the simulation box, as a constant value, or as a variable.  If
EDGE is used, then the corresponding boundary of the current
simulation box is used.  If a numeric constant is specified then the
wall is placed at that position in the appropriate dimension (x, y, or
z).  In both the EDGE and constant cases, the wall will never move.
If the wall position is a variable, it should be specified as v\_name,
where name is an :doc:`equal-style variable <variable>` name.  In this
case the variable is evaluated each timestep and the result becomes
the current position of the reflecting wall.  Equal-style variables
can specify formulas with various mathematical functions, and include
:doc:`thermo_style <thermo_style>` command keywords for the simulation
box parameters and timestep and elapsed time.  Thus it is easy to
specify a time-dependent wall position.

.. note::

   Because the trajectory of the SRD particle is tracked as it
   collides with the wall, you must insure that r = distance of the
   particle from the wall, is always > 0 for SRD particles, or LAMMPS
   will generate an error.  This means you cannot start your simulation
   with SRD particles at the wall position *coord* (r = 0) or with
   particles on the wrong side of the wall (r < 0).

.. note::

   If you have 2 or more walls that come together at an edge or
   corner (e.g. walls in the x and y dimensions), then be sure to set the
   *overlap* keyword to *yes* in the :doc:`fix srd <fix_srd>` command,
   since the walls effectively overlap when SRD particles collide with
   them.  LAMMPS will issue a warning if you do not do this.

.. note::

   The walls of this fix only interact with SRD particles, as
   defined by the :doc:`fix srd <fix_srd>` command.  If you are simulating
   a mixture containing other kinds of particles, then you should
   typically use :doc:`another wall command <fix_wall>` to act on the other
   particles.  Since SRD particles will be colliding both with the walls
   and the other particles, it is important to insure that the other
   particle's finite extent does not overlap an SRD wall.  If you do not
   do this, you may generate errors when SRD particles end up "inside"
   another particle or a wall at the beginning of a collision step.

The *units* keyword determines the meaning of the distance units used
to define a wall position, but only when a numeric constant is used.
It is not relevant when EDGE or a variable is used to specify a face
position.

A *box* value selects standard distance units as defined by the
:doc:`units <units>` command, e.g. Angstroms for units = real or metal.
A *lattice* value means the distance units are in lattice spacings.
The :doc:`lattice <lattice>` command must have been previously used to
define the lattice spacings.


----------


Here are examples of variable definitions that move the wall position
in a time-dependent fashion using equal-style
:doc:`variables <variable>`.


.. parsed-literal::

   variable ramp equal ramp(0,10)
   fix 1 all wall/srd xlo v_ramp

   variable linear equal vdisplace(0,20)
   fix 1 all wall/srd xlo v_linear

   variable wiggle equal swiggle(0.0,5.0,3.0)
   fix 1 all wall/srd xlo v_wiggle

   variable wiggle equal cwiggle(0.0,5.0,3.0)
   fix 1 all wall/srd xlo v_wiggle

The ramp(lo,hi) function adjusts the wall position linearly from lo to
hi over the course of a run.  The displace(c0,velocity) function does
something similar using the equation position = c0 + velocity\*delta,
where delta is the elapsed time.

The swiggle(c0,A,period) function causes the wall position to
oscillate sinusoidally according to this equation, where omega = 2 PI
/ period:


.. parsed-literal::

   position = c0 + A sin(omega\*delta)

The cwiggle(c0,A,period) function causes the wall position to
oscillate sinusoidally according to this equation, which will have an
initial wall velocity of 0.0, and thus may impose a gentler
perturbation on the particles:


.. parsed-literal::

   position = c0 + A (1 - cos(omega\*delta))


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global array of values which can be accessed by
various :doc:`output commands <Howto_output>`.  The number of rows in
the array is equal to the number of walls defined by the fix.  The
number of columns is 3, for the x,y,z components of force on each
wall.

Note that an outward normal force on a wall will be a negative value
for *lo* walls and a positive value for *hi* walls.  The array values
calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


Any dimension (xyz) that has an SRD wall must be non-periodic.

Related commands
""""""""""""""""

:doc:`fix srd <fix_srd>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
