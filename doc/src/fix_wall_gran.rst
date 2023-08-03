.. index:: fix wall/gran
.. index:: fix wall/gran/kk

fix wall/gran command
=====================

Accelerator Variants: *wall/gran/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID wall/gran fstyle fstyle_params wallstyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/gran = style name of this fix command
* fstyle = style of force interactions between particles and wall

  .. parsed-literal::

       possible choices: hooke, hooke/history, hertz/history, granular

* fstyle_params = parameters associated with force interaction style

  .. parsed-literal::

       For *hooke*, *hooke/history*, and *hertz/history*, *fstyle_params* are:
             Kn = elastic constant for normal particle repulsion (force/distance units or pressure units - see discussion below)
             Kt = elastic constant for tangential contact (force/distance units or pressure units - see discussion below)
             gamma_n = damping coefficient for collisions in normal direction (1/time units or 1/time-distance units - see discussion below)
             gamma_t = damping coefficient for collisions in tangential direction (1/time units or 1/time-distance units - see discussion below)
             xmu = static yield criterion (unitless value between 0.0 and 1.0e4)
             dampflag = 0 or 1 if tangential damping force is excluded or included
             optional keyword = *limit_damping*, limit damping to prevent attractive interaction

  .. parsed-literal::

       For *granular*, *fstyle_params* are set using the same syntax as for the *pair_coeff* command of :doc:`pair_style granular <pair_granular>`

* wallstyle = *xplane* or *yplane* or *zplane* or *zcylinder*
* args = list of arguments for a particular style

  .. parsed-literal::

       *xplane* or *yplane* or *zplane* args = lo hi
         lo,hi = position of lower and upper plane (distance units), either can be NULL)
       *zcylinder* args = radius
         radius = cylinder radius (distance units)

* zero or more keyword/value pairs may be appended to args
* keyword = *wiggle* or *shear* or *contacts* or *temperature*

  .. parsed-literal::

       *wiggle* values = dim amplitude period
         dim = *x* or *y* or *z*
         amplitude = size of oscillation (distance units)
         period = time of oscillation (time units)
       *shear* values = dim vshear
         dim = *x* or *y* or *z*
         vshear = magnitude of shear velocity (velocity units)
      *contacts* value = none
         generate contact information for each particle
      *temperature* value = temperature
         specify temperature of wall


Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all wall/gran hooke  200000.0 NULL 50.0 NULL 0.5 0 xplane -10.0 10.0
   fix 1 all wall/gran hooke/history 200000.0 NULL 50.0 NULL 0.5 0 zplane 0.0 NULL
   fix 2 all wall/gran hooke 100000.0 20000.0 50.0 30.0 0.5 1 zcylinder 15.0 wiggle z 3.0 2.0
   fix 3 all wall/gran granular hooke 1000.0 50.0 tangential linear_nohistory 1.0 0.4 damping velocity region myBox
   fix 4 all wall/gran granular jkr 1e5 1500.0 0.3 10.0 tangential mindlin NULL 1.0 0.5 rolling sds 500.0 200.0 0.5 twisting marshall region myCone
   fix 5 all wall/gran granular dmt 1e5 0.2 0.3 10.0 tangential mindlin NULL 1.0 0.5 rolling sds 500.0 200.0 0.5 twisting marshall damping tsuji heat 10 region myCone temperature 1.0
   fix 6 all wall/gran hooke  200000.0 NULL 50.0 NULL 0.5 0 xplane -10.0 10.0 contacts

Description
"""""""""""

Bound the simulation domain of a granular system with a frictional
wall.  All particles in the group interact with the wall when they are
close enough to touch it.

The nature of the wall/particle interactions are determined by the
*fstyle* setting.  It can be any of the styles defined by the
:doc:`pair_style gran/\* <pair_gran>` or the more general
:doc:`pair_style granular <pair_granular>` commands.  Currently the
options are *hooke*, *hooke/history*, or *hertz/history* for the
former, and *granular* with all the possible options of the associated
*pair_coeff* command for the latter.  The equation for the force
between the wall and particles touching it is the same as the
corresponding equation on the :doc:`pair_style gran/\* <pair_gran>` and
:doc:`pair_style granular <pair_granular>` doc pages, in the limit of
one of the two particles going to infinite radius and mass (flat wall).
Specifically, delta = radius - r = overlap of particle with wall, m_eff
= mass of particle, and the effective radius of contact = RiRj/Ri+Rj is
set to the radius of the particle.

The parameters *Kn*, *Kt*, *gamma_n*, *gamma_t*, *xmu*, *dampflag*,
and the optional keyword *limit_damping*
have the same meaning and units as those specified with the
:doc:`pair_style gran/\* <pair_gran>` commands.  This means a NULL can be
used for either *Kt* or *gamma_t* as described on that page.  If a
NULL is used for *Kt*, then a default value is used where *Kt* = 2/7
*Kn*\ .  If a NULL is used for *gamma_t*, then a default value is used
where *gamma_t* = 1/2 *gamma_n*.

All the model choices for cohesion, tangential friction, rolling
friction and twisting friction supported by the :doc:`pair_style granular <pair_granular>` through its *pair_coeff* command are also
supported for walls. These are discussed in greater detail on the doc
page for :doc:`pair_style granular <pair_granular>`.

Note that you can choose a different force styles and/or different
values for the wall/particle coefficients than for particle/particle
interactions.  E.g. if you wish to model the wall as a different
material.

.. note::

   As discussed on the page for :doc:`pair_style gran/\* <pair_gran>`,
   versions of LAMMPS before 9Jan09 used a different equation for
   Hertzian interactions.  This means Hertizian wall/particle
   interactions have also changed.  They now include a sqrt(radius) term
   which was not present before.  Also the previous versions used Kn and
   Kt from the pairwise interaction and hardwired dampflag to 1, rather
   than letting them be specified directly.  This means you can set the
   values of the wall/particle coefficients appropriately in the current
   code to reproduce the results of a previous Hertzian monodisperse
   calculation.  For example, for the common case of a monodisperse
   system with particles of diameter 1, Kn, Kt, gamma_n, and gamma_s
   should be set sqrt(2.0) larger than they were previously.

The effective mass *m_eff* in the formulas listed on the :doc:`pair_style granular <pair_gran>` page is the mass of the particle for
particle/wall interactions (mass of wall is infinite).  If the
particle is part of a rigid body, its mass is replaced by the mass of
the rigid body in those formulas.  This is determined by searching for
a :doc:`fix rigid <fix_rigid>` command (or its variants).

The *wallstyle* can be planar or cylindrical.  The 3 planar options
specify a pair of walls in a dimension.  Wall positions are given by
*lo* and *hi*\ .  Either of the values can be specified as NULL if a
single wall is desired.  For a *zcylinder* wallstyle, the cylinder's
axis is at x = y = 0.0, and the radius of the cylinder is specified.

Optionally, the wall can be moving, if the *wiggle* or *shear*
keywords are appended.  Both keywords cannot be used together.

For the *wiggle* keyword, the wall oscillates sinusoidally, similar to
the oscillations of particles which can be specified by the :doc:`fix move <fix_move>` command.  This is useful in packing simulations of
granular particles.  The arguments to the *wiggle* keyword specify a
dimension for the motion, as well as it's *amplitude* and *period*\ .
Note that if the dimension is in the plane of the wall, this is
effectively a shearing motion.  If the dimension is perpendicular to
the wall, it is more of a shaking motion.  A *zcylinder* wall can only
be wiggled in the z dimension.

Each timestep, the position of a wiggled wall in the appropriate *dim*
is set according to this equation:

.. parsed-literal::

   position = coord + A - A cos (omega \* delta)

where *coord* is the specified initial position of the wall, *A* is
the *amplitude*, *omega* is 2 PI / *period*, and *delta* is the time
elapsed since the fix was specified.  The velocity of the wall is set
to the derivative of this expression.

For the *shear* keyword, the wall moves continuously in the specified
dimension with velocity *vshear*\ .  The dimension must be tangential to
walls with a planar *wallstyle*, e.g. in the *y* or *z* directions for
an *xplane* wall.  For *zcylinder* walls, a dimension of *z* means the
cylinder is moving in the z-direction along it's axis.  A dimension of
*x* or *y* means the cylinder is spinning around the z-axis, either in
the clockwise direction for *vshear* > 0 or counter-clockwise for
*vshear* < 0.  In this case, *vshear* is the tangential velocity of
the wall at whatever *radius* has been defined.

The *temperature* keyword is used to assign a temperature to the wall.
The following value can either be a numeric value or an equal-style
:doc:`variable <variable>`.  If the value is a variable, it should be
specified as v_name, where name is the variable name.  In this case, the
variable will be evaluated each timestep, and its value used to determine
the temperature. This option must be used in conjunction with a heat
conduction model defined in :doc:`pair_style granular <pair_granular>`,
:doc:`fix property/atom <fix_property_atom>` to store temperature and a
heat flow, and :doc:`fix heat/flow <fix_heat_flow>` to integrate heat
flow.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the shear friction state of atoms interacting with the
wall to :doc:`binary restart files <restart>`, so that a simulation can
continue correctly if granular potentials with shear "history" effects
are being used.  See the :doc:`read_restart <read_restart>` command for
info on how to re-specify a fix in an input script that reads a
restart file, so that the operation of the fix continues in an
uninterrupted fashion.

If the :code:`contacts` option is used, this fix generates a per-atom array
with 8 columns as output, containing the contact information for owned
particles (nlocal on each processor). All columns in this per-atom array will
be zero if no contact has occurred.  The values of these columns are listed in
the following table:

+-------+----------------------------------------------------+----------------+
| Index | Value                                              | Units          |
+=======+====================================================+================+
|     1 | 1.0 if particle is in contact with wall,           |                |
|       | 0.0 otherwise                                      |                |
+-------+----------------------------------------------------+----------------+
|     2 | Force :math:`f_x` exerted by the wall              | force units    |
+-------+----------------------------------------------------+----------------+
|     3 | Force :math:`f_y` exerted by the wall              | force units    |
+-------+----------------------------------------------------+----------------+
|     4 | Force :math:`f_z` exerted by the wall              | force units    |
+-------+----------------------------------------------------+----------------+
|     5 | :math:`x`-coordinate of contact point on wall      | distance units |
+-------+----------------------------------------------------+----------------+
|     6 | :math:`y`-coordinate of contact point on wall      | distance units |
+-------+----------------------------------------------------+----------------+
|     7 | :math:`z`-coordinate of contact point on wall      | distance units |
+-------+----------------------------------------------------+----------------+
|     8 | Radius :math:`r` of atom                           | distance units |
+-------+----------------------------------------------------+----------------+

None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.
No parameter of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command. This fix is not invoked during :doc:`energy
minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Any dimension (xyz) that has a granular wall must be non-periodic.

Related commands
""""""""""""""""

:doc:`fix move <fix_move>`,
:doc:`fix wall/gran/region <fix_wall_gran_region>`,
:doc:`pair_style gran/\* <pair_gran>`
:doc:`pair_style granular <pair_granular>`

Default
"""""""

none
