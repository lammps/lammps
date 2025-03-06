.. index:: fix gravity
.. index:: fix gravity/omp
.. index:: fix gravity/kk

fix gravity command
===================

Accelerator Variants: *gravity/omp*, *gravity/kk*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group gravity magnitude style args

* ID, group are documented in :doc:`fix <fix>` command
* gravity = style name of this fix command
* magnitude = size of acceleration (force/mass units)
* magnitude can be a variable (see below)
* style = *chute* or *spherical* or *gradient* or *vector*

  .. parsed-literal::

       *chute* args = angle
         angle = angle in +x away from -z or -y axis in 3d/2d (in degrees)
         angle can be a variable (see below)
       *spherical* args = phi theta
         phi = azimuthal angle from +x axis (in degrees)
         theta = angle from +z or +y axis in 3d/2d (in degrees)
         phi or theta can be a variable (see below)
       *vector* args = x y z
         x y z = vector direction to apply the acceleration
         x or y or z can be a variable (see below)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all gravity 1.0 chute 24.0
   fix 1 all gravity v_increase chute 24.0
   fix 1 all gravity 1.0 spherical 0.0 -180.0
   fix 1 all gravity 10.0 spherical v_phi v_theta
   fix 1 all gravity 100.0 vector 1 1 0

Description
"""""""""""

Impose an additional acceleration on each particle in the group.  This
fix is typically used with granular systems to include a "gravity"
term acting on the macroscopic particles.  More generally, it can
represent any kind of driving field, e.g. a pressure gradient inducing
a Poiseuille flow in a fluid.  Note that this fix operates differently
than the :doc:`fix addforce <fix_addforce>` command.  The addforce fix
adds the same force to each atom, independent of its mass.  This
command imparts the same acceleration to each atom (force/mass).

The *magnitude* of the acceleration is specified in force/mass units.
For granular systems (LJ units) this is typically 1.0.  See the
:doc:`units <units>` command for details.

Style *chute* is typically used for simulations of chute flow where
the specified *angle* is the chute angle, with flow occurring in the +x
direction.  For 3d systems, the tilt is away from the z axis; for 2d
systems, the tilt is away from the y axis.

Style *spherical* allows an arbitrary 3d direction to be specified for
the acceleration vector.  *Phi* and *theta* are defined in the usual
spherical coordinates.  Thus for acceleration acting in the -z
direction, *theta* would be 180.0 (or -180.0).  *Theta* = 90.0 and
*phi* = -90.0 would mean acceleration acts in the -y direction.  For
2d systems, *phi* is ignored and *theta* is an angle in the xy plane
where *theta* = 0.0 is the y-axis.

Style *vector* imposes an acceleration in the vector direction given
by (x,y,z).  Only the direction of the vector is important; it's
length is ignored.  For 2d systems, the *z* component is ignored.

Any of the quantities *magnitude*, *angle*, *phi*, *theta*, *x*, *y*,
*z* which define the gravitational magnitude and direction, can be
specified as an equal-style :doc:`variable <variable>`.  If the value is
a variable, it should be specified as v_name, where name is the
variable name.  In this case, the variable will be evaluated each
timestep, and its value used to determine the quantity.  You should
ensure that the variable calculates a result in the appropriate units,
e.g. force/mass or degrees.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent gravitational
field.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the gravitational potential energy of the system to
the global potential energy of the system as part of
:doc:`thermodynamic output <thermo_style>`.  The default setting for
this fix is :doc:`fix_modify energy no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its forces. Default is the
outermost level.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  This scalar is the
gravitational potential energy of the particles in the defined field,
namely mass \* (g dot x) for each particles, where x and mass are the
particles position and mass, and g is the gravitational field.  The
scalar value calculated by this fix is "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`atom_style sphere <atom_style>`, :doc:`fix addforce <fix_addforce>`

Default
"""""""

none
