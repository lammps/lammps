.. index:: fix wall/piston

fix wall/piston command
=======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID wall/piston face ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/piston = style name of this fix command
* face = *zlo*
* zero or more keyword/value pairs may be appended
* keyword = *pos* or *vel* or *ramp* or *units*
  
  .. parsed-literal::
  
       *pos* args = z
         z = z coordinate at which the piston begins (distance units)
       *vel* args = vz
         vz = final velocity of the piston (velocity units)
       *ramp* = use a linear velocity ramp from 0 to vz
       *temp* args = target damp seed extent
         target = target velocity for region immediately ahead of the piston
         damp = damping parameter (time units)
         seed = random number seed for langevin kicks
         extent = extent of thermostatted region (distance units)
       *units* value = *lattice* or *box*
         *lattice* = the wall position is defined in lattice units
         *box* = the wall position is defined in simulation box units



Examples
""""""""


.. parsed-literal::

   fix xwalls all wall/piston zlo
   fix walls all wall/piston zlo pos 1.0 vel 10.0 units box
   fix top all wall/piston zlo vel 10.0 ramp

Description
"""""""""""

Bound the simulation with a moving wall which reflect particles in the
specified group and drive the system with an effective infinite-mass
piston capable of driving shock waves.

A momentum mirror technique is used, which means that if an atom (or
the wall) moves such that an atom is outside the wall on a timestep by
a distance delta (e.g. due to :doc:`fix nve <fix_nve>`), then it is put
back inside the face by the same delta, and the velocity relative to
the moving wall is flipped in z.  For instance, a stationary particle
hit with a piston wall with velocity vz, will end the timestep with a
velocity of 2\*vz.

Currently the *face* keyword can only be *zlo*\ .  This creates a piston
moving in the positive z direction.  Particles with z coordinate less
than the wall position are reflected to a z coordinate greater than
the wall position.  If the piston velocity is vpz and the particle
velocity before reflection is vzi, the particle velocity after
reflection is -vzi + 2\*vpz.

The initial position of the wall can be specified by the *pos* keyword.

The final velocity of the wall can be specified by the *vel* keyword

The *ramp* keyword will cause the wall/piston to adjust the velocity
linearly from zero velocity to *vel* over the course of the run. If
the *ramp* keyword is omitted then the wall/piston moves at a constant
velocity defined by *vel*\ .

The *temp* keyword will cause the region immediately in front of the
wall/piston to be thermostatted with a Langevin thermostat.  This
region moves with the piston.  The damping and kicking are measured in
the reference frame of the piston.  So, a temperature of zero would
mean all particles were moving at exactly the speed of the
wall/piston.

The *units* keyword determines the meaning of the distance units used
to define a wall position, but only when a numeric constant is used.

A *box* value selects standard distance units as defined by the
:doc:`units <units>` command, e.g. Angstroms for units = real or metal.
A *lattice* value means the distance units are in lattice spacings.
The :doc:`lattice <lattice>` command must have been previously used to
define the lattice spacings.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix style is part of the SHOCK package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

The face that has the wall/piston must be boundary type 's'
(shrink-wrapped). The opposing face can be
any boundary type other than periodic.

A wall/piston should not be used with rigid bodies such as those
defined by a "fix rigid" command.  This is because the wall/piston
displaces atoms directly rather than exerting a force on them.

Related commands
""""""""""""""""

:doc:`fix wall/reflect <fix_wall>` command, :doc:`fix append/atoms <fix_append_atoms>` command

Default
"""""""

The keyword defaults are pos = 0, vel = 0, units = lattice.
