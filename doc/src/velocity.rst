.. index:: velocity

velocity command
================

Syntax
""""""

.. parsed-literal::

   velocity group-ID style args keyword value ...

* group-ID = ID of group of atoms whose velocity will be changed
* style = *create* or *set* or *scale* or *ramp* or *zero*

  .. parsed-literal::

       *create* args = temp seed
         temp = temperature value (temperature units)
         seed = random # seed (positive integer)
       *set* args = vx vy vz
         vx,vy,vz = velocity value or NULL (velocity units)
         any of vx,vy,vz van be a variable (see below)
       *scale* arg = temp
         temp = temperature value (temperature units)
       *ramp* args = vdim vlo vhi dim clo chi
         vdim = *vx* or *vy* or *vz*
         vlo,vhi = lower and upper velocity value (velocity units)
         dim = *x* or *y* or *z*
         clo,chi = lower and upper coordinate bound (distance units)
       *zero* arg = *linear* or *angular*
         *linear* = zero the linear momentum
         *angular* = zero the angular momentum

* zero or more keyword/value pairs may be appended
* keyword = *dist* or *sum* or *mom* or *rot* or *temp* or *bias* or *loop* or *units*

  .. parsed-literal::

       *dist* value = *uniform* or *gaussian*
       *sum* value = *no* or *yes*
       *mom* value = *no* or *yes*
       *rot* value = *no* or *yes*
       *temp* value = temperature compute ID
       *bias* value = *no* or *yes*
       *loop* value = *all* or *local* or *geom*
       *rigid* value = fix-ID
         fix-ID = ID of rigid body fix
       *units* value = *box* or *lattice*

Examples
""""""""

.. code-block:: LAMMPS

   velocity all create 300.0 4928459 rot yes dist gaussian
   velocity border set NULL 4.0 v_vz sum yes units box
   velocity flow scale 300.0
   velocity flow ramp vx 0.0 5.0 y 5 25 temp mytemp
   velocity all zero linear

Description
"""""""""""

Set or change the velocities of a group of atoms in one of several
styles.  For each style, there are required arguments and optional
keyword/value parameters.  Not all options are used by each style.
Each option has a default as listed below.

The *create* style generates an ensemble of velocities using a random
number generator with the specified seed at the specified temperature.

The *set* style sets the velocities of all atoms in the group to the
specified values.  If any component is specified as NULL, then it is
not set.  Any of the vx,vy,vz velocity components can be specified as
an equal-style or atom-style :doc:`variable <variable>`.  If the value
is a variable, it should be specified as v_name, where name is the
variable name.  In this case, the variable will be evaluated, and its
value used to determine the velocity component.  Note that if a
variable is used, the velocity it calculates must be in box units, not
lattice units; see the discussion of the *units* keyword below.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters or other parameters.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent
velocity field.

The *scale* style computes the current temperature of the group of
atoms and then rescales the velocities to the specified temperature.

The *ramp* style is similar to that used by the :doc:`compute temp/ramp <compute_temp_ramp>` command.  Velocities ramped
uniformly from vlo to vhi are applied to dimension vx, or vy, or vz.
The value assigned to a particular atom depends on its relative
coordinate value (in dim) from clo to chi.  For the example above, an
atom with y-coordinate of 10 (1/4 of the way from 5 to 25), would be
assigned a x-velocity of 1.25 (1/4 of the way from 0.0 to 5.0).  Atoms
outside the coordinate bounds (less than 5 or greater than 25 in this
case), are assigned velocities equal to vlo or vhi (0.0 or 5.0 in this
case).

The *zero* style adjusts the velocities of the group of atoms so that
the aggregate linear or angular momentum is zero.  No other changes
are made to the velocities of the atoms.  If the *rigid* option is
specified (see below), then the zeroing is performed on individual
rigid bodies, as defined by the :doc:`fix rigid or fix rigid/small <fix_rigid>` commands.  In other words, zero linear
will set the linear momentum of each rigid body to zero, and zero
angular will set the angular momentum of each rigid body to zero.
This is done by adjusting the velocities of the atoms in each rigid
body.

All temperatures specified in the velocity command are in temperature
units; see the :doc:`units <units>` command.  The units of velocities and
coordinates depend on whether the *units* keyword is set to *box* or
*lattice*\ , as discussed below.

For all styles, no atoms are assigned z-component velocities if the
simulation is 2d; see the :doc:`dimension <dimension>` command.

----------

The keyword/value options are used in the following ways by the
various styles.

The *dist* keyword is used by *create*\ .  The ensemble of generated
velocities can be a *uniform* distribution from some minimum to
maximum value, scaled to produce the requested temperature.  Or it can
be a *gaussian* distribution with a mean of 0.0 and a sigma scaled to
produce the requested temperature.

The *sum* keyword is used by all styles, except *zero*\ .  The new
velocities will be added to the existing ones if sum = yes, or will
replace them if sum = no.

The *mom* and *rot* keywords are used by *create*\ .  If mom = yes, the
linear momentum of the newly created ensemble of velocities is zeroed;
if rot = yes, the angular momentum is zeroed.

----------

If specified, the *temp* keyword is used by *create* and *scale* to
specify a :doc:`compute <compute>` that calculates temperature in a
desired way, e.g. by first subtracting out a velocity bias, as
discussed on the :doc:`Howto thermostat <Howto_thermostat>` doc page.
If this keyword is not specified, *create* and *scale* calculate
temperature using a compute that is defined internally as follows:

.. code-block:: LAMMPS

   compute velocity_temp group-ID temp

where group-ID is the same ID used in the velocity command. i.e. the
group of atoms whose velocity is being altered.  This compute is
deleted when the velocity command is finished.  See the :doc:`compute temp <compute_temp>` command for details.  If the calculated
temperature should have degrees-of-freedom removed due to fix
constraints (e.g. SHAKE or rigid-body constraints), then the
appropriate fix command must be specified before the velocity command
is issued.

The *bias* keyword with a *yes* setting is used by *create* and
*scale*\ , but only if the *temp* keyword is also used to specify a
:doc:`compute <compute>` that calculates temperature in a desired way.
If the temperature compute also calculates a velocity bias, the
bias is subtracted from atom velocities before the *create* and
*scale* operations are performed.  After the operations, the bias is
added back to the atom velocities.  See the :doc:`Howto thermostat <Howto_thermostat>` doc page for more discussion of
temperature computes with biases.  Note that the velocity bias is only
applied to atoms in the temperature compute specified with the *temp*
keyword.

As an example, assume atoms are currently streaming in a flow
direction (which could be separately initialized with the *ramp*
style), and you wish to initialize their thermal velocity to a desired
temperature.  In this context thermal velocity means the per-particle
velocity that remains when the streaming velocity is subtracted.  This
can be done using the *create* style with the *temp* keyword
specifying the ID of a :doc:`compute temp/ramp <compute_temp_ramp>` or
:doc:`compute temp/profile <compute_temp_profile>` command, and the
*bias* keyword set to a *yes* value.

----------

The *loop* keyword is used by *create* in the following ways.

If loop = all, then each processor loops over all atoms in the
simulation to create velocities, but only stores velocities for atoms
it owns.  This can be a slow loop for a large simulation.  If atoms
were read from a data file, the velocity assigned to a particular atom
will be the same, independent of how many processors are being used.
This will not be the case if atoms were created using the
:doc:`create_atoms <create_atoms>` command, since atom IDs will likely
be assigned to atoms differently.

If loop = local, then each processor loops over only its atoms to
produce velocities.  The random number seed is adjusted to give a
different set of velocities on each processor.  This is a fast loop,
but the velocity assigned to a particular atom will depend on which
processor owns it.  Thus the results will always be different when a
simulation is run on a different number of processors.

If loop = geom, then each processor loops over only its atoms.  For
each atom a unique random number seed is created, based on the atom's
xyz coordinates.  A velocity is generated using that seed.  This is a
fast loop and the velocity assigned to a particular atom will be the
same, independent of how many processors are used.  However, the set
of generated velocities may be more correlated than if the *all* or
*local* keywords are used.

Note that the *loop geom* keyword will not necessarily assign
identical velocities for two simulations run on different machines.
This is because the computations based on xyz coordinates are
sensitive to tiny differences in the double-precision value for a
coordinate as stored on a particular machine.

----------

The *rigid* keyword only has meaning when used with the *zero* style.
It allows specification of a fix-ID for one of the :doc:`rigid-body fix <fix_rigid>` variants which defines a set of rigid bodies.  The
zeroing of linear or angular momentum is then performed for each rigid
body defined by the fix, as described above.

The *units* keyword is used by *set* and *ramp*\ .  If units = box,
the velocities and coordinates specified in the velocity command are
in the standard units described by the :doc:`units <units>` command
(e.g. Angstroms/fmsec for real units).  If units = lattice, velocities
are in units of lattice spacings per time (e.g. spacings/fmsec) and
coordinates are in lattice spacings.  The :doc:`lattice <lattice>`
command must have been previously used to define the lattice spacing.

----------

Restrictions
""""""""""""

Assigning a temperature via the *create* style to a system with :doc:`rigid bodies <fix_rigid>` or :doc:`SHAKE constraints <fix_shake>` may not
have the desired outcome for two reasons.  First, the velocity command
can be invoked before all of the relevant fixes are created and
initialized and the number of adjusted degrees of freedom (DOFs) is
known.  Thus it is not possible to compute the target temperature
correctly.  Second, the assigned velocities may be partially canceled
when constraints are first enforced, leading to a different
temperature than desired.  A workaround for this is to perform a :doc:`run 0 <run>` command, which insures all DOFs are accounted for
properly, and then rescale the temperature to the desired value before
performing a simulation.  For example:

.. code-block:: LAMMPS

   velocity all create 300.0 12345
   run 0                             # temperature may not be 300K
   velocity all scale 300.0          # now it should be

Related commands
""""""""""""""""

:doc:`fix rigid <fix_rigid>`, :doc:`fix shake <fix_shake>`,
:doc:`lattice <lattice>`

Default
"""""""

The keyword defaults are dist = uniform, sum = no, mom = yes, rot =
no, bias = no, loop = all, and units = lattice.  The temp and rigid
keywords are not defined by default.
