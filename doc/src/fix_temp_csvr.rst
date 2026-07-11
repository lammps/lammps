.. index:: fix temp/csvr
.. index:: fix temp/csld

fix temp/csvr command
=====================

fix temp/csld command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID temp/csvr Tstart Tstop Tdamp seed
   fix ID group-ID temp/csld Tstart Tstop Tdamp seed

* ID, group-ID are documented in :doc:`fix <fix>` command
* temp/csvr or temp/csld = style name of this fix command
* Tstart,Tstop = desired temperature at start/end of run

  .. parsed-literal::

       Tstart can be a variable (see below)

* Tdamp = temperature damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all temp/csvr 300.0 300.0 100.0 54324
   fix 1 all temp/csld 100.0 300.0 10.0 123321

Description
"""""""""""

Adjust the temperature with a canonical sampling thermostat that uses
global velocity rescaling with Hamiltonian dynamics (\ *temp/csvr*\ )
:ref:`(Bussi1) <Bussi1>`, or Langevin dynamics (\ *temp/csld*\ )
:ref:`(Bussi2) <Bussi2>`.  In the case of *temp/csvr* the thermostat is
similar to the empirical Berendsen thermostat in
:doc:`temp/berendsen <fix_temp_berendsen>`, but chooses the actual
scaling factor from a suitably chosen (gaussian) distribution rather
than having it determined from the time constant directly. In the case
of *temp/csld* the velocities are updated to a linear combination of
the current velocities with a gaussian distribution of velocities at
the desired temperature.  Both thermostats are applied every timestep.

The thermostat is applied to only the translational degrees of freedom
for the particles, which is an important consideration for finite-size
particles which have rotational degrees of freedom are being
thermostatted with these fixes.  The translational degrees of freedom
can also have a bias velocity removed from them before thermostatting
takes place; see the description below.

The desired temperature at each timestep is a ramped value during the
run from *Tstart* to *Tstop*\ .  The *Tdamp* parameter is specified in
time units and determines how rapidly the temperature is relaxed.  For
example, a value of 100.0 means to relax the temperature in a timespan
of (roughly) 100 time units (tau or fs or ps - see the
:doc:`units <units>` command).

*Tstart* can be specified as an equal-style :doc:`variable <variable>`.
In this case, the *Tstop* setting is ignored.  If the value is a
variable, it should be specified as v_name, where name is the variable
name.  In this case, the variable will be evaluated each timestep, and
its value used to determine the target temperature.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent temperature.

.. note::

   Unlike the :doc:`fix nvt <fix_nh>` command which performs
   Nose/Hoover thermostatting AND time integration, these fixes do NOT
   perform time integration. They only modify velocities to effect
   thermostatting.  Thus you must use a separate time integration fix,
   like :doc:`fix nve <fix_nve>` to actually update the positions of atoms
   using the modified velocities.  Likewise, these fixes should not
   normally be used on atoms that also have their temperature controlled
   by another fix - e.g. by :doc:`fix nvt <fix_nh>` or :doc:`fix langevin <fix_langevin>` commands.

See the :doc:`Howto thermostat <Howto_thermostat>` page for a
discussion of different ways to compute temperature and perform
thermostatting.

These fixes compute a temperature each timestep.  To do this, the fix
creates its own compute of style "temp", as if this command had been
issued:

.. code-block:: LAMMPS

   compute fix-ID_temp group-ID temp

See the :doc:`compute temp <compute_temp>` command for details.  Note
that the ID of the new compute is the fix-ID + underscore + "temp",
and the group for the new compute is the same as the fix group.

Note that this is NOT the compute used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo_temp*.
This means you can change the attributes of this fix's temperature
(e.g. its degrees-of-freedom) via the
:doc:`compute_modify <compute_modify>` command or print this temperature
during thermodynamic output via the :doc:`thermo_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo_temp* will have no
effect on this fix.

Like other fixes that perform thermostatting, this fix can be used
with :doc:`compute commands <compute>` that remove a "bias" from the
atom velocities.  E.g. to apply the thermostat only to atoms within a
spatial :doc:`region <region>`, or to remove the center-of-mass
velocity from a group of atoms, or to remove the x-component of
velocity from the calculation.

This is not done by default, but only if the :doc:`fix_modify
<fix_modify>` command is used to assign a temperature compute to this
fix that includes such a bias term.  See the doc pages for individual
:doc:`compute temp commands <compute>` to determine which ones include
a bias.  In this case, the thermostat works in the following manner:
bias is removed from each atom, thermostatting is performed on the
remaining thermal degrees of freedom, and the bias is added back in.

An important feature of these thermostats is that they have an
associated effective energy that is a constant of motion.  The
effective energy is the total energy (kinetic + potential) plus the
accumulated kinetic energy changes due to the thermostat. The latter
quantity is the global scalar computed by these fixes. This feature is
useful to check the integration of the equations of motion against
discretization errors. In other words, the conservation of the
effective energy can be used to choose an appropriate integration
:doc:`timestep <timestep>`. This is similar to the usual paradigm of
checking the conservation of the total energy in the microcanonical
ensemble.


----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These fixes write the cumulative global energy change and the random
number generator states to :doc:`binary restart files <restart>`.  See
the :doc:`read_restart <read_restart>` command for info on how to
re-specify a fix in an input script that reads a restart file, so that
the selected fix continues in an uninterrupted fashion.  The random
number generator state can only be restored when the number of
processors remains unchanged from what is recorded in the restart
file.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by these
fixes.  You can use it to assign a temperature :doc:`compute
<compute>` you have defined to these fixes which will be used in its
thermostatting procedure, as described above.  For consistency, the
group used by these fixes and by the compute should be the same.

The cumulative energy change in the system imposed by these fixes is
included in the :doc:`thermodynamic output <thermo_style>` keywords
*ecouple* and *econserve*.  See the :doc:`thermo_style <thermo_style>`
doc page for details.

These fixes compute a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the same
cumulative energy change due to this fix described in the previous
paragraph.  The scalar value calculated by this fix is "extensive".

These fixes can ramp their target temperature over multiple runs,
using the *start* and *stop* keywords of the :doc:`run <run>` command.
See the :doc:`run <run>` command for details of how to do this.

These fixes are not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

Fix *temp/csld* is not compatible with :doc:`fix shake <fix_shake>`.

These fixes are part of the EXTRA-FIX package.  They are only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

These fixes can be used with dynamic groups as defined by the
:doc:`group <group>` command.  Likewise it can be used with groups to
which atoms are added or deleted over time, e.g. a deposition
simulation.  However, the conservation properties of the thermostat
and barostat are defined for systems with a static set of atoms.  You
may observe odd behavior if the atoms in a group vary dramatically
over time or the atom count becomes very small.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`, :doc:`fix nvt <fix_nh>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix langevin <fix_langevin>`,
:doc:`fix_modify <fix_modify>`, :doc:`compute temp <compute_temp>`,
:doc:`fix temp/berendsen <fix_temp_berendsen>`

Default
"""""""

none

----------

.. _Bussi1:

.. _Bussi2:

**(Bussi1)** Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101(2007)

**(Bussi2)** Bussi and Parrinello, Phys. Rev. E 75, 056707 (2007)
