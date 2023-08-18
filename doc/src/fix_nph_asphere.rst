.. index:: fix nph/asphere
.. index:: fix nph/asphere/omp

fix nph/asphere command
=======================

Accelerator Variants: *nph/asphere/omp*

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nph/asphere args keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nph/asphere = style name of this fix command
* additional barostat related keyword/value pairs from the :doc:`fix nph <fix_nh>` command can be appended

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nph/asphere iso 0.0 0.0 1000.0
   fix 2 all nph/asphere x 5.0 5.0 1000.0
   fix 2 all nph/asphere x 5.0 5.0 1000.0 drag 0.2
   fix 2 water nph/asphere aniso 0.0 0.0 1000.0 dilate partial

Description
"""""""""""

Perform constant NPH integration to update position, velocity,
orientation, and angular velocity each timestep for aspherical or
ellipsoidal particles in the group using a Nose/Hoover pressure
barostat.  P is pressure; H is enthalpy.  This creates a system
trajectory consistent with the isenthalpic ensemble.

This fix differs from the :doc:`fix nph <fix_nh>` command, which assumes
point particles and only updates their position and velocity.

Additional parameters affecting the barostat are specified by keywords
and values documented with the :doc:`fix nph <fix_nh>` command.  See,
for example, discussion of the *aniso*, and *dilate* keywords.

The particles in the fix group are the only ones whose velocities and
positions are updated by the velocity/position update portion of the
NPH integration.

Regardless of what particles are in the fix group, a global pressure is
computed for all particles.  Similarly, when the size of the simulation
box is changed, all particles are re-scaled to new positions, unless the
keyword *dilate* is specified with a value of *partial*, in which case
only the particles in the fix group are re-scaled.  The latter can be
useful for leaving the coordinates of particles in a solid substrate
unchanged and controlling the pressure of a surrounding fluid.

----------

This fix computes a temperature and pressure each timestep.  To do
this, the fix creates its own computes of style "temp/asphere" and
"pressure", as if these commands had been issued:

.. code-block:: LAMMPS

   compute fix-ID_temp all temp/asphere
   compute fix-ID_press all pressure fix-ID_temp

See the :doc:`compute temp/asphere <compute_temp_asphere>` and :doc:`compute pressure <compute_pressure>` commands for details.  Note that the
IDs of the new computes are the fix-ID + underscore + "temp" or fix_ID
+ underscore + "press", and the group for the new computes is "all"
since pressure is computed for the entire system.

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo_style <thermo_style>` command) with ID = *thermo_temp*
and *thermo_press*.  This means you can change the attributes of this
fix's temperature or pressure via the
:doc:`compute_modify <compute_modify>` command or print this temperature
or pressure during thermodynamic output via the :doc:`thermo_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo_temp* or
*thermo_press* will have no effect on this fix.

----------

.. include:: accel_styles.rst

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the Nose/Hoover barostat to :doc:`binary
restart files <restart>`.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* and *press* options are
supported by this fix.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its thermostatting or barostatting procedure.  If you do this, note
that the kinetic energy derived from the compute temperature should be
consistent with the virial term computed using all atoms for the
pressure.  LAMMPS will warn you if you choose to compute temperature
on a subset of atoms.

The cumulative energy change in the system imposed by this fix is
included in the :doc:`thermodynamic output <thermo_style>` keywords
*ecouple* and *econserve*.  See the :doc:`thermo_style <thermo_style>`
doc page for details.

This fix computes the same global scalar and global vector of
quantities as does the :doc:`fix nph <fix_nh>` command.

This fix can ramp its target pressure over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the ASPHERE package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the :doc:`atom_style ellipsoid <atom_style>`
command.

All particles in the group must be finite-size.  They cannot be point
particles, but they can be aspherical or spherical as defined by their
shape attribute.

Related commands
""""""""""""""""

:doc:`fix nph <fix_nh>`, :doc:`fix nve_asphere <fix_nve_asphere>`, :doc:`fix nvt_asphere <fix_nvt_asphere>`, :doc:`fix npt_asphere <fix_npt_asphere>`, :doc:`fix_modify <fix_modify>`

Default
"""""""

none
