.. index:: fix nvt/body

fix nvt/body command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nvt/body keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nvt/body = style name of this fix command
* additional thermostat related keyword/value pairs from the :doc:`fix nvt <fix_nh>` command can be appended

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt/body temp 300.0 300.0 100.0
   fix 1 all nvt/body temp 300.0 300.0 100.0 drag 0.2

Description
"""""""""""

Perform constant NVT integration to update position, velocity,
orientation, and angular velocity each timestep for body
particles in the group using a Nose/Hoover temperature
thermostat.  V is volume; T is temperature.  This creates a system
trajectory consistent with the canonical ensemble.

This fix differs from the :doc:`fix nvt <fix_nh>` command, which
assumes point particles and only updates their position and velocity.

The thermostat is applied to both the translational and rotational
degrees of freedom for the body particles, assuming a compute is
used which calculates a temperature that includes the rotational
degrees of freedom (see below).  The translational degrees of freedom
can also have a bias velocity removed from them before thermostatting
takes place; see the description below.

Additional parameters affecting the thermostat are specified by
keywords and values documented with the :doc:`fix nvt <fix_nh>`
command.  See, for example, discussion of the *temp* and *drag*
keywords.

This fix computes a temperature each timestep.  To do this, the fix
creates its own compute of style "temp/body", as if this command
had been issued:

.. code-block:: LAMMPS

   compute fix-ID_temp group-ID temp/body

See the :doc:`compute temp/body <compute_temp_body>` command for
details.  Note that the ID of the new compute is the fix-ID +
underscore + "temp", and the group for the new compute is the same as
the fix group.

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

----------

.. include:: accel_styles.rst

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the Nose/Hoover thermostat to
:doc:`binary restart files <restart>`.  See the :doc:`read_restart
<read_restart>` command for info on how to re-specify a fix in an
input script that reads a restart file, so that the operation of the
fix continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a :doc:`compute <compute>` you have
defined to this fix which will be used in its thermostatting
procedure.

The cumulative energy change in the system imposed by this fix is
included in the :doc:`thermodynamic output <thermo_style>` keywords
*ecouple* and *econserve*.  See the :doc:`thermo_style <thermo_style>`
doc page for details.

This fix computes the same global scalar and global vector of
quantities as does the :doc:`fix nvt <fix_nh>` command.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BODY package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the :doc:`atom_style body <atom_style>`
command.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix nve_body <fix_nve_body>`, :doc:`fix npt_body <fix_npt_body>`, :doc:`fix_modify <fix_modify>`

Default
"""""""

none
