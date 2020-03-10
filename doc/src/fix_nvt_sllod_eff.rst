.. index:: fix nvt/sllod/eff

fix nvt/sllod/eff command
=========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nvt/sllod/eff keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nvt/sllod/eff = style name of this fix command
* additional thermostat related keyword/value pairs from the :doc:`fix nvt/eff <fix_nh_eff>` command can be appended

Examples
""""""""

.. parsed-literal::

   fix 1 all nvt/sllod/eff temp 300.0 300.0 0.1
   fix 1 all nvt/sllod/eff temp 300.0 300.0 0.1 drag 0.2

Description
"""""""""""

Perform constant NVT integration to update positions and velocities
each timestep for nuclei and electrons in the group for the :doc:`electron force field <pair_eff>` model, using a Nose/Hoover temperature
thermostat.  V is volume; T is temperature.  This creates a system
trajectory consistent with the canonical ensemble.

The operation of this fix is exactly like that described by the :doc:`fix nvt/sllod <fix_nvt_sllod>` command, except that the radius and
radial velocity of electrons are also updated and thermostatted.
Likewise the temperature calculated by the fix, using the compute it
creates (as discussed in the :doc:`fix nvt, npt, and nph <fix_nh>` doc
page), is performed with a :doc:`compute temp/deform/eff <compute_temp_deform_eff>` command that includes
the eFF contribution to the temperature from the electron radial
velocity.

**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix writes the state of the Nose/Hoover thermostat to :doc:`binary restart files <restart>`.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a :doc:`compute <compute>` you have
defined to this fix which will be used in its thermostatting
procedure.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Nose/Hoover thermostatting to
the system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes the same global scalar and global vector of
quantities as does the :doc:`fix nvt/eff <fix_nh_eff>` command.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the USER-EFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix works best without Nose-Hoover chain thermostats, i.e. using
tchain = 1.  Setting tchain to larger values can result in poor
equilibration.

Related commands
""""""""""""""""

:doc:`fix nve/eff <fix_nve_eff>`, :doc:`fix nvt/eff <fix_nh_eff>`, :doc:`fix langevin/eff <fix_langevin_eff>`, :doc:`fix nvt/sllod <fix_nvt_sllod>`, :doc:`fix_modify <fix_modify>`, :doc:`compute temp/deform/eff <compute_temp_deform_eff>`

Default
"""""""

Same as :doc:`fix nvt/eff <fix_nh_eff>`, except tchain = 1.

----------

.. _Tuckerman2:

**(Tuckerman)** Tuckerman, Mundy, Balasubramanian, Klein, J Chem Phys,
106, 5615 (1997).
