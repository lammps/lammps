.. index:: fix temp/rescale/eff

fix temp/rescale/eff command
============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID temp/rescale/eff N Tstart Tstop window fraction

* ID, group-ID are documented in :doc:`fix <fix>` command
* temp/rescale/eff = style name of this fix command
* N = perform rescaling every N steps
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* window = only rescale if temperature is outside this window (temperature units)
* fraction = rescale to target temperature by this fraction

Examples
""""""""


.. parsed-literal::

   fix 3 flow temp/rescale/eff 10 1.0 100.0 0.02 1.0

Description
"""""""""""

Reset the temperature of a group of nuclei and electrons in the
:doc:`electron force field <pair_eff>` model by explicitly rescaling
their velocities.

The operation of this fix is exactly like that described by the :doc:`fix temp/rescale <fix_temp_rescale>` command, except that the rescaling
is also applied to the radial electron velocity for electron
particles.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  You can use it to assign a temperature :doc:`compute <compute>`
you have defined to this fix which will be used in its thermostatting
procedure, as described above.  For consistency, the group used by
this fix and by the compute should be the same.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change implied by a velocity rescaling to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the cumulative
energy change due to this fix.  The scalar value calculated by this
fix is "extensive".

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the USER-EFF package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix langevin/eff <fix_langevin_eff>`, :doc:`fix nvt/eff <fix_nh_eff>`, :doc:`fix_modify <fix_modify>`,
:doc:`fix temp rescale <fix_temp_rescale>`,

**Default:** none


