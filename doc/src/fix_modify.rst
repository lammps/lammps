.. index:: fix_modify

fix_modify command
==================

Syntax
""""""


.. code-block:: LAMMPS

   fix_modify fix-ID keyword value ...

* fix-ID = ID of the fix to modify
* one or more keyword/value pairs may be appended
* keyword = *temp* or *press* or *energy* or *virial* or *respa* or *dynamic/dof* or *bodyforces*
  
  .. parsed-literal::
  
       *temp* value = compute ID that calculates a temperature
       *press* value = compute ID that calculates a pressure
       *energy* value = *yes* or *no*
       *virial* value = *yes* or *no*
       *respa* value = *1* to *max respa level* or *0* (for outermost level)
       *dynamic/dof* value = *yes* or *no*
         yes/no = do or do not re-compute the number of degrees of freedom (DOF) contributing to the temperature
       *bodyforces* value = *early* or *late*
         early/late = compute rigid-body forces/torques early or late in the timestep



Examples
""""""""


.. code-block:: LAMMPS

   fix_modify 3 temp myTemp press myPress
   fix_modify 1 energy yes
   fix_modify tether respa 2

Description
"""""""""""

Modify one or more parameters of a previously defined fix.  Only
specific fix styles support specific parameters.  See the doc pages
for individual fix commands for info on which ones support which
fix\_modify parameters.

The *temp* keyword is used to determine how a fix computes
temperature.  The specified compute ID must have been previously
defined by the user via the :doc:`compute <compute>` command and it must
be a style of compute that calculates a temperature.  All fixes that
compute temperatures define their own compute by default, as described
in their documentation.  Thus this option allows the user to override
the default method for computing T.

The *press* keyword is used to determine how a fix computes pressure.
The specified compute ID must have been previously defined by the user
via the :doc:`compute <compute>` command and it must be a style of
compute that calculates a pressure.  All fixes that compute pressures
define their own compute by default, as described in their
documentation.  Thus this option allows the user to override the
default method for computing P.

The *energy* keyword can be used with fixes that support it.
*energy yes* adds a contribution to the potential energy of the
system. The fix's global and per-atom
energy is included in the calculation performed by the :doc:`compute pe <compute_pe>` or :doc:`compute pe/atom <compute_pe_atom>`
commands.  See the :doc:`thermo_style <thermo_style>` command for info
on how potential energy is output.  For fixes that tally a global
energy, it can be printed by using the keyword f\_ID in the
thermo\_style custom command, where ID is the fix-ID of the appropriate
fix.

.. note::

   You must also specify the *energy yes* setting for a fix if you
   are using it when performing an :doc:`energy minimization <minimize>`
   and if you want the energy and forces it produces to be part of the
   optimization criteria.

The *virial* keyword can be used with fixes that support it.
*virial yes* adds a contribution to the virial of the
system. The fix's global and per-atom
virial is included in the calculation performed by the :doc:`compute pressure <compute_pressure>` or
:doc:`compute stress/atom <compute_stress_atom>`
commands.  See the :doc:`thermo_style <thermo_style>` command for info
on how pressure is output.

.. note::

   You must specify the *virial yes* setting for a fix if you
   are doing :doc:`box relaxation <fix_box_relax>` and
   if you want virial contribution of the fix to be part of the
   relaxation criteria, although this seems unlikely.

.. note::

   This option is only supported by fixes that explicitly say
   so. For some of these (e.g. the :doc:`fix shake <fix_shake>` command)
   the default setting is *virial yes*\ , for others it is *virial no*\ .

For fixes that set or modify forces, it may be possible to select at
which :doc:`r-RESPA <run_style>` level the fix operates via the *respa*
keyword. The RESPA level at which the fix is active can be selected.
This is a number ranging from 1 to the number of levels. If the RESPA
level is larger than the current maximum, the outermost level will be
used, which is also the default setting. This default can be restored
using a value of *0* for the RESPA level. The affected fix has to be
enabled to support this feature; if not, *fix\_modify* will report an
error. Active fixes with a custom RESPA level setting are reported
with their specified level at the beginning of a r-RESPA run.

The *dynamic/dof* keyword determines whether the number of atoms N in
the fix group and their associated degrees of freedom are re-computed
each time a temperature is computed.  Only fix styles that calculate
their own internal temperature use this option.  Currently this is
only the :doc:`fix rigid/nvt/small <fix_rigid>` and :doc:`fix rigid/npt/small <fix_rigid>` commands for the purpose of
thermostatting rigid body translation and rotation.  By default, N and
their DOF are assumed to be constant.  If you are adding atoms or
molecules to the system (see the :doc:`fix pour <fix_pour>`, :doc:`fix deposit <fix_deposit>`, and :doc:`fix gcmc <fix_gcmc>` commands) or
expect atoms or molecules to be lost (e.g. due to exiting the
simulation box or via :doc:`fix evaporate <fix_evaporate>`), then
this option should be used to insure the temperature is correctly
normalized.

.. note::

   Other thermostatting fixes, such as :doc:`fix nvt <fix_nh>`, do
   not use the *dynamic/dof* keyword because they use a temperature
   compute to calculate temperature.  See the :doc:`compute_modify dynamic/dof <compute_modify>` command for a similar way to insure
   correct temperature normalization for those thermostats.

The *bodyforces* keyword determines whether the forces and torques
acting on rigid bodies are computed *early* at the post-force stage of
each timestep (right after per-atom forces have been computed and
communicated among processors), or *late* at the final-integrate stage
of each timestep (after any other fixes have finished their post-force
tasks).  Only the rigid-body integration fixes use this option, which
includes :doc:`fix rigid <fix_rigid>` and :doc:`fix rigid/small <fix_rigid>`, and their variants, and also :doc:`fix poems <fix_poems>`.

The default is *late*\ .  If there are other fixes that add forces to
individual atoms, then the rigid-body constraints will include these
forces when time-integrating the rigid bodies.  If *early* is
specified, then new fixes can be written that use or modify the
per-body force and torque, before time-integration of the rigid bodies
occurs.  Note however this has the side effect, that fixes such as
:doc:`fix addforce <fix_addforce>`, :doc:`fix setforce <fix_setforce>`,
:doc:`fix spring <fix_spring>`, which add forces to individual atoms
will have no effect on the motion of the rigid bodies if they are
specified in the input script after the fix rigid command.  LAMMPS
will give a warning if that is the case.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`fix <fix>`, :doc:`compute temp <compute_temp>`, :doc:`compute pressure <compute_pressure>`, :doc:`thermo_style <thermo_style>`

Default
"""""""

The option defaults are temp = ID defined by fix, press = ID defined
by fix, energy = no, virial = different for each fix style, respa = 0,
bodyforce = late.
