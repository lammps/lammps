.. index:: fix uvt

fix uvt command
===============

.. versionadded:: TBD

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID uvt temp Tstart Tstop Tdamp mu Mu Udamp ne Ne0 dedn source [keyword value ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* *temp* values = Tstart Tstop Tdamp (temperature, temperature, and time units)
* *mu* value = constant chemical potential
* *Udamp* value = positive electronic damping time (time units)
* *ne* value = initial electron number
* *dedn* value = global scalar or vector reference used as dE/dN

Optional keywords are the Nose-Hoover keywords supported by
:doc:`fix nvt <fix_nh>`. The additional *ne_velocity* keyword sets the
initial electron-number velocity (electron number per unit time); its
default is zero.

Examples
""""""""

.. code-block:: LAMMPS

   variable dEdN equal 5.0*(f_cp[1]-1.0)
   fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 Udamp 0.5 ne 1.8 dedn v_dEdN

Description
"""""""""""

*fix uvt* performs classical molecular dynamics at constant chemical
potential by adding one global electronic coordinate, *Ne*, to the
standard Nose-Hoover NVT equations of motion.  The electronic coordinate
is propagated with an extended-mass variable and is driven by
*-dE/dN + mu*.  The method is described in :ref:`(Fu2026) <uvt-Fu2026>`.

A single Nose-Hoover chain thermostats the nuclear velocities and the
electron-number velocity together.  The fix creates
:doc:`compute temp/uvt <compute_temp_uvt>` with ID *fix-ID_temp* to provide
the combined temperature and :math:`g+1` degrees of freedom, where
:math:`g` is the nuclear DOF.  The electronic mass is
:math:`W_\mathrm{e}=g k_B T_\mathrm{target} \, \mathrm{Udamp}^2`.

The standard thermo temperature remains nuclear-only unless explicitly
changed.  The first two entries of the global vector of this fix report
the electron number and the combined temperature, so no explicit
temperature compute is required to monitor them.  For the example
above::

   thermo_style custom step f_cp[1] temp f_cp[2]

prints the timestep, the electron number, the nuclear temperature, and
the combined temperature.

The *dedn* source may be an equal-style variable (*v_name*), a
global compute (*c_ID*), or a global fix (*f_ID*).  If the source
provides a global vector, an entry can be selected with the usual
*[index]* syntax.

The *mu* target is constant throughout a run. *Udamp* has time units
and must be positive. Chemical-potential ramps are not supported.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the Nose-Hoover thermostat and of the
electronic degree of freedom to :doc:`binary restart files <restart>`.
See the :doc:`read_restart <read_restart>` command for info on how to
re-specify a fix in an input script that reads a restart file, so that
the operation of the fix continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *temp* option is supported by this
fix.  The selected compute must be a :doc:`compute temp/uvt
<compute_temp_uvt>` that references this fix and uses the same group.
The automatically created compute can also be referenced directly as
*c_ID_temp*, where *ID* is the fix ID, and the usual
:doc:`compute_modify <compute_modify>` options for degrees of freedom
can be applied to it; the nuclear DOF must remain positive.

The cumulative energy change in the system imposed by this fix, via
thermostatting and the electronic degree of freedom, is included in the
:doc:`thermodynamic output <thermo_style>` keywords *ecouple* and
*econserve*.  See the :doc:`thermo_style <thermo_style>` page for
details.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the thermostat
energy of :doc:`fix nvt <fix_nh>` plus the electronic kinetic energy and
the electronic potential contribution *-mu*Ne*.  The scalar value
calculated by this fix is "extensive".

This fix also computes a global vector of quantities, which can be
accessed by various :doc:`output commands <Howto_output>`.  The vector
values are "intensive".  The first six entries are:

#. electron number *Ne*
#. combined instantaneous temperature *T_ins*
#. electron-number velocity *Ne_dot*
#. energy derivative *dEdN*
#. electronic kinetic energy
#. electronic potential contribution *-mu*Ne*

They are followed by the Nose-Hoover thermostat entries of :doc:`fix nvt
<fix_nh>`, starting at entry 7, so the vector has 6+4*M entries for a
thermostat chain of length M.  The first six entries do not depend on
M.  The same six quantities appear in the same order at entries 11
through 16 of :doc:`fix pimd/uvt <fix_pimd>`.

This fix supports automatically generated thermo column names when using
:doc:`thermo_modify colname auto <thermo_modify>`.  The names of the
first six entries are "f\_", followed by the fix ID, followed by a
colon and one of *Ne*, *T_ins*, *Ne_dot*, *dEdN*, *ke_Ne*, and *pe_mu*.
The thermostat entries use the names of :doc:`fix nvt <fix_nh>`.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the :doc:`run <run>` command.  See the
:doc:`run <run>` command for details of how to do this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix supports temperature control only.  Pressure control keywords
are not allowed.

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  It also requires the EXTRA-COMPUTE
package, because it creates a :doc:`compute temp/uvt <compute_temp_uvt>`
for its thermostat.  See the :doc:`Build package <Build_package>` page
for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix pimd/uvt <fix_pimd>`,
:doc:`compute temp/uvt <compute_temp_uvt>`,
:doc:`read_restart <read_restart>`

Default
"""""""

The default for *ne_velocity* is 0.0.  The optional Nose-Hoover keywords
have the same defaults as for :doc:`fix nvt <fix_nh>`.

----------

.. _uvt-Fu2026:

**(Fu2026)** L. Fu, Y. Li, and S. Xu, "fix uvt and fix pimd/uvt: A unified
LAMMPS framework for constant-potential constant-temperature molecular
dynamics," arXiv:2609.21935 (2026).
