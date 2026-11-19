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
* *ne* value = initial electronic coordinate
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
*-dE/dN + mu*.

A single Nose-Hoover chain thermostats the nuclear velocities and the
electron-number velocity together.  The fix creates
:doc:`compute temp/uvt <compute_temp_uvt>` with ID *fix-ID_temp* to provide
the combined temperature and :math:`g+1` degrees of freedom, where
:math:`g` is the nuclear DOF.  The electronic mass is
:math:`W_\mathrm{e}=g k_B T_\mathrm{target} \, \mathrm{Udamp}^2`.

The standard thermo temperature remains nuclear-only unless explicitly
changed.  The first two fix-vector entries report the electron number and combined
temperature, independently of the thermostat chain length. Use::

   thermo_style custom step f_tp[1] temp f_tp[2]

Here *tp* is the fix ID; the columns are timestep, electron number,
nuclear temperature, and combined temperature. No explicit temperature
compute is required. The automatically created compute can alternatively
be referenced as *c_ID_temp*, where *ID* is the fix ID. A temperature selected
with *fix_modify temp*
must be a temp/uvt compute for this fix and group.  The usual
:doc:`compute_modify <compute_modify>` DOF options can be applied to that
compute; the nuclear DOF must remain positive.

The *dedn* source may be an equal-style variable (*v_name*), a
global compute (*c_ID*), or a global fix (*f_ID*).  If the source
provides a global vector, an entry can be selected with the usual
*[index]* syntax.

The *mu* target is constant throughout a run. *Udamp* has time units
and must be positive. Chemical-potential ramps are not supported.

The global fix vector starts with six physical outputs:

#. electron number *Ne*
#. combined instantaneous temperature *T_ins*
#. electron-number velocity *Ne_dot*
#. energy derivative *dEdN*
#. electronic kinetic energy
#. electronic potential contribution *-mu*Ne*

The same six quantities appear in the same order at entries 11 through 16
of :doc:`fix pimd/uvt <fix_pimd>`, following its ten nuclear outputs.
For PIMD, the derivative is reported separately for each bead.

The regular Nose-Hoover vector follows, starting at entry 7. The total
vector length is 6+4*M for a thermostat chain of length M. All six physical
indices are independent of M. The temperature output is intensive and is
not normalized by atom count.

The fix writes restart data for the Nose-Hoover state and electronic
degree of freedom, so simulations may be continued with
:doc:`read_restart <read_restart>` followed by the same *fix uvt*
command.

Restrictions
""""""""""""

This fix supports temperature control only.  Pressure control keywords
are not allowed.

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with both EXTRA-FIX and EXTRA-COMPUTE, which provides
the required :doc:`compute temp/uvt <compute_temp_uvt>`.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix pimd/uvt <fix_pimd>`,
:doc:`read_restart <read_restart>`
