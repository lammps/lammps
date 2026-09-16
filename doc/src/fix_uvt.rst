.. index:: fix uvt

fix uvt command
===============

.. versionadded:: TBD

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID uvt temp Tstart Tstop Tdamp mu Mustart Mustop Mudamp ne Ne0 dedn source [keyword value ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* ``temp`` values = Tstart Tstop Tdamp
* ``mu`` values = Mustart Mustop Mudamp
* ``ne`` value = initial electronic coordinate
* ``dedn`` value = global scalar or vector reference used as dE/dN

Optional keywords are the Nose-Hoover keywords supported by
:doc:`fix nvt <fix_nh>`, plus ``ne_velocity``.

Examples
""""""""

.. code-block:: LAMMPS

   variable dEdN equal 5.0*(f_cp[1]-1.0)
   fix cp all uvt temp 1.0 1.0 0.5 mu 2.0 2.0 0.5 ne 1.8 dedn v_dEdN

Description
"""""""""""

``fix uvt`` performs classical molecular dynamics at constant chemical
potential by adding one global electronic coordinate, ``Ne``, to the
standard Nose-Hoover NVT equations of motion.  The electronic coordinate
is propagated with an extended-mass variable and is driven by
``-dE/dN + mu``.

A single Nose-Hoover chain thermostats the nuclear velocities and the
electron-number velocity together.  The fix creates
:doc:`compute temp/uvt <compute_temp_uvt>` with ID *fix-ID_temp* to provide
the combined temperature and :math:`g+1` degrees of freedom, where
:math:`g` is the nuclear DOF.  The electronic mass remains
:math:`W_\mathrm{e}=g k_B T_\mathrm{target} \, \mathrm{Mudamp}^2`.
It is initialized before the first combined temperature evaluation and
updated with the temperature target during integration.

The standard thermo temperature remains nuclear-only unless explicitly
changed.  The first two fix-vector entries report the electron number and combined
temperature, independently of the thermostat chain length. Use::

   thermo_style custom step f_tp[1] temp f_tp[2]

Here ``tp`` is the fix ID; the columns are timestep, electron number,
nuclear temperature, and combined temperature. No explicit temperature
compute is required. The automatically created compute can alternatively
be referenced as *c_fix-ID_temp*.  A temperature selected with *fix_modify temp*
must be a temp/uvt compute for this fix and group.  The usual
:doc:`compute_modify <compute_modify>` DOF options can be applied to that
compute; the nuclear DOF must remain positive.

The ``dedn`` source may be an equal-style variable (``v_name``), a
global compute (``c_ID``), or a global fix (``f_ID``).  If the source
provides a global vector, an entry can be selected with the usual
``[index]`` syntax.

The derivative is read once during setup, after the initial force calculation,
and then in ``post_force`` after each force calculation.  The cached value
is used by the second electronic velocity half-step and the first half-step
of the next timestep.  This applies to both analytical variables and
sources whose derivative is produced during force evaluation.  There is no
``dedn_defer`` option.

If the source is produced by another fix in its ``post_force`` callback,
define that fix before ``fix uvt`` so it updates the derivative first.
The provider must also make the initial derivative available during its
``setup`` callback.  The same ordering requirement applies to indirect
compute or variable references to such a fix.

With r-RESPA, the derivative is refreshed and the electronic velocity is
kicked at the outermost level; the electronic coordinate drifts with the
particle coordinates at the innermost level.  Derivative providers must
supply the complete derivative at the outermost force stage.

The global fix vector starts with seven physical outputs:

#. electron number ``Ne``
#. combined instantaneous temperature ``T_ins``
#. electron-number velocity ``Ne_dot``
#. energy derivative ``dEdN``
#. electrochemical-potential target ``mu``
#. electronic kinetic energy
#. electronic potential contribution ``-mu*Ne``

The regular Nose-Hoover vector follows, starting at entry 8. The total
vector length is 7+4*M for a thermostat chain of length M. All seven physical
indices are independent of M. The temperature output is intensive and is
not normalized by atom count.

This replaces the earlier layout that appended the electronic outputs after
the chain. Update input scripts accordingly, including derivative-source
references to the electron number (formerly entry 13 for ``tchain 3``, now
entry 1).

The first Nose-Hoover potential-energy entry includes all
:math:`g+1` controlled DOF.  The fix scalar includes the combined chain
energy plus the electronic kinetic and potential terms, without counting
the electronic thermostat contribution twice.

The fix writes restart data for the Nose-Hoover state and electronic
degree of freedom, so simulations may be continued with
:doc:`read_restart <read_restart>` followed by the same ``fix uvt``
command.

Restrictions
""""""""""""

This fix supports temperature control only.  Pressure control keywords
are not allowed.

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix pimd/uvt <fix_pimd_uvt>`,
:doc:`read_restart <read_restart>`
