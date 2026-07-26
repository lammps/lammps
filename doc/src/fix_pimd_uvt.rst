.. index:: fix pimd/uvt

fix pimd/uvt command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID pimd/uvt method nmpimd thermostat NHC temp T Tdamp tau_T mu Mu Udamp tau_U ne Ne0 dedn source [keyword value ...]

* ID, group-ID are documented in :doc:`fix <fix>` command
* ``method`` value = ``nmpimd``
* ``thermostat`` value = ``NHC``
* ``temp`` value = target temperature
* ``Tdamp`` value = nuclear thermostat damping parameter
* ``mu`` value = target chemical potential
* ``Udamp`` value = electronic thermostat damping parameter
* ``ne`` value = initial electronic coordinate
* ``dedn`` value = global scalar or vector reference used as dE/dN

Optional keywords include ``integrator``, ``drag``, ``fmass``,
``fmmode``, ``sp``, ``lj``, ``ne_velocity``, ``removecom``,
``tchain``, ``tloop``, and ``tdof``.

Examples
""""""""

.. code-block:: LAMMPS

   variable dEdN equal 3.0*(f_cp[23]-1.2)
   fix cp all pimd/uvt method nmpimd thermostat NHC temp 0.8 Tdamp 0.2 &
          tchain 3 tloop 1 mu 1.5 Udamp 0.2 ne 1.8 dedn v_dEdN

   mpirun -np 2 lmp -partition 2x1 -in in.electron_pimd_uvt

Description
"""""""""""

``fix pimd/uvt`` performs constant-potential path-integral molecular
dynamics.  It extends the :doc:`fix pimd <fix_pimd>` Nose-Hoover
hierarchy with one global electronic coordinate shared by all beads.
The number of path-integral beads is set by the number of LAMMPS
partitions, so the job must be launched with one partition per bead.

The ``dedn`` source may be an equal-style variable (``v_name``), a
global compute (``c_ID``), or a global fix (``f_ID``).  If the source
provides a global vector, an entry can be selected with the usual
``[index]`` syntax.  ``dEdN`` is bead averaged before updating the
global electronic coordinate.

The electronic coordinate is exposed through the global fix vector.
The first 10 entries are the path-integral energy and estimator outputs,
followed by ``4*tchain`` thermostat entries and then nine uVT entries:
``Ne``, ``Ne_dot``, ``dEdN``, ``mu``, the electronic kinetic energy,
the electronic potential contribution, the nuclear thermostat work, the
electronic thermostat work, and their sum.  For ``tchain 3``, ``Ne``,
``Ne_dot``, ``dEdN``, and ``mu`` are ``f_ID[23]`` through ``f_ID[26]``.

The fix writes restart data for the nuclear thermostat state and the
electronic degree of freedom, so simulations may be continued with
:doc:`read_restart <read_restart>` followed by the same ``fix pimd/uvt``
command.

Restrictions
""""""""""""

This fix currently supports only 3d systems.

This fix requires an atom map, e.g. ``atom_modify map yes``.

This fix requires one partition per bead.  In other words, the bead
count must match ``universe->nworlds``.

This fix is a complete time integrator.  No other time-integration fix
should be applied to the same atoms.

Pressure control is not supported by this fix.

This fix is part of the REPLICA package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix pimd <fix_pimd>`, :doc:`fix uvt <fix_uvt>`,
:doc:`read_restart <read_restart>`, :doc:`run_style <run_style>`
