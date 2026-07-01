.. index:: fix tprpmd

fix tprpmd command
==================

.. versionadded:: TBD

Syntax
""""""

.. parsed-literal::

   fix ID group-ID tprpmd method = tp-pimd or tp-rpmd or tprpmd ensemble uvt thermostat NHC temp T Tdamp tau_T \
       tchain Nc tloop Nloop mu Mustart Mustop Mudamp ne Ne0 dedn source [keyword value ...]

where the supported optional keywords are::

   keyword = integrator or drag or fmass or fmmode or sp or lj or tau or ne_velocity or removecom

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all tprpmd method tp-rpmd ensemble uvt thermostat NHC temp 300.0 Tdamp 100.0 \
          tchain 3 tloop 1 mu 0.0 0.0 100.0 ne 1.0 dedn v_dEdN

   mpirun -np 8 lmp -partition 8x1 -in in.electron_tprpmd

Description
"""""""""""

``fix tprpmd`` implements thermostatted ring-polymer molecular dynamics
with a constant-potential electronic degree of freedom.  The path-integral
replicas are distributed across LAMMPS partitions, so the number of beads
is set implicitly by ``universe->nworlds`` and the job must be launched
with the :doc:`-partition <Run_options>` command-line switch.

This fix supports the uVT ensemble and a Nose-Hoover chain thermostat.
The required keywords are:

* ``method`` = ``tp-pimd`` or ``tp-rpmd`` (``tprpmd`` accepted as alias)
* ``ensemble uvt``
* ``thermostat NHC``
* ``temp`` and ``Tdamp`` for the ionic thermostat
* ``tchain`` and ``tloop`` for the thermostat-chain length and sub-cycling
* ``mu``, ``ne``, and ``dedn`` for the constant-potential electronic coordinate

``tp-pimd`` selects thermodynamic path-integral molecular dynamics
with thermostats applied to all normal modes, including the centroid.
``tp-rpmd`` selects thermostatted path-integral RPMD-style dynamics
where the centroid mode is left unthermostatted and thermostats are
applied only to internal modes.  ``tprpmd`` is accepted as an alias
for ``tp-rpmd`` for compatibility.
For ``tp-rpmd`` with a single bead, the unique nuclear centroid mode
remains unthermostatted, while the shared electronic uVT variable
continues to use the Nose-Hoover chain.

The ``dedn`` source may be supplied as an equal-style variable
(``v_name``), a global compute (``c_ID``), or a global fix (``f_ID``).
If the source provides a global vector, an entry can be selected with
the usual ``[index]`` syntax.

The ``integrator`` keyword selects either ``obabo`` or ``baoab``.
The ``fmmode`` and ``fmass`` keywords control the fictitious-mass
preconditioning used for the ring-polymer normal modes.  The ``sp``
keyword scales Planck's constant.  For ``lj`` units, the ``lj`` keyword
supplies the reduced-unit conversion factors in the order
``epsilon sigma mass hplanck mvv2e``.

The electronic coordinate is exposed through the global fix vector.
The first 10 entries are the path-integral energy and estimator outputs,
followed by ``4*tchain`` thermostat entries and then six uVT entries:
``Ne``, ``Ne_dot``, ``dEdN``, ``mu``, the electronic kinetic energy, and
the electronic potential contribution.  Thus, for ``tchain 3``,
``Ne``, ``Ne_dot``, ``dEdN``, and ``mu`` are ``f_ID[23]`` through
``f_ID[26]``.
These six electronic outputs are reported as global values and are not
normalized by the number of atoms when used in thermo output.  The
``dEdN`` entry is the bead-averaged value used to update the global
electronic coordinate.

The fix writes restart data for the thermostat state and electronic
degree of freedom, so simulations may be continued with
:doc:`read_restart <read_restart>` followed by the same ``fix tprpmd``
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

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix pimd <fix_pimd>`, :doc:`fix uvt <fix_uvt>`,
:doc:`read_restart <read_restart>`, :doc:`run_style <run_style>`
