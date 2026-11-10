.. index:: fix langevin/tspin

fix langevin/tspin command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID langevin/tspin Tstart Tstop Tdamp seed keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin/tspin = style name of this fix command
* Tstart, Tstop = desired spin temperature at start and end of the run
  (temperature units, K in metal units)
* Tdamp = spin temperature damping parameter (time units)
* seed = random number seed to use for the white noise (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *zero*

  .. parsed-literal::

       *zero* value = *no* or *yes*
         no = the net random force on the spin system need not be zero
         yes = subtract the mean random force so that it sums to zero

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all langevin/tspin 300.0 300.0 0.05 48279
   fix 2 all langevin/tspin 10.0 600.0 0.1 48279 zero yes

Description
"""""""""""

.. versionadded:: TBD

Used with :doc:`fix nve/tspin <fix_nve_tspin>`, this fix generates
inertial Langevin dynamics of the full Cartesian spin vectors for
equilibrium sampling.  These auxiliary dynamics do not represent physical
spin precession.  Friction and random forces modify the generalized spin
force on each atom in the group with positive spin mass,

.. math::

   \vec{F}^{m}_i \rightarrow \vec{F}^{m}_i
   - \frac{m_{s,i}}{T_{\mathrm{damp}}}\vec{v}^{s}_i
   + \sqrt{\frac{2m_{s,i}k_B T_{\mathrm{target}}(t)}
                 {T_{\mathrm{damp}}\,\Delta t}}\,\boldsymbol{\xi}_i.

Here :math:`m_{s,i}` is the spin mass, :math:`\Delta t` is the timestep,
and :math:`\boldsymbol{\xi}_i` is a dimensionless random vector.
The target temperature :math:`T_{\mathrm{target}}(t)` varies linearly
from *Tstart* to *Tstop*.  With *zero yes*, the mean random force over
the active spins is subtracted after sampling.

Only spin velocities are thermostatted.  Atomic velocities can be
thermostatted separately if needed.  Do not thermostat the same spins
with :doc:`fix nvt/tspin or fix npt/tspin <fix_nvt_tspin>`, or combine
this bath with :doc:`fix langevin/spin <fix_langevin_spin>`.  See
:doc:`fix nve/tspin <fix_nve_tspin>` for the longitudinal-potential
requirement.

For example, after defining the atoms, atomic masses, and initial spins,
and assigning spin masses:

.. code-block:: LAMMPS

   fix 1 all spring/tspin 1.0 2.2
   fix 2 all langevin/tspin 300.0 300.0 0.05 48279
   fix 3 all nve/tspin lattice frozen

Report the auxiliary kinetic spin temperature and spin kinetic energy
with :doc:`compute temp/tspin <compute_temp_tspin>`.  The bath exchanges
energy with the spins, so *etotal* plus the spin kinetic energy is not
conserved.  This fix does not tally that exchange.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The random-number generator state is not written to :doc:`binary
restart files <restart>`.  Restarts preserve statistical sampling but
do not reproduce the same trajectory.  Each MPI rank uses a separate
random-number stream, so changing the number of ranks changes the
trajectory.

No :doc:`fix_modify <fix_modify>` options or :doc:`output quantities
<Howto_output>` are supported.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The *langevin/tspin* fix is part of the SPIN package.  This style is
only enabled if LAMMPS was built with this package.  See the
:doc:`Build package <Build_package>` page for more info.

This fix requires :doc:`atom_style spin <atom_style>`.  The time
integration has to be performed by :doc:`fix nve/tspin
<fix_nve_tspin>`.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`compute temp/tspin <compute_temp_tspin>`,
:doc:`fix langevin <fix_langevin>`

Default
""""""""

The default value for the *zero* keyword is *no*.
