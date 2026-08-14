.. index:: fix baoab

fix baoab command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID baoab Tstart Tstop damp seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* baoab = style name of this fix command
* Tstart, Tstop = desired temperature at start/end of run (temperature units)
* damp = damping parameter (time units)
* seed = random number seed (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *zero*

  .. parsed-literal::

       *zero* value = *no* or *yes*
         *no* = do not zero the net random momentum each step (default)
         *yes* = zero the net random momentum each step

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all baoab 300.0 300.0 100.0 12345
   fix 1 all baoab 1.0 1.0 10.0 48279 zero yes
   fix 2 all baoab 300.0 400.0 200.0 77777

.. versionadded:: 4Jul2026

Description
"""""""""""

Apply the BAOAB Langevin integrator of :ref:`(Leimkuhler) <Leimkuhler1>` to
the atoms in the fix group.  This fix performs **complete** time integration
of Newton's equations of motion coupled to a Langevin heat bath -- do
**not** combine it with :doc:`fix nve <fix_nve>`.

The Langevin equation of motion integrated by this fix is

.. math::

   m\,d\mathbf{v} = \mathbf{F}_c\,dt
   - \frac{m}{\tau}\,\mathbf{v}\,dt
   + \sqrt{\frac{2\,m\,k_B T}{\tau}}\,d\mathbf{W}_t,

where :math:`\mathbf{F}_c` is the conservative force, :math:`\tau` is the
damping time (*damp*), :math:`T` is the target temperature, and
:math:`d\mathbf{W}_t` is a Wiener increment.

The BAOAB splitting implements the following sub-steps for each MD timestep:

* **B** -- half-step velocity kick from conservative forces
* **A** -- half-step position drift
* **O** -- full-step exact Ornstein-Uhlenbeck (OU) update (thermostat)
* **A** -- half-step position drift
* **B** -- half-step velocity kick from conservative forces (using new forces)

The O step applies the exact solution of the OU process,

.. math::

   \mathbf{v} \leftarrow e^{-\gamma\,dt}\,\mathbf{v}
   + \sqrt{\frac{k_B T}{m}\bigl(1 - e^{-2\gamma dt}\bigr)}\;\mathbf{R},
   \qquad \gamma = 1/\mathrm{damp},

where :math:`\mathbf{R}` is a vector of independent standard Gaussian
random numbers.  This is exact for any timestep, unlike a first-order
Euler discretization of the friction and noise.

BAOAB achieves second-order accuracy in configuration space: the
stationary distribution of positions converges to the Boltzmann
distribution with error :math:`O(dt^2)`, making it significantly more
accurate than first-order Langevin integrators (such as the
Euler-Maruyama scheme) at the same timestep.  This is particularly
important for computing configurational averages, free energies, and
structural properties.

The target temperature can be ramped linearly from *Tstart* to *Tstop*
over the course of a :doc:`run <run>` by using the *start* and *stop*
keywords of the :doc:`run <run>` command.

The *damp* parameter has units of time and sets the relaxation time of
the Langevin thermostat.  Smaller values couple the system more strongly
to the bath (faster thermalization, more friction); larger values give
weaker coupling (less friction, closer to NVE dynamics).  A typical
starting value is 100 times the MD timestep.

If the *zero* keyword is set to *yes*, the net random momentum injected
by the O step is subtracted each timestep so that the total momentum of
the fix group does not drift.  This is useful for bulk systems without
periodic boundary conditions or when the group does not span the full
system.

The cumulative energy exchanged between the atoms and the Langevin
reservoir is always tracked and available via the :doc:`ecouple thermo
<thermo_style>` keyword.  The sign convention is that a positive value
means energy has flowed from the system into the reservoir (cooling),
and negative means the reservoir has heated the system.  Correspondingly
the ``econserve`` thermo keyword reports the sum of ``etotal`` and
``ecouple`` and thus the conserved quantity of a constant temperature
simulation.

This fix supports both per-type masses (``mass`` command) and per-atom
masses (atom styles such as ``sphere``).

This fix can be used with dynamic atom groups.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The cumulative thermostat energy is written to and read from
:doc:`binary restart files <restart>` so that energy accounting is
continuous across restarts.

This fix produces a global scalar (the cumulative thermostat energy) that
is extensive and accessible to :doc:`thermo_style custom <thermo_style>`
via the ``f_ID`` keyword.

This fix can ramp its target temperature over multiple runs using the
*start* and *stop* keywords of the :doc:`run <run>` command.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This fix cannot be combined with :doc:`fix shake <fix_shake>` or
:doc:`fix rattle <fix_shake>`.

Do not combine this fix with :doc:`fix nve <fix_nve>` or any other
time-integration fix on the same group of atoms.  LAMMPS will print
a warning in this case.

Related commands
""""""""""""""""

:doc:`fix langevin <fix_langevin>`,
:doc:`fix nvt <fix_nh>`,
:doc:`fix gjf <fix_gjf>`,
:doc:`fix gle <fix_gle>`,
:doc:`fix gld <fix_gld>`

Default
"""""""

The option defaults are zero = no.

----------

.. _Leimkuhler1:

**(Leimkuhler)** B. Leimkuhler and C. Matthews, "Rational Construction of
Stochastic Numerical Methods for Molecular Sampling", Appl. Math. Res.
Express, 2013(1), 34-56 (2013).
https://doi.org/10.1093/amrx/abs010
