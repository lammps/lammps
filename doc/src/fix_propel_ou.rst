.. index:: fix propel/ou

fix propel/ou command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID propel/ou magnitude tau seed

* ID, group-ID are documented in :doc:`fix <fix>` command
* propel/ou = style name of this fix command
* magnitude = root mean square value of each component of the active force (force units)
* tau = correlation time of the active force (time units)
* seed = random number seed to use for the active force (positive integer)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all propel/ou 2.0 0.5 12345

Description
"""""""""""

.. versionadded:: TBD

Add a fluctuating active force :math:`\mathbf{f}^a_i` to each atom in the
group, whose components are independent Ornstein-Uhlenbeck processes.
Particles driven by such a force are known as active Ornstein-Uhlenbeck
particles (AOUP) :ref:`(Szamel2014) <Szamel2014>`, :ref:`(Fodor2016)
<Fodor2016>`, :ref:`(Martin2021) <Martin2021>`, a popular minimal model
of active matter, since the persistent force is Gaussian and can be
treated analytically.  In contrast to :doc:`fix propel/self
<fix_propel_self>`, both the direction and the magnitude of the active
force fluctuate in time, and no orientation vector is required, so this
fix can be used with any atom style.  See the :doc:`Howto active matter
<Howto_active>` page for an overview of active matter models in LAMMPS.

Each Cartesian component :math:`f` of the active force of an atom
follows the stochastic differential equation

.. math::

   \tau \frac{d f}{d t} = - f + f_0 \sqrt{2 \tau}\, \xi(t)

where :math:`\tau` is the correlation time *tau*, :math:`f_0` is the
*magnitude*, and :math:`\xi(t)` is Gaussian white noise with zero mean
and :math:`\langle \xi(t) \xi(t') \rangle = \delta(t - t')`.  In the
stationary state, each force component is Gaussian distributed with zero
mean and the root mean square value :math:`f_0`, and the components are
exponentially correlated in time:

.. math::

   \langle f_\alpha(t) f_\beta(t') \rangle = f_0^2\, \delta_{\alpha\beta}\, e^{-|t-t'|/\tau}

In 2d simulations only the x and y components of the active force are
non-zero.

The active force is advanced by one timestep using the exact solution of
the Ornstein-Uhlenbeck process over a time interval :math:`\Delta t`
:ref:`(Gillespie1996) <Gillespie1996>`,

.. math::

   f(t + \Delta t) = f(t)\, e^{-\Delta t/\tau} + f_0 \sqrt{1 - e^{-2\Delta t/\tau}}\; \eta

with :math:`\eta` a Gaussian random number with zero mean and unit
variance, so that the statistical properties of the active force do not
depend on the size of the timestep.  The initial active force of each
atom in the group is drawn from the stationary distribution when the fix
is defined.  Atoms created later (e.g. by :doc:`create_atoms
<create_atoms>` or :doc:`fix deposit <fix_deposit>`) receive an initial
active force in the same way if they belong to the fix group at the time
of their creation, and otherwise start with an active force of zero.

When the translational motion is overdamped, e.g. when using :doc:`fix
brownian <fix_brownian>` with the translational friction coefficient
:math:`\gamma_t`, the active force results in an active velocity
:math:`\mathbf{v}^a = \mathbf{f}^a / \gamma_t` with the root mean square
value :math:`f_0/\gamma_t` per component and the persistence time
:math:`\tau`.  At times much longer than :math:`\tau`, the active force
then contributes an active diffusion coefficient :math:`D_a = f_0^2
\tau/\gamma_t^2` to the diffusion of the particles, i.e. their mean
squared displacement grows as :math:`2 d\, (D_t + D_a)\, t` with
:math:`D_t` the thermal diffusion coefficient and :math:`d` the
dimensionality of the system.

Along with adding a force contribution, this fix can also contribute to
the virial (pressure) of the system, defined as :math:`\sum_i \langle
\mathbf{f}^a_i \cdot \mathbf{r}_i \rangle/(d V)`, where
:math:`\mathbf{r}_i` is the *unwrapped* coordinate of particle *i* in
the case of periodic boundary conditions.  See the discussion of the
active pressure in the documentation of :doc:`fix propel/self
<fix_propel_self>` for the limitations of this definition.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the per-atom active forces and the state of its random
number generator to :doc:`binary restart files <restart>`, so that a
simulation continued from a restart reproduces the same stochastic
trajectory (provided the number of MPI processes is unchanged).

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added forces on atoms to the
system's virial as part of :doc:`thermodynamic output <thermo_style>`.
The default is *virial yes*.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix.  This allows to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its forces.  Default is the
outermost level.

This fix computes a per-atom array with 3 columns, which contains the x,
y, and z components of the current active force of each atom in force
units, and which can be accessed by various :doc:`output commands
<Howto_output>`.  The values are zero for atoms outside the fix group.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the BROWNIAN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix propel/self <fix_propel_self>`, :doc:`fix brownian <fix_brownian>`,
:doc:`fix langevin <fix_langevin>`, :doc:`fix tumble <fix_tumble>`

Default
"""""""

none

----------

.. _Szamel2014:

**(Szamel2014)** G. Szamel, Self-propelled particle in an external potential: Existence of an effective temperature, Phys. Rev. E 90, 012111 (2014).

.. _Fodor2016:

**(Fodor2016)** E. Fodor, C. Nardini, M. E. Cates, J. Tailleur, P. Visco, and F. van Wijland, How Far from Equilibrium Is Active Matter?, Phys. Rev. Lett. 117, 038103 (2016).

.. _Martin2021:

**(Martin2021)** D. Martin, J. O'Byrne, M. E. Cates, E. Fodor, C. Nardini, J. Tailleur, and F. van Wijland, Statistical mechanics of active Ornstein-Uhlenbeck particles, Phys. Rev. E 103, 032607 (2021).

.. _Gillespie1996:

**(Gillespie1996)** D. T. Gillespie, Exact numerical simulation of the Ornstein-Uhlenbeck process and its integral, Phys. Rev. E 54, 2084 (1996).
