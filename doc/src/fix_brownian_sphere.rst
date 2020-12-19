.. index:: fix brownian/sphere

fix brownian/sphere command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID brownian/sphere gamma_t gamma_r diff_t diff_r seed keyword args

* ID, group-ID are documented in :doc:`fix <fix>` command
* brownian/sphere = style name of this fix command
* gamma_t = translational friction coefficient
* gamma_r = rotational friction coefficient
* diff_t = translational diffusion coefficient
* diff_r = rotational diffusion coefficient
* zero or more keyword/value pairs may be appended
* keyword = *rng* or *dipole*

  .. parsed-literal::

        *rng* value = *uniform* or *gaussian* or *none*
         *uniform* = use uniform random number generator
         *gaussian* = use gaussian random number generator
         *none* = turn off noise
        *dipole* value = none = update orientation of dipoles during integration

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all brownian/sphere 1.0 1.0 1.0 1.0 1294019
   fix 1 all brownian/sphere 1.0 1.0 1.0 1.0 19581092 rng none dipole
   fix 1 all brownian/sphere 1.0 1.0 1.0 1.0 19581092 rng uniform
   fix 1 all brownian/sphere 1.0 1.0 1.0 1.0 19581092 dipole rng gaussian


Description
"""""""""""

Perform Brownian Dynamics integration to update position, velocity, 
angular velocity, and dipole moment for finite-size spherical particles
in the group each timestep. Brownian Dynamics uses Newton's laws of
motion in the limit that inertial forces are negligible compared to
viscous forces. The stochastic equations of motion are

.. math::

   dr = \frac{F}{\gamma_t}dt+\sqrt{2D_t}dW_t, \\
   d\Omega = \frac{T}{\gamma_r}dt + \sqrt{2D_r}dW_r,

where :math:`d\Omega` is an infinitesimal rotation vector (see e.g.
Chapter 4 of :ref:`(Goldstein) <GoldsteinCM>`), :math:`dW_t` and
:math:`dW_r` are Wiener processes (see e.g. :ref:`(Gardiner) <GardinerC>`).
The dipole vectors :math:`e_i` are updated using the rotation matrix

.. math::

   e_i(t+dt) = e^{\theta_X} e_i(t),\\

where :math:`\omega=d\Omega/dt` is the angular velocity,
:math:`\Delta\theta=|\omega|dt` is the rotation angle about
the :math:`\omega` axis, and
:math:`(\theta_X)_{ij}=-\epsilon_{ijk}d\Omega_k` is the
infinitesimal rotation matrix (see e.g. :ref:`(Callegari) <Callegari1>`,
section 7.4).

.. note::
   This integrator is designed for generic non-equilibrium
   simulations with additive noise. There are two important cases which
   (conceptually) reduce the number of free parameters in this fix.
   (a) In equilibrium simulations
   (where fluctuation dissipation theorems are obeyed), one can define
   the thermal energy :math:`k_bT=D_t\gamma_t=D_r\gamma_r`.
   (b) When a no-slip boundary condition is expected between the spheres and
   the surrounding medium, the condition
   :math:`\gamma_t=3\gamma_r/\sigma^2` must be explicitly
   accounted for (e.g. by setting *gamma_t* to 3 and *gamma_r* to 1) where
   :math:`sigma` is the particle diameter.
   If both (a) and (b) are true, then one must ensure this explicitly via
   the above relationships.

---------

.. note::
   Temperature computation using the :doc:`compute temp <compute_temp>`
   will not correctly compute temperature of these overdamped dynamics
   since we are explicitly neglecting inertial effects.
   See e.g. chapter 6 of :ref:`(Doi) <Doi1>` for more details on this.
   Temperature is instead defined in terms of the note above (for
   equilibrium systems).

---------

.. note::
   The diffusion coefficient :math:`D_t` is measured
   in units of (length*length)/time and the diffusion coefficient
   :math:`D_r` is measured in units of 1/time, where time and length
   are in the units specified on the :doc:`units <units>` page. Similarly,
   :math:`\gamma_t` and :math:`\gamma_r` are measured in
   units of mass/time and (mass*length*length)/(time).

---------

If the *rng* keyword is used with the *uniform* value, then the noise
is generated from a uniform distribution (see
:ref:`(Dunweg) <Dunweg6>` for why this works). This is the same method
of noise generation as used in :doc:`fix_langevin <fix_langevin>`.

If the *rng* keyword is used with the *gaussian* value, then the noise
is generated from a gaussian distribution. Typically this added
complexity is unnecessary, and one should be fine using the *uniform*
value for reasons argued in :ref:`(Dunweg) <Dunweg6>`.

If the *rng* keyword is used with the *none* value, then the noise
terms are set to zero.

If the *dipole* keyword is used, then the dipole moments of the particles
are updated as described above.

----------

.. include:: accel_styles.rst

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added forces on atoms to the
system's virial as part of :doc:`thermodynamic output <thermo_style>`.
The default is *virial no*.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix requires that atoms store torque and angular velocity (omega)
as defined by the :doc:`atom_style sphere <atom_style>` command.
If the *dipole* keyword is used, they must also store a dipole moment
as defined by the :doc:`atom_style dipole <atom_style>` command.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.


Related commands
""""""""""""""""

:doc:`fix langevin <fix_langevin>`, :doc:`fix nve/sphere <fix_nve_sphere>`,
:doc:`atom style <atom_style>`

Default
"""""""

The default for *rng* is *uniform*.

----------

.. _GoldsteinCM:

**(Goldstein)** Goldstein, Poole, and Safko, Classical Mechanics, 3rd Ed. (2001).

.. _GardinerC:

**(Gardiner)** Gardiner, A Handbook for the Natural and Social Sciences 4th Ed. (2009).

.. _Callegari1:

**(Callegari)** Callegari and Volpe, *Numerical Simulations of Active Brownian
Particles*, Flowing Matter, 211-238 (2019).

.. _Doi1:

**(Doi)** Doi, Soft Matter Physics (2013).

.. _Dunweg6:

**(Dunweg)** Dunweg and Paul, Int J of Modern Physics C, 2, 817-27 (1991).


