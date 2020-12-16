.. index:: fix bd/asphere

fix bd/asphere command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID bd/asphere gamma_t gamma_r diff_t diff_r seed keyword args

* ID, group-ID are documented in :doc:`fix <fix>` command
* bd/asphere = style name of this fix command
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

   fix 1 all bd/asphere 1.0 1.0 1.0 1.0 1294019
   fix 1 all bd/asphere 1.0 1.0 1.0 1.0 19581092 rng none dipole
   fix 1 all bd/asphere 1.0 1.0 1.0 1.0 19581092 rng uniform
   fix 1 all bd/asphere 1.0 1.0 1.0 1.0 19581092 dipole rng gaussian


Description
"""""""""""

Perform Brownian Dynamics integration to update position, velocity, 
angular velocity, particle orientation, and dipole moment for
finite-size elipsoidal particles in the group each timestep.
Brownian Dynamics uses Newton's laws of
motion in the limit that inertial forces are negligible compared to
viscous forces. The stochastic equations of motion are

.. math::

   dr = \frac{F}{\gamma_t}dt+\sqrt{2D_t}dW_t, \\
   d\Omega = \frac{T}{\gamma_r}dt + \sqrt{2D_r}dW_r,

where :math:`d\Omega` is an infinitesimal rotation vector (see e.g.
Chapter 4 of :ref:`(Goldstein) <GoldsteinCM1>`), :math:`dW_t` and
:math:`dW_r` are Wiener processes (see e.g. :ref:`(Gardiner) <GardinerC1>`).
The quaternions :math:`q` of the ellipsoid are updated each timestep from
the angular velocity vector.

See :doc:`fix bd/sphere <fix_bd_sphere>` for discussion on the
values of :math:`\gamma_t`, :math:`\gamma_r`, :math:`D_t`,
and :math:`D_r` when simulating equilibrium systems.


If the *rng* keyword is used with the *uniform* value, then the noise
is generated from a uniform distribution (see
:ref:`(Dunweg) <Dunweg7>` for why this works). This is the same method
of noise generation as used in :doc:`fix_langevin <fix_langevin>`.

If the *rng* keyword is used with the *gaussian* value, then the noise
is generated from a gaussian distribution. Typically this added
complexity is unnecessary, and one should be fine using the *uniform*
value for reasons argued in :ref:`(Dunweg) <Dunweg7>`.

If the *rng* keyword is used with the *none* value, then the noise
terms are set to zero.

If the *dipole* keyword is used, then the dipole moments of the particles
are updated by setting them along the x axis of the ellipsoidal frames of
reference.

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
as defined by the :doc:`atom_style sphere <atom_style>` command, as well
as atoms which have a definite orientation as defined by the
:doc:`atom_style ellipsoid <atom_style>` command.
Optionally, they can also store a dipole moment as defined by the
:doc:`atom_style dipole <atom_style>` command.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

All particles in the group must be finite-size ellipsoids.  They cannot
be point particles.

Related commands
""""""""""""""""

:doc:`fix bd/sphere <fix_bd_sphere>`, :doc:`fix langevin <fix_langevin>`,
:doc:`fix nve/asphere <fix_nve_asphere>`, :doc:`atom style <atom_style>`

Default
"""""""

The default for *rng* is *uniform*.

----------

.. _GoldsteinCM1:

**(Goldstein)** Goldstein, Poole, and Safko, Classical Mechanics, 3rd Ed. (2001).

.. _GardinerC1:

**(Gardiner)** Gardiner, A Handbook for the Natural and Social Sciences 4th Ed. (2009).

.. _Dunweg7:

**(Dunweg)** Dunweg and Paul, Int J of Modern Physics C, 2, 817-27 (1991).


