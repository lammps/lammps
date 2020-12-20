.. index:: fix brownian

fix brownian command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID brownian gamma_t diff_t seed keyword args

* ID, group-ID are documented in :doc:`fix <fix>` command
* brownian/sphere = style name of this fix command
* gamma_t = translational friction coefficient
* diff_t = translational diffusion coefficient
* zero or more keyword/value pairs may be appended
* keyword = *rng* 

  .. parsed-literal::

        *rng* value = *uniform* or *gaussian* or *none*
         *uniform* = use uniform random number generator
         *gaussian* = use gaussian random number generator
         *none* = turn off noise

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all brownian 1.0 3.0 1294019
   fix 1 all brownian 1.0 3.0 19581092 rng none 
   fix 1 all brownian 1.0 3.0 19581092 rng uniform
   fix 1 all brownian 1.0 3.0 19581092 rng gaussian


Description
"""""""""""

Perform Brownian Dynamics integration to update position and velocity
of atoms in the group each timestep. Brownian Dynamics uses Newton's laws of
motion in the limit that inertial forces are negligible compared to
viscous forces. The stochastic equation of motion is

.. math::

   dr = \frac{F}{\gamma_t}dt+\sqrt{2D_t}dW_t, \\

where :math:`dW_t` is a Wiener processes (see e.g. :ref:`(Gardiner) <GardinerC1>`).

.. note::
   This integrator is designed for generic non-equilibrium
   simulations with additive noise. There are two important cases which
   (conceptually) reduce the number of free parameters in this fix.
   (a) In equilibrium simulations
   (where fluctuation dissipation theorems are obeyed), one can define
   the thermal energy :math:`k_bT=D_t\gamma_t`.
   
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
   in units of (length*length)/time, where time and length
   are in the units specified on the :doc:`units <units>` page.
   Similarly, :math:`\gamma_t` is measured in
   units of mass/time.

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

This fix is part of the USER-BROWNIAN package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.


Related commands
""""""""""""""""

:doc:`fix langevin <fix_langevin>`, :doc:`fix nve/sphere <fix_nve_sphere>`,
:doc:`fix brownian/sphere <fix_brownian_sphere>`,
:doc:`fix brownian/asphere <fix_brownian_asphere>`

Default
"""""""

The default for *rng* is *uniform*.

----------

.. _GardinerC1:

**(Gardiner)** Gardiner, A Handbook for the Natural and Social Sciences 4th Ed. (2009).

.. _Doi1:

**(Doi)** Doi, Soft Matter Physics (2013).

.. _Dunweg6:

**(Dunweg)** Dunweg and Paul, Int J of Modern Physics C, 2, 817-27 (1991).


