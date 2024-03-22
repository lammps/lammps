.. index:: fix hmc

fix hmc command
===============

Syntax
""""""
.. code-block:: LAMMPS

   fix ID group-ID hmc N seed temp integrator keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* hmc = style name of this fix command
* N = invoke a Monto Carlo step every N steps
* seed = random # seed (positive integer)
* temp = temperature for assigning velocities
* integrator = integrator fix: flexible (for nve) or rigid (for rigid/small)
* keyword = *mom* or *ra*

  .. parsed-literal::

       *mom* value = *yes* or *no*
       *ra* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all hmc 10 123 500 flexible
   fix hmc_water water hmc 100 123 298.15 rigid
   fix 2 all hmc 10 12345 300 flexible mom no ra yes

Description
"""""""""""

.. versionadded:: TBD

This fix performs the the Hybrid/Hamiltonian Monte Carlo (HMC) algorithm
in line with the following order of steps:

The new particle configuration (positions and velocities) is calculated
by invoking the velocity-Verlet time integration algorithm.
Before these configuration changes are performed, the proposed change
in the Hamiltonian, :math:`\Delta{H}` is calculated following the equation:

.. math::

   \Delta{H} = H(q',p') -  H(q,p)

This new proposed configuration is then accepted/rejected according to
the Metropolis criterion with probability:

.. math::

   p^{acc} = min(1,e^{\frac{-\Delta{H}}{k_B T}})

Upon acceptance, the new proposed particle configuration positions and
velocities are updated. Upon rejection, the old particle configuration
is kept, and particle momenta (and therefore velocities) are randomly
resampled from a normal distribution:

.. math::

   p_{x,y,z} = \textbf{N}(0,1) \sqrt{\frac{k_B T}{2 m^2}}

The algorithm then continues, proposing a new configuration of particles
and velocities N integration steps later.

----------

The keyword/value options are used in the following ways.

The *mom* keyword sets the linear momentum of the ensemble of particles.
If mom = yes, the linear momentum of the ensemble of velocities is
zeroed. If mom = no, the linear momentum of the ensemble of velocities
is not zeroed.

The *ra* keyword decides whether velocities are resampled upon acceptance.
If ra = yes, velocities are resampled upon acceptance. If ra = no,
velocities are not resampled upon acceptance.

----------

Ouput info
""""""""""

This fix computes a global scalar and global vector of length 5,
which can be accessed by various :doc:`output commands
<Howto_output>`.  The scalar is the fraction of attempted MC moves which have been accepted.  The vector stores the
following quantities:

* 1 = number of accepted moves
* 2 = number of rejected moves
* 3 = change in potential energy
* 4 = change in kinetic energy
* 5 = change in total energy (kinetic + potential energy)

These values are calculated every N timesteps

Restrictions
""""""""""""

This fix is part of the MC package and requires the RIGID package to
be installed. It is only enabled if LAMMPS was built with both packages.
See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix gcmc <fix_gcmc>`, :doc:`fix tfmc <fix_tfmc>`

Default
"""""""

The option default is mom = yes, ra = no.

----------

**(Watkins)** Watkins and Jorgensen, J Phys Chem A, 105, 4118-4125 (2001).

**(Betancourt)** Betancourt, M. A Conceptual Introduction to Hamiltonian Monte Carlo, 2018.

**(Duane)** Duane, S.; Kennedy, A. D.; Pendleton, B. J.; Roweth, D. Hybrid Monte Carlo. Physics Letters B 1987, 195 (2), 216-222. https://doi.org/10.1016/0370-2693(87)91197-X.

**(Metropolis)** Metropolis, N.; Rosenbluth, A. W.; Rosenbluth, M. N.; Teller, A. H.; Teller, E. The journal of chemical physics 1953, 21, 1087-1092.
