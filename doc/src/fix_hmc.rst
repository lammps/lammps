fix hmc command
==============================

Accelerator Variants: *None*

Syntax
""""""
.. code-block:: LAMMPS

   fix ID group-ID hmc N seed temp integrator

* ID, group-ID are documented in :doc:`fix <fix>` command
* hmc = style name of this fix command
* N = invoke this fix every N steps
* seed = random # seed (positive integer)
* temp = temperature for assigning velocities
* integrator = integrator fix: flexible (for nve) or rigid (for rigid/small)

  .. parsed-literal::

     keyword = *mom*
       *mom* value = *no* or *yes*
Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all hmc 10 123 500 flexible
   fix hmc_water water hmc 100 123 298.15 rigid

Description
"""""""""""
This fix performs the the Hybrid/Hamiltonian Monte Carlo (HMC) algorithm in line with the following order of steps:

The new particle positions and velocities are calculated by invoking the velocity-Verlet time integration algorithm.
Before these position and velocity changes are performed, the proposed change in the Hamiltonian,
:math:`\Delta{H}`
is calculated following the equation:

.. math::

   \Delta{H} = H(q′,p′) - H(q,p)


This new proposed configuration is then accepted/rejected according to the Metropolis criterion with probability:

.. math::

   p^{acc} = min(1,e^{\frac{-\Delta{H}}{k_B T}})

Upon acceptance, the new particle configuration positions and velocities (momenta) are updated.

Upon rejection, the old particle configuration is kept, and particle momenta are randomly resampled from a normal distribution:

.. math::

   p_{x,y,z} = \textbf{N}(0,1) \sqrt{\frac{k_B T}{2 m^2}}

The algorithm then continues, proposing a new MD step until a configuration is accepted.

Restrictions
""""""""""""

This fix is part of the MC package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix gcmc <fix_gcmc>`, :doc:`fix tfmc <fix_tfmc>`

Default
"""""""

The option default is mom = yes

----------

**(Watkins)** Watkins and Jorgensen, J Phys Chem A, 105, 4118-4125 (2001).

**(Betancourt)** Betancourt, M. A Conceptual Introduction to Hamiltonian Monte Carlo, 2018.

**(Duane)** Duane, S.; Kennedy, A. D.; Pendleton, B. J.; Roweth, D. Hybrid Monte Carlo. Physics Letters B 1987, 195 (2), 216–222. https://doi.org/10.1016/0370-2693(87)91197-X.

**(Metropolis)** Metropolis, N.; Rosenbluth, A. W.; Rosenbluth, M. N.; Teller, A. H.; Teller, E. The journal of chemical physics
1953, 21, 1087–1092.

LAMMPS Developers Issue 565: [Brief description of the issue] GitHub issue, https://github.com/lammps/
lammps/issues/565.

LAMMPS Development Team LAMMPS Documentation: Modify Requirements https://docs.lammps.org/
Modify_requirements.html.

