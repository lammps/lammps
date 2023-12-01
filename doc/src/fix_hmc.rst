fix hmc command
==============================

Accelerator Variants: *None*

Syntax
""""""
.. code-block:: LAMMPS

   fix ID group-ID hmc N seed temp integrator

* ID = user-assigned name for the fix
* group-ID = ID of the group of atoms to apply the fix to
* hmc = style name of this fix command
* N = number of integrator timesteps between HMC evaluations
* seed = random number seed
* temp = temperature for assigning velocities
* integrator = integrator fix: flexible (for nve) or rigid (for rigid/small)

  .. parsed-literal::

     keyword = *mom*
       *mom* value = *no* or *yes*
Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all hmc 10 123 500 flexible
   fix hmc_water water 100 123 298.15 rigid

Description
"""""""""""
Perform the the Hybrid/Hamiltonian Monte Carlo (HMC) algorithm in line with the following order of steps:

The new particle positions and velocities are calculated by invoking the velocity form of the Stoermer-Verlet time integration algorithm (velocity-Verlet).
Before these position and velocity changes are performed, the proposed change in the Hamiltonian
:math:`\Delta{H}`
is calculated following the equation:

.. math::

   \Delta{H} = H(q′,p′) - H(q,p)


This new proposed configuration is then accepted/rejected with the probability:

.. math::

   p^{acc} = min(1,e^{\frac{-\Delta{H}}{k_B T}})

If the configuration is accepted, the positions and velocities are updated.

If the configuration is rejected, particle momenta are randomly resampled from a normal distribution:

.. math::

   p_{x,y,z} = \textbf{N}(0,1) \sqrt{\frac{k_B T}{2 m^2}}

the simulation is then continued, where a new MD step is proposed, and the procedure is repeated.

Related commands
""""""""""""""""

fix tfmc, fix gcmc, fix nve

Default
"""""""

none

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

