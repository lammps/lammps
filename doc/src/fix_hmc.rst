fix hmc command
==============================

Accelerator Variants: *None*

Syntax
""""""
::

   fix ID group-ID hmc

* ID = user-assigned name for the fix
* group-ID = ID of the group of atoms to apply the fix to
* hmc = style name of this fix command

Examples
""""""""

::

   fix 1 all hmc
   fix hmc_water water hmc

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

   p = min(1,e^{-\Delta{H}})

If the configuration is accepted, the positions and velocities are updated.

If the configuration is rejected, the previous positions and velocities are kept.


