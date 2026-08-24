.. index:: fix spring/tspin

fix spring/tspin command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID spring/tspin K S0

* ID, group-ID are documented in :doc:`fix <fix>` command
* spring/tspin = style name of this fix command
* K = spring constant of the modulus restoring potential (energy units)
* S0 = equilibrium value of the spin modulus (adim)

Examples
""""""""

.. code-block:: LAMMPS

   fix pin all spring/tspin 1.0 2.2
   fix landau all spring/tspin 1000.0 1.72

Description
"""""""""""

.. versionadded:: TBD

Apply a harmonic restoring potential to the modulus of each magnetic
spin in the group,

.. math::

   U_L = \frac{1}{2} K \left( \left| \vec{S}_i \right| - S_0 \right)^2

which adds the radial magnetic force

.. math::

   \vec{F}^{m}_{i} \rightarrow \vec{F}^{m}_{i}
   - K \left( \left| \vec{S}_i \right| - S_0 \right) \hat{s}_i

This fix supplies the longitudinal energy scale that inertial spin
dynamics requires.  The SPIN pair styles provide a Landau-Lifshitz
effective field, whose component parallel to :math:`\vec{S}_i` does no
work on a fixed-modulus spin but acts as a radial driving force once
the modulus is a dynamical degree of freedom.  A simulation using
:doc:`fix nve/tspin <fix_nve_tspin>`, :doc:`fix nvt/tspin
<fix_nvt_tspin>` or :doc:`fix langevin/tspin <fix_langevin_tspin>`
therefore has to include this fix, or some other longitudinal
potential, otherwise the spin modulus is unbounded.

In thermal equilibrium at spin temperature :math:`T_s` and for a stiff
spring, the modulus fluctuates about :math:`S_0` with variance
:math:`k_B T_s / K`, which is a convenient check that the longitudinal
degree of freedom is correctly thermostatted.

This fix has no effect on the direction of the spins, and does not
affect :doc:`fix nve/spin <fix_nve_spin>` simulations, in which the
modulus is constant.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy stored in the modulus springs to the global
potential energy of the system as part of :doc:`thermodynamic output
<thermo_style>`.

This fix computes a global scalar, the energy stored in the modulus
springs, which can be accessed by various :doc:`output commands
<Howto_output>`.  The scalar value is "extensive".

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

Restrictions
""""""""""""

The *spring/tspin* fix is part of the SPIN package.  This style is only
enabled if LAMMPS was built with this package.  See the :doc:`Build
package <Build_package>` page for more info.

This fix requires :doc:`atom_style tspin <atom_style>`.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`,
:doc:`fix precession/spin <fix_precession_spin>`

Default
""""""""

none
