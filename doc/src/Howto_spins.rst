Magnetic spins
==============

The magnetic spin simulations are enabled by the SPIN package, whose
implementation is detailed in :ref:`Tranchida <Tranchida>`.

The model represents the simulation of atomic magnetic spins coupled
to lattice vibrations. The dynamics of those magnetic spins can be used
to simulate a broad range a phenomena related to magneto-elasticity, or
or to study the influence of defects on the magnetic properties of
materials.

The magnetic spins are interacting with each others and with the
lattice via pair interactions. Typically, the magnetic exchange
interaction can be defined using the
:doc:`pair/spin/exchange <pair_spin_exchange>` command. This exchange
applies a magnetic torque to a given spin, considering the orientation
of its neighboring spins and their relative distances.
It also applies a force on the atoms as a function of the spin
orientations and their associated inter-atomic distances.

The command :doc:`fix precession/spin <fix_precession_spin>` allows to
apply a constant magnetic torque on all the spins in the system. This
torque can be an external magnetic field (Zeeman interaction), and an
uniaxial or cubic magnetic anisotropy.

A Langevin thermostat can be applied to those magnetic spins using
:doc:`fix langevin/spin <fix_langevin_spin>`. Typically, this thermostat
can be coupled to another Langevin thermostat applied to the atoms
using :doc:`fix langevin <fix_langevin>` in order to simulate
thermostatted spin-lattice systems.

The magnetic damping can also be applied
using :doc:`fix langevin/spin <fix_langevin_spin>`.
It allows to either dissipate the thermal energy of the Langevin
thermostat, or to perform a relaxation of the magnetic configuration
toward an equilibrium state.

The command :doc:`fix setforce/spin <fix_setforce>` allows to set the
components of the magnetic precession vectors (while erasing and
replacing the previously computed magnetic precession vectors on
the atom).
This command can be used to freeze the magnetic moment of certain
atoms in the simulation by zeroing their precession vector.

The command :doc:`fix nve/spin <fix_nve_spin>` can be used to
perform a symplectic integration of the combined dynamics of spins
and atomic motions.

The minimization style :doc:`min/spin <min_spin>` can be applied
to the spins to perform a minimization of the spin configuration.

Inertial spin dynamics (TSPIN)
------------------------------

.. versionadded:: TBD

TSPIN provides auxiliary inertial dynamics for equilibrium sampling of
spin-lattice systems with variable spin magnitudes.  Both the direction
and magnitude evolve under a Newtonian equation of motion.  The spin
velocities and masses are auxiliary sampling variables; their trajectories
do not represent physical Landau-Lifshitz precession or relaxation times.

TSPIN uses the existing :doc:`atom_style spin <atom_style>` and requires
magnetic interactions that provide the full Cartesian spin force, including
its longitudinal component.  See :doc:`fix nve/tspin <fix_nve_tspin>` for
the force convention and compatibility requirements.

The following commands support TSPIN simulations:

* :doc:`fix nve/tspin <fix_nve_tspin>` integrates the auxiliary spin
  and lattice equations of motion.
* :doc:`fix nvt/tspin, npt/tspin and nph/tspin <fix_nvt_tspin>` provide
  Nose-Hoover temperature and/or pressure control; *nph/tspin* does
  not thermostat the spins or lattice.
* :doc:`fix langevin/tspin <fix_langevin_tspin>` thermostats spin
  velocities when combined with *nve/tspin*.
* :doc:`fix spring/tspin <fix_spring_tspin>` provides a longitudinal
  restoring potential when needed.
* :doc:`velocity/tspin <velocity_tspin>` initializes auxiliary spin
  velocities and can assign spin masses.
* :doc:`compute temp/tspin <compute_temp_tspin>` reports the auxiliary
  spin temperature and kinetic energy.

The following fixed-lattice example assumes that the simulation box,
atoms, masses, and initial spins have been defined.  It demonstrates
commands using a harmonic modulus potential without spin-spin interactions.

.. code-block:: LAMMPS

   pair_style zero 4.0
   pair_coeff * *
   fix 1 all spring/tspin 1.0 2.2
   fix_modify 1 energy yes
   fix 2 all nve/tspin lattice frozen
   velocity/tspin all create 300.0 12345 spinmass 0.0075
   fix 3 all langevin/tspin 300.0 300.0 0.05 48279
   compute 1 all temp/tspin
   thermo_style custom step pe c_1 c_1[1]

Magnetic output
^^^^^^^^^^^^^^^

:doc:`compute spin <compute_spin>` reports averages of unit spin
directions, selected magnetic energies, and a torque-based temperature.
For variable-magnitude spins, these directional averages are not
modulus-weighted magnetic moments, and this temperature is not the
TSPIN kinetic temperature.  Use :doc:`compute property/atom
<compute_property_atom>` to access spin magnitudes and directions.

The usual *ke*, *etotal*, and *econserve* exclude auxiliary spin kinetic
energy.  See :doc:`compute temp/tspin <compute_temp_tspin>` for energy
conservation checks and the kinetic spin temperature.

Do not combine inertial and fixed-modulus spin integrators on the same
atoms.  With SPIN pair styles, this restriction also applies to disjoint
groups.  Choose spin masses and timesteps for stable integration and
converged equilibrium averages.

.. note::

   Do not :doc:`unfix <unfix>` the automatically created *TSPIN_STATE*
   fix while TSPIN styles are in use.

----------

.. _Tranchida:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
