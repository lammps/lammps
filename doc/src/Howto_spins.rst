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
torque can be an external magnetic field (Zeeman interaction), or an
uniaxial magnetic anisotropy.

A Langevin thermostat can be applied to those magnetic spins using
:doc:`fix langevin/spin <fix_langevin_spin>`. Typically, this thermostat
can be coupled to another Langevin thermostat applied to the atoms
using :doc:`fix langevin <fix_langevin>` in order to simulate
thermostatted spin-lattice systems.

The magnetic Gilbert damping can also be applied using :doc:`fix langevin/spin <fix_langevin_spin>`. It allows to either dissipate
the thermal energy of the Langevin thermostat, or to perform a
relaxation of the magnetic configuration toward an equilibrium state.

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

All the computed magnetic properties can be output by two main
commands. The first one is :doc:`compute spin <compute_spin>`, that
enables to evaluate magnetic averaged quantities, such as the total
magnetization of the system along x, y, or z, the spin temperature, or
the magnetic energy. The second command is :doc:`compute property/atom <compute_property_atom>`. It enables to output all the
per atom magnetic quantities. Typically, the orientation of a given
magnetic spin, or the magnetic force acting on this spin.


----------


.. _Tranchida:



**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
