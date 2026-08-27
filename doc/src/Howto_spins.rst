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

Inertial spin dynamics
======================

.. versionadded:: TBD

All of the above describes fixed-modulus spin dynamics, in which the
magnitude of each spin is a constant of the motion and only its
direction evolves.  The SPIN package also supports inertial spin
dynamics, in which the spin modulus itself is a dynamical degree of
freedom.  Each spin then carries a spin velocity and a spin mass, and
obeys a Newtonian second-order equation of motion instead of the
first-order Landau-Lifshitz precession equation.

It runs on the same :doc:`atom_style spin <atom_style>` and reuses the
per-atom *sp* and *fm* arrays, so existing SPIN data and restart files do not
need a new atom style.  It does not, however, accept every fixed-modulus SPIN
potential.  An inertial integrator needs all three components of the derivative
with respect to the unconstrained spin vector, including the component that
changes its modulus.

* :doc:`fix nve/tspin <fix_nve_tspin>` integrates the equations of
  motion.
* :doc:`fix nvt/tspin, fix npt/tspin and fix nph/tspin <fix_nvt_tspin>`
  do the same and add a Nose-Hoover chain on the spin velocities on top of
  the usual lattice thermostat and barostat, while
  :doc:`fix langevin/tspin <fix_langevin_tspin>` provides a Langevin
  bath for the spin degrees of freedom alone.
* :doc:`fix spring/tspin <fix_spring_tspin>` supplies a longitudinal
  potential on the spin modulus.  The modulus is a free coordinate, so
  unless the magnetic potential itself restores it, this fix or an
  equivalent one is required.
* :doc:`velocity/tspin <velocity_tspin>` initializes the spin velocities
  at a given spin temperature.
* :doc:`compute ke/tspin <compute_ke_tspin>` reports the kinetic energy
  of the spin degrees of freedom, which is not part of the
  thermodynamic keyword *ke*.

A minimal thermostatted example on a fixed lattice is:

.. code-block:: LAMMPS

   atom_style      spin
   pair_style      zero 4.0
   pair_coeff      * *

   fix             pin  all spring/tspin 1.0 2.2
   fix             1    all nve/tspin lattice frozen spinmass 0.0075
   velocity/tspin  all create 300.0 12345
   fix             bath all langevin/tspin 300.0 300.0 0.05 48279

   compute         ske all ke/tspin
   thermo_style    custom step pe f_pin c_ske

The per-atom array *fm* keeps the rad.THz units used by the SPIN package.  A
TSPIN-compatible interaction must encode the full magnetic force

.. math::

   \vec{F}^{m}_i = -\frac{\partial U}{\partial \vec{S}_i}

in that array as

.. math::

   \vec{fm}_i = \frac{|\vec{S}_i|}{\hbar}\vec{F}^{m}_i.

The inertial integrator multiplies *fm* by
:math:`\hbar/|\vec{S}_i|` to recover the force.  This is a stronger contract
than the one needed by fixed-modulus Landau-Lifshitz dynamics, where any
component of *fm* parallel to the spin disappears from the cross product.  The
existing SPIN pair styles only guarantee the resulting torque and are therefore
rejected by the TSPIN integrators.  A variable-moment potential that explicitly
uses the full-gradient encoding, such as *pair_style deepspin* of the DeePMD-kit
package, can be used without changing its force output.

The inertial and the fixed-modulus styles describe different physics
and must not be combined on the same group of atoms.  Because the
resulting spin dynamics is second order in time, its characteristic
frequencies scale as the inverse square root of the spin mass and are
not the Landau-Lifshitz precession frequencies.  Static thermodynamic
averages should be verified to be independent of the spin mass.

----------

All the computed magnetic properties can be output by two main
commands. The first one is :doc:`compute spin <compute_spin>`, that
enables to evaluate magnetic averaged quantities, such as the total
magnetization of the system along x, y, or z, the spin temperature, or
the magnetic energy. The second command
is :doc:`compute property/atom <compute_property_atom>`.
It enables to output all the per atom magnetic quantities. Typically,
the orientation of a given magnetic spin, or the magnetic force
acting on this spin.

----------

.. _Tranchida:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
