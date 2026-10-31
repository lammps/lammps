.. index:: fix nve/spin

fix nve/spin command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/spin keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/spin = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *lattice* or *selfconsistent*

  .. parsed-literal::

       *lattice* value = *moving* or *frozen*
         moving = integrate both spin and atomic degrees of freedom
         frozen = integrate spins on a fixed lattice
       *selfconsistent* value = N
         N = number of iterations for the field of each single spin update (>= 0)

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all nve/spin lattice moving
   fix 1 all nve/spin lattice frozen
   fix 1 all nve/spin lattice frozen selfconsistent 2

Description
"""""""""""

Perform a symplectic integration for the spin or spin-lattice system.

The *lattice* keyword defines if the spins are integrated on a lattice
of fixed atoms (lattice = frozen), or if atoms are moving
(lattice = moving).
The first case corresponds to a spin dynamics calculation, and
the second to a spin-lattice calculation.
By default a spin-lattice integration is performed (lattice = moving).

The *nve/spin* fix applies a Suzuki-Trotter decomposition to
the equations of motion of the spin lattice system, following the scheme:

.. image:: JPG/fix_integration_spin_stdecomposition.jpg
   :align: center

according to the implementation reported in :ref:`(Omelyan) <Omelyan1>`.

A sectoring method enables this scheme for parallel calculations.
The implementation of this sectoring algorithm is reported
in :ref:`(Tranchida) <Tranchida1>`.

.. versionadded:: TBD

Each spin is rotated about its precession field, which is computed just
before its update and held constant during the rotation.  If this field
depends on the spin itself, as for magnetic anisotropy (:doc:`fix
precession/spin <fix_precession_spin>` with *anisotropy* or *cubic*),
:doc:`biquadratic exchange <pair_spin_exchange>`, or :doc:`pair spin/neel
<pair_spin_neel>`, this update is no longer time-symmetric, and the
energy drifts steadily, at a rate proportional to the timestep.  With
the *selfconsistent* keyword, the field of each spin update is instead
evaluated at the midpoint between the old and the new spin direction,
found by *N* iterations, each of which computes the field on that spin
once more.  This removes the drift.  In tests with 1024 spins of bcc
iron at 300 K on a frozen lattice, with timesteps of 0.1 to 0.4 fs,
*N* = 2 reduced the energy drift by a factor of 10\ :sup:`5` or more
for anisotropy and biquadratic exchange, and of 10\ :sup:`3` or more
for *spin/neel*.  The energy was then conserved to round-off for
uniaxial anisotropy and, with *N* = 4, for biquadratic exchange.  Each
iteration adds one field evaluation per spin update: with exchange and
anisotropy only, *N* = 2 made the run about 2.4 times slower, and with
an EAM potential and a moving lattice about 1.15 times slower.  For
interactions that do not depend on
the spin itself (e.g. Heisenberg exchange, DMI, Zeeman), the results do
not change.  The default *N* = 0 gives the original update.

----------

Restrictions
""""""""""""

This fix style can only be used if LAMMPS was built with the SPIN
package.  See the :doc:`Build package <Build_package>` page for more
info.

The *selfconsistent* keyword cannot be used with :doc:`fix langevin/spin
<fix_langevin_spin>` at a non-zero temperature, since the random field
would be drawn again at each iteration.

To use the spin algorithm, it is necessary to define a map with
the atom_modify command. Typically, by adding the command:

.. code-block:: LAMMPS

   atom_modify map array

before you create the simulation box. Note that the keyword "hash"
instead of "array" is also valid.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`fix nve <fix_nve>`

Default
"""""""

The option defaults are lattice = moving and selfconsistent = 0.

----------

.. _Omelyan1:

**(Omelyan)** Omelyan, Mryglod, and Folk. Phys. Rev. Lett.
86(5), 898. (2001).

.. _Tranchida1:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
