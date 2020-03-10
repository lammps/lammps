Drude induced dipoles
=====================

The thermalized Drude model represents induced dipoles by a pair of
charges (the core atom and the Drude particle) connected by a harmonic
spring.  See the :doc:`Howto polarizable <Howto_polarizable>` doc page
for a discussion of all the polarizable models available in LAMMPS.

The Drude model has a number of features aimed at its use in
molecular systems (:ref:`Lamoureux and Roux <howto-Lamoureux>`):

* Thermostatting of the additional degrees of freedom associated with the
  induced dipoles at very low temperature, in terms of the reduced
  coordinates of the Drude particles with respect to their cores. This
  makes the trajectory close to that of relaxed induced dipoles.
* Consistent definition of 1-2 to 1-4 neighbors. A core-Drude particle
  pair represents a single (polarizable) atom, so the special screening
  factors in a covalent structure should be the same for the core and
  the Drude particle.  Drude particles have to inherit the 1-2, 1-3, 1-4
  special neighbor relations from their respective cores.
* Stabilization of the interactions between induced dipoles. Drude
  dipoles on covalently bonded atoms interact too strongly due to the
  short distances, so an atom may capture the Drude particle of a
  neighbor, or the induced dipoles within the same molecule may align
  too much. To avoid this, damping at short range can be done by Thole
  functions (for which there are physical grounds). This Thole damping
  is applied to the point charges composing the induced dipole (the
  charge of the Drude particle and the opposite charge on the core, not
  to the total charge of the core atom).

A detailed tutorial covering the usage of Drude induced dipoles in
LAMMPS is on the :doc:`Howto drude2e <Howto_drude2>` doc page.

As with the core-shell model, the cores and Drude particles should
appear in the data file as standard atoms. The same holds for the
springs between them, which are described by standard harmonic bonds.
The nature of the atoms (core, Drude particle or non-polarizable) is
specified via the :doc:`fix drude <fix_drude>` command.  The special
list of neighbors is automatically refactored to account for the
equivalence of core and Drude particles as regards special 1-2 to 1-4
screening. It may be necessary to use the *extra/special/per/atom*
keyword of the :doc:`read_data <read_data>` command. If using :doc:`fix shake <fix_shake>`, make sure no Drude particle is in this fix
group.

There are two ways to thermostat the Drude particles at a low
temperature: use either :doc:`fix langevin/drude <fix_langevin_drude>`
for a Langevin thermostat, or :doc:`fix drude/transform/\* <fix_drude_transform>` for a Nose-Hoover
thermostat. The former requires use of the command :doc:`comm_modify vel yes <comm_modify>`. The latter requires two separate integration
fixes like *nvt* or *npt*\ . The correct temperatures of the reduced
degrees of freedom can be calculated using the :doc:`compute temp/drude <compute_temp_drude>`. This requires also to use the
command *comm\_modify vel yes*.

Short-range damping of the induced dipole interactions can be achieved
using Thole functions through the :doc:`pair style thole <pair_thole>` in :doc:`pair_style hybrid/overlay <pair_hybrid>`
with a Coulomb pair style. It may be useful to use *coul/long/cs* or
similar from the CORESHELL package if the core and Drude particle come
too close, which can cause numerical issues.

----------

.. _howto-Lamoureux:

**(Lamoureux and Roux)** G. Lamoureux, B. Roux, J. Chem. Phys 119, 3025 (2003)
