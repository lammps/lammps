Granular models
===============

Granular system are composed of spherical particles with a diameter,
as opposed to point particles.  This means they have an angular
velocity and torque can be imparted to them to cause them to rotate.

To run a simulation of a granular model, you will want to use
the following commands:

* :doc:`atom\_style sphere <atom_style>`
* :doc:`fix nve/sphere <fix_nve_sphere>`
* :doc:`fix gravity <fix_gravity>`

This compute

* :doc:`compute erotate/sphere <compute_erotate_sphere>`

calculates rotational kinetic energy which can be :doc:`output with thermodynamic info <Howto_output>`.

Use one of these 3 pair potentials, which compute forces and torques
between interacting pairs of particles:

* :doc:`pair\_style <pair_style>` gran/history
* :doc:`pair\_style <pair_style>` gran/no\_history
* :doc:`pair\_style <pair_style>` gran/hertzian

These commands implement fix options specific to granular systems:

* :doc:`fix freeze <fix_freeze>`
* :doc:`fix pour <fix_pour>`
* :doc:`fix viscous <fix_viscous>`
* :doc:`fix wall/gran <fix_wall_gran>`

The fix style *freeze* zeroes both the force and torque of frozen
atoms, and should be used for granular system instead of the fix style
*setforce*\ .

For computational efficiency, you can eliminate needless pairwise
computations between frozen atoms by using this command:

* :doc:`neigh\_modify <neigh_modify>` exclude

.. note::

   By default, for 2d systems, granular particles are still modeled
   as 3d spheres, not 2d discs (circles), meaning their moment of inertia
   will be the same as in 3d.  If you wish to model granular particles in
   2d as 2d discs, see the note on this topic on the :doc:`Howto 2d <Howto_2d>`
   doc page, where 2d simulations are discussed.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
