Granular models
===============

Granular systems are composed of finite-size spherical or aspherical
particles with explicit diameters, as opposed to point particles.
This means they have an angular velocity and torque can be imparted
to them to cause them to rotate.

The various atom, pair, fix, and compute styles listed below on this
page are useful for creating granular models with such particles.

You can also define granular "surfaces" which are a collection of line
segments (2d systems) or triangles (3d systems), which act as
boundaries interacting with the particles.  This is particularly
useful for defining a complex wall or boundary geometry.
Particle/surface interactions can be specified with similar options as
the pair styles listed below for particle/particle interactions.

These doc pages describe how to define and use granular surfaces:

* :doc:`Howto granular surfaces <Howto_granular_surfaces>`
* :doc:`fix surface/global <fix_surface_global>`
* :doc:`fix surface/local <fix_surface_local>`
* :doc:`pair_style surf/granular <pair_surf_granular>`

----------

To run a simulation of particles with a granular model, you will
typically want to use the following commands:

* :doc:`atom_style sphere <atom_style>`
* :doc:`fix nve/sphere <fix_nve_sphere>`
* :doc:`fix gravity <fix_gravity>`

Aspherical granular particles can be simulated by creating clusters of spherical
particles using either the :doc:`rigid <fix_rigid>` or :doc:`BPM <Howto_bpm>`
package or by using :doc:`superellipsoids <pair_granular_superellipsoid>`.

This compute

* :doc:`compute erotate/sphere <compute_erotate_sphere>`

calculates rotational kinetic energy which can be :doc:`output with
thermodynamic info <Howto_output>`.

The compute

* :doc:`compute fabric <compute_fabric>`

calculates various versions of the fabric tensor for granular and
non-granular pair styles.

Use one of these 4 pair potentials, which compute forces and torques
between interacting pairs of particles:

* :doc:`pair_style gran/history <pair_gran>`
* :doc:`pair_style gran/no_history <pair_gran>`
* :doc:`pair_style gran/hertzian <pair_gran>`
* :doc:`pair_style granular <pair_granular>`

The first 3 are the older, original granular pair styles implemented
in LAMMPS.  The 4th :doc:`pair_style gran/no_history <pair_gran>` is
the newest and most flexible.  It's pairwise interaction models
include options similar to the first 3 pair styles, plus many
additional models.

To add your own custom granular contact model to the :doc:`pair_style
granular <pair_granular>` command, see the :doc:`Modifying granular
sub-models <Modify_gran_sub_mod>` doc page.

----------

These commands implement fix options specific to granular systems:

* :doc:`fix freeze <fix_freeze>`
* :doc:`fix pour <fix_pour>`
* :doc:`fix viscous <fix_viscous>`
* :doc:`fix wall/gran <fix_wall_gran>`
* :doc:`fix wall/gran/region <fix_wall_gran_region>`

The fix style *freeze* zeroes both the force and torque of frozen
atoms, and should be used for granular systems instead of the fix style
*setforce*\ .

To model heat conduction, one must add the temperature and heatflow
atom variables with:

* :doc:`fix property/atom <fix_property_atom>`

plus a temperature integration fix

* :doc:`fix heat/flow <fix_heat_flow>`

and a heat conduction option defined in both these commands:

* :doc:`pair_style granular <pair_granular>`
* :doc:`fix wall/gran <fix_wall_gran>`

For computational efficiency, you can eliminate needless pairwise
computations between frozen atoms by using this command:

* :doc:`neigh_modify <neigh_modify>` exclude

.. note::

   By default, for 2d systems, granular particles are still modeled
   as 3d spheres, not 2d discs (circles), meaning their moment of inertia
   will be the same as in 3d.  If you wish to model granular particles in
   2d as 2d discs, see the note on this topic on the :doc:`Howto 2d <Howto_2d>`
   doc page, where 2d simulations are discussed.
