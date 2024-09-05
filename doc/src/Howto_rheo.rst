Reproducing hydrodynamics and elastic objects (RHEO)
====================================================

The RHEO package is a hybrid implementation of smoothed particle
hydrodynamics (SPH) for fluid flow, which can couple to the :doc:`BPM package
<Howto_bpm>` to model solid elements. RHEO combines these methods to enable
mesh-free modeling of multi-phase material systems. Its SPH solver supports
many advanced options including reproducing kernels, particle shifting, free
surface identification, and solid surface reconstruction. To model fluid-solid
systems, the status of particles can dynamically change between a fluid and
solid state, e.g. during melting/solidification, which determines how they
interact and their physical behavior. The package is designed with modularity
in mind, so one can easily turn various features on/off, adjust physical
details of the system, or develop new capabilities. For instance, the numerics
associated with calculating gradients, reproducing kernels, etc. are separated
into distinct classes to simplify the development of new integration schemes
which can call these calculations. Additional numerical details can be found in
:ref:`(Clemmer) <howto_rheo_clemmer>`. Example movies illustrating some of these
capabilities are found at https://www.lammps.org/movies.html#rheopackage.

Note, if you simply want to run a traditional SPH simulation, the :ref:`SPH package
<PKG-SPH>` package is likely better suited for your application. It has fewer advanced
features and therefore benefits from improved performance. The :ref:`MACHDYN
<PKG-MACHDYN>` package for solids may also be relevant for fluid-solid problems.

----------

At the core of the package is :doc:`fix rheo <fix_rheo>` which integrates
particle trajectories and controls many optional features (e.g. the use
of reproducing kernels). In conjunction to fix rheo, one must specify an
instance of :doc:`fix rheo/pressure <fix_rheo_pressure>` and
:doc:`fix rheo/viscosity <fix_rheo_viscosity>` to define a pressure equation
of state and viscosity model, respectively. Optionally, one can model
a heat equation with :doc:`fix rheo/thermal <fix_rheo_thermal>`, which also
allows the user to specify equations for a particle's thermal conductivity,
specific heat, latent heat, and melting temperature. The ordering of these
fixes in an an input script matters. Fix rheo must be defined prior to all
other RHEO fixes.

Typically, RHEO requires atom style rheo. In addition to typical atom
properties like positions and forces, particles store a local density,
viscosity, pressure, and status. If thermal evolution is modeled, one must
use atom style rheo/thermal which also includes a local energy, temperature, and
conductivity. Note that the temperature is always derived from the energy.
This implies the *temperature* attribute of :doc:`the set command <set>` does not
affect particles. Instead, one should use the *sph/e* attribute.

The status variable uses bit-masking to track various properties of a particle
such as its current state of matter (fluid or solid) and its location relative
to a surface. Some of these properties (and others) can be accessed using
:doc:`compute rheo/property/atom <compute_rheo_property_atom>`. The *status*
attribute in :doc:`the set command <set>` only allows control over the first bit
which sets the state of matter, 0 is fluid and 1 is solid.

Fluid interactions, including pressure forces, viscous forces, and heat exchange,
are calculated using :doc:`pair rheo <pair_rheo>`. Unlike typical pair styles,
pair rheo ignores the :doc:`special bond <special_bonds>` settings. Instead,
it determines whether to calculate forces based on the status of particles: e.g.,
hydrodynamic forces are only calculated if a fluid particle is involved.

----------

To model elastic objects, there are currently two mechanisms in RHEO, one designed
for bulk solid bodies and the other for thin shells. Both mechanisms rely on
introducing bonded forces between particles and therefore require a hybrid of atom
style bond and rheo (or rheo/thermal).

To create an elastic solid body, one has to (a) change the status of constituent
particles to solid (e.g. with the :doc:`set <set>` command), (b) create bpm
bonds between the particles (see the :doc:`bpm howto <Howto_bpm>` page for
more details), and (c) use :doc:`pair rheo/solid <pair_rheo_solid>` to
apply repulsive contact forces between distinct solid bodies. Akin to pair rheo,
pair rheo/solid considers a particles fluid/solid phase to determine whether to
apply forces. However, unlike pair rheo, pair rheo/solid does obey special bond
settings such that contact forces do not have to be calculated between two bonded
solid particles in the same elastic body.

In systems with thermal evolution, fix rheo/thermal can optionally set a
melting/solidification temperature allowing particles to dynamically swap their
state between fluid and solid when the temperature exceeds or drops below the
critical temperature, respectively. Using the *react* option, one can specify a maximum
bond length and a bond type. Then, when solidifying, particles will search their
local neighbors and automatically create bonds with any neighboring solid particles
in range. For BPM bond styles, bonds will then use the immediate position of the two
particles to calculate a reference state. When melting, particles will delete any
bonds of the specified type when reverting to a fluid state. Special bonds are updated
as bonds are created/broken.

The other option for elastic objects is an elastic shell that is nominally much
thinner than a particle diameter, e.g. a oxide skin which gradually forms over time
on the surface of a fluid. Currently, this is implemented using
:doc:`fix rheo/oxidation <fix_rheo_oxidation>` and bond style
:doc:`rheo/shell <bond_rheo_shell>`. Essentially, fix rheo/oxidation creates candidate
bonds of a specified type between surface fluid particles within a specified distance.
a newly created rheo/shell bond will then start a timer. While the timer is counting
down, the bond will delete itself if particles move too far apart or move away from the
surface. However, if the timer reaches a user-defined threshold, then the bond will
activate and apply additional forces to the fluid particles. Bond style rheo/shell
then operates very similarly to a BPM bond style, storing a reference length and
breaking if stretched too far. Unlike the above method, this option does not remove
the underlying fluid interactions (although particle shifting is turned off) and does
not modify special bond settings of particles.

While these two options are not expected to be appropriate for every system,
either framework can be modified to create more suitable models (e.g. by changing the
criteria for creating/deleting a bond or altering force calculations).

----------

.. _howto_rheo_clemmer:

**(Clemmer)** Clemmer, Pierce, O'Connor, Nevins, Jones, Lechman, Tencer, Appl. Math. Model., 130, 310-326 (2024).
