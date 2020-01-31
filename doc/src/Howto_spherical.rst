Finite-size spherical and aspherical particles
==============================================

Typical MD models treat atoms or particles as point masses.  Sometimes
it is desirable to have a model with finite-size particles such as
spheroids or ellipsoids or generalized aspherical bodies.  The
difference is that such particles have a moment of inertia, rotational
energy, and angular momentum.  Rotation is induced by torque coming
from interactions with other particles.

LAMMPS has several options for running simulations with these kinds of
particles.  The following aspects are discussed in turn:

* atom styles
* pair potentials
* time integration
* computes, thermodynamics, and dump output
* rigid bodies composed of finite-size particles

Example input scripts for these kinds of models are in the body,
colloid, dipole, ellipse, line, peri, pour, and tri directories of the
:doc:`examples directory <Examples>` in the LAMMPS distribution.

Atom styles
-----------

There are several :doc:`atom styles <atom_style>` that allow for
definition of finite-size particles: sphere, dipole, ellipsoid, line,
tri, peri, and body.

The sphere style defines particles that are spheroids and each
particle can have a unique diameter and mass (or density).  These
particles store an angular velocity (omega) and can be acted upon by
torque.  The "set" command can be used to modify the diameter and mass
of individual particles, after then are created.

The dipole style does not actually define finite-size particles, but
is often used in conjunction with spherical particles, via a command
like


.. parsed-literal::

   atom_style hybrid sphere dipole

This is because when dipoles interact with each other, they induce
torques, and a particle must be finite-size (i.e. have a moment of
inertia) in order to respond and rotate.  See the :doc:`atom_style dipole <atom_style>` command for details.  The "set" command can be
used to modify the orientation and length of the dipole moment of
individual particles, after then are created.

The ellipsoid style defines particles that are ellipsoids and thus can
be aspherical.  Each particle has a shape, specified by 3 diameters,
and mass (or density).  These particles store an angular momentum and
their orientation (quaternion), and can be acted upon by torque.  They
do not store an angular velocity (omega), which can be in a different
direction than angular momentum, rather they compute it as needed.
The "set" command can be used to modify the diameter, orientation, and
mass of individual particles, after then are created.  It also has a
brief explanation of what quaternions are.

The line style defines line segment particles with two end points and
a mass (or density).  They can be used in 2d simulations, and they can
be joined together to form rigid bodies which represent arbitrary
polygons.

The tri style defines triangular particles with three corner points
and a mass (or density).  They can be used in 3d simulations, and they
can be joined together to form rigid bodies which represent arbitrary
particles with a triangulated surface.

The peri style is used with :doc:`Peridynamic models <pair_peri>` and
defines particles as having a volume, that is used internally in the
:doc:`pair_style peri <pair_peri>` potentials.

The body style allows for definition of particles which can represent
complex entities, such as surface meshes of discrete points,
collections of sub-particles, deformable objects, etc.  The body style
is discussed in more detail on the :doc:`Howto body <Howto_body>` doc
page.

Note that if one of these atom styles is used (or multiple styles via
the :doc:`atom_style hybrid <atom_style>` command), not all particles in
the system are required to be finite-size or aspherical.

For example, in the ellipsoid style, if the 3 shape parameters are set
to the same value, the particle will be a sphere rather than an
ellipsoid.  If the 3 shape parameters are all set to 0.0 or if the
diameter is set to 0.0, it will be a point particle.  In the line or
tri style, if the lineflag or triflag is specified as 0, then it
will be a point particle.

Some of the pair styles used to compute pairwise interactions between
finite-size particles also compute the correct interaction with point
particles as well, e.g. the interaction between a point particle and a
finite-size particle or between two point particles.  If necessary,
:doc:`pair_style hybrid <pair_hybrid>` can be used to insure the correct
interactions are computed for the appropriate style of interactions.
Likewise, using groups to partition particles (ellipsoids versus
spheres versus point particles) will allow you to use the appropriate
time integrators and temperature computations for each class of
particles.  See the doc pages for various commands for details.

Also note that for :doc:`2d simulations <dimension>`, atom styles sphere
and ellipsoid still use 3d particles, rather than as circular disks or
ellipses.  This means they have the same moment of inertia as the 3d
object.  When temperature is computed, the correct degrees of freedom
are used for rotation in a 2d versus 3d system.

Pair potentials
---------------

When a system with finite-size particles is defined, the particles
will only rotate and experience torque if the force field computes
such interactions.  These are the various :doc:`pair styles <pair_style>` that generate torque:

* :doc:`pair_style gran/history <pair_gran>`
* :doc:`pair_style gran/hertzian <pair_gran>`
* :doc:`pair_style gran/no\_history <pair_gran>`
* :doc:`pair_style dipole/cut <pair_dipole>`
* :doc:`pair_style gayberne <pair_gayberne>`
* :doc:`pair_style resquared <pair_resquared>`
* :doc:`pair_style brownian <pair_brownian>`
* :doc:`pair_style lubricate <pair_lubricate>`
* :doc:`pair_style line/lj <pair_line_lj>`
* :doc:`pair_style tri/lj <pair_tri_lj>`
* :doc:`pair_style body/nparticle <pair_body_nparticle>`

The granular pair styles are used with spherical particles.  The
dipole pair style is used with the dipole atom style, which could be
applied to spherical or ellipsoidal particles.  The GayBerne and
REsquared potentials require ellipsoidal particles, though they will
also work if the 3 shape parameters are the same (a sphere).  The
Brownian and lubrication potentials are used with spherical particles.
The line, tri, and body potentials are used with line segment,
triangular, and body particles respectively.

Time integration
----------------

There are several fixes that perform time integration on finite-size
spherical particles, meaning the integrators update the rotational
orientation and angular velocity or angular momentum of the particles:

* :doc:`fix nve/sphere <fix_nve_sphere>`
* :doc:`fix nvt/sphere <fix_nvt_sphere>`
* :doc:`fix npt/sphere <fix_npt_sphere>`

Likewise, there are 3 fixes that perform time integration on
ellipsoidal particles:

* :doc:`fix nve/asphere <fix_nve_asphere>`
* :doc:`fix nvt/asphere <fix_nvt_asphere>`
* :doc:`fix npt/asphere <fix_npt_asphere>`

The advantage of these fixes is that those which thermostat the
particles include the rotational degrees of freedom in the temperature
calculation and thermostatting.  The :doc:`fix langevin <fix_langevin>`
command can also be used with its *omgea* or *angmom* options to
thermostat the rotational degrees of freedom for spherical or
ellipsoidal particles.  Other thermostatting fixes only operate on the
translational kinetic energy of finite-size particles.

These fixes perform constant NVE time integration on line segment,
triangular, and body particles:

* :doc:`fix nve/line <fix_nve_line>`
* :doc:`fix nve/tri <fix_nve_tri>`
* :doc:`fix nve/body <fix_nve_body>`

Note that for mixtures of point and finite-size particles, these
integration fixes can only be used with :doc:`groups <group>` which
contain finite-size particles.

Computes, thermodynamics, and dump output
-----------------------------------------

There are several computes that calculate the temperature or
rotational energy of spherical or ellipsoidal particles:

* :doc:`compute temp/sphere <compute_temp_sphere>`
* :doc:`compute temp/asphere <compute_temp_asphere>`
* :doc:`compute erotate/sphere <compute_erotate_sphere>`
* :doc:`compute erotate/asphere <compute_erotate_asphere>`

These include rotational degrees of freedom in their computation.  If
you wish the thermodynamic output of temperature or pressure to use
one of these computes (e.g. for a system entirely composed of
finite-size particles), then the compute can be defined and the
:doc:`thermo_modify <thermo_modify>` command used.  Note that by default
thermodynamic quantities will be calculated with a temperature that
only includes translational degrees of freedom.  See the
:doc:`thermo_style <thermo_style>` command for details.

These commands can be used to output various attributes of finite-size
particles:

* :doc:`dump custom <dump>`
* :doc:`compute property/atom <compute_property_atom>`
* :doc:`dump local <dump>`
* :doc:`compute body/local <compute_body_local>`

Attributes include the dipole moment, the angular velocity, the
angular momentum, the quaternion, the torque, the end-point and
corner-point coordinates (for line and tri particles), and
sub-particle attributes of body particles.

Rigid bodies composed of finite-size particles
----------------------------------------------

The :doc:`fix rigid <fix_rigid>` command treats a collection of
particles as a rigid body, computes its inertia tensor, sums the total
force and torque on the rigid body each timestep due to forces on its
constituent particles, and integrates the motion of the rigid body.

If any of the constituent particles of a rigid body are finite-size
particles (spheres or ellipsoids or line segments or triangles), then
their contribution to the inertia tensor of the body is different than
if they were point particles.  This means the rotational dynamics of
the rigid body will be different.  Thus a model of a dimer is
different if the dimer consists of two point masses versus two
spheroids, even if the two particles have the same mass.  Finite-size
particles that experience torque due to their interaction with other
particles will also impart that torque to a rigid body they are part
of.

See the "fix rigid" command for example of complex rigid-body models
it is possible to define in LAMMPS.

Note that the :doc:`fix shake <fix_shake>` command can also be used to
treat 2, 3, or 4 particles as a rigid body, but it always assumes the
particles are point masses.

Also note that body particles cannot be modeled with the :doc:`fix rigid <fix_rigid>` command.  Body particles are treated by LAMMPS
as single particles, though they can store internal state, such as a
list of sub-particles.  Individual body particles are typically treated
as rigid bodies, and their motion integrated with a command like :doc:`fix nve/body <fix_nve_body>`.  Interactions between pairs of body
particles are computed via a command like :doc:`pair_style body/nparticle <pair_body_nparticle>`.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
