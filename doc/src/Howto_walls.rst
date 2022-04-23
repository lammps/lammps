Walls
=====

Walls in an MD simulation are typically used to bound particle motion,
i.e. to serve as a boundary condition.

Walls in LAMMPS can be of rough (made of particles) or idealized
surfaces.  Ideal walls can be smooth, generating forces only in the
normal direction, or frictional, generating forces also in the
tangential direction.

Rough walls, built of particles, can be created in various ways.  The
particles themselves can be generated like any other particle, via the
:doc:`lattice <lattice>` and :doc:`create_atoms <create_atoms>` commands,
or read in via the :doc:`read_data <read_data>` command.

Their motion can be constrained by many different commands, so that
they do not move at all, move together as a group at constant velocity
or in response to a net force acting on them, move in a prescribed
fashion (e.g. rotate around a point), etc.  Note that if a time
integration fix like :doc:`fix nve <fix_nve>` or :doc:`fix nvt <fix_nh>`
is not used with the group that contains wall particles, their
positions and velocities will not be updated.

* :doc:`fix aveforce <fix_aveforce>` - set force on particles to average value, so they move together
* :doc:`fix setforce <fix_setforce>` - set force on particles to a value, e.g. 0.0
* :doc:`fix freeze <fix_freeze>` - freeze particles for use as granular walls
* :doc:`fix nve/noforce <fix_nve_noforce>` - advect particles by their velocity, but without force
* :doc:`fix move <fix_move>` - prescribe motion of particles by a linear velocity, oscillation, rotation, variable

The :doc:`fix move <fix_move>` command offers the most generality, since
the motion of individual particles can be specified with
:doc:`variable <variable>` formula which depends on time and/or the
particle position.

For rough walls, it may be useful to turn off pairwise interactions
between wall particles via the :doc:`neigh_modify exclude <neigh_modify>` command.

Rough walls can also be created by specifying frozen particles that do
not move and do not interact with mobile particles, and then tethering
other particles to the fixed particles, via a :doc:`bond <bond_style>`.
The bonded particles do interact with other mobile particles.

Idealized walls can be specified via several fix commands.  :doc:`Fix wall/gran <fix_wall_gran>` creates frictional walls for use with
granular particles; all the other commands create smooth walls.

* :doc:`fix wall/reflect <fix_wall_reflect>` - reflective flat walls
* :doc:`fix wall/lj93 <fix_wall>` - flat walls, with Lennard-Jones 9/3 potential
* :doc:`fix wall/lj126 <fix_wall>` - flat walls, with Lennard-Jones 12/6 potential
* :doc:`fix wall/colloid <fix_wall>` - flat walls, with :doc:`pair_style colloid <pair_colloid>` potential
* :doc:`fix wall/harmonic <fix_wall>` - flat walls, with repulsive harmonic spring potential
* :doc:`fix wall/morse <fix_wall>` - flat walls, with Morse potential
* :doc:`fix wall/region <fix_wall_region>` - use region surface as wall
* :doc:`fix wall/gran <fix_wall_gran>` - flat or curved walls with :doc:`pair_style granular <pair_gran>` potential

The *lj93*, *lj126*, *colloid*, *harmonic*, and *morse* styles all
allow the flat walls to move with a constant velocity, or oscillate in
time.  The :doc:`fix wall/region <fix_wall_region>` command offers the
most generality, since the region surface is treated as a wall, and
the geometry of the region can be a simple primitive volume (e.g. a
sphere, or cube, or plane), or a complex volume made from the union
and intersection of primitive volumes.  :doc:`Regions <region>` can also
specify a volume "interior" or "exterior" to the specified primitive
shape or *union* or *intersection*\ .  :doc:`Regions <region>` can also be
"dynamic" meaning they move with constant velocity, oscillate, or
rotate.

The only frictional idealized walls currently in LAMMPS are flat or
curved surfaces specified by the :doc:`fix wall/gran <fix_wall_gran>`
command.  At some point we plan to allow region surfaces to be used as
frictional walls, as well as triangulated surfaces.
