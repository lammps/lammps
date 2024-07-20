Granular surfaces
=================

As explained on the :doc:`Howto granular <Howto_granular>` doc page,
granular systems are composed of spherical particles with a diameter,
as opposed to point particles.  This means they have an angular
velocity and torque can be imparted to them to cause them to rotate.

The :doc:`Howto granular <Howto_granular>` doc page lists various
atom, pair, fix, and compute styles useful for creaeting granular
models.

This page explains how you can also define granular surfaces which are
a collection of triangles (3d systems) or line segments (2d systems),
which act as boundaries interacting with the particles.  Different
kinds of particle/surface interactions can be specified with similar
options as the pair styles listed on the :doc:`Howto granular
<Howto_granular>` doc page.

----------

Global versus local surfaces
""""""""""""""""""""""""""""
A key point to understand is that LAMMPS supports two forms of
granular surfaces.  You cannot use both in the same simulation.

The first is *global* which means that each processor stores a copy of
all the triangles/lines.  This is suitable when a modest number of
triangles/lines is needed.  They can be large triangles/lines, any of
which span a significant fraction of the simulation box size in one or
more dimensions.

The second is *local* which means that the collection of
triangles/lines is distributed across processers the same what that
particles are distributed.  Each processor is assigned to sub-domain
of the simulation box and owns whichever triangles/lines have their
center point in the processor's sub-domain.  Similar to particles,
each processor may also own ghost copies of triangles/lines whose
finite size overlaps with the processor's sub-domain.  The number of
triangles/lines can now be very large.  For effective distribution and
minimal communication all the triangles/lines should be small, no more
than a few particle diameters in size.  If even one larger triangle or
line is defined then the neighbor list cutoff and communication cutoff
will be set correspondingly larger, which can slow down the
simulation.

One of these two commands must be specified to use *global* or *local*
surfaces in your granular simulation:

* :doc:`fix surface/global <fix_surface_global>`
* :doc:`fix surface/local <fix_surface_local>`

The :doc:`fix surface/global <fix_surface_global>` command reads in
the global surfaces in one of two ways.  The first option is from a
molecule file(s) previously defined by the :doc:`molecule <molecule>`
command.  The file should define triangles or lines with header
keywords and a Triangles or Lines section.  The second option is from
a text or binary STL (strereolithogray) file which defines a set of
triangles.  It can only be used with 3d simulations.

The :doc:`fix surface/local <fix_surface_local>` command defines local
surface in one of three ways.  The first two options are the same
molecule and STL files explained in the previous paragraph.  In this
case, the list of triangles/lines is distributed across processors
based on the center point of each triangle/line.  The third option is
to include them in a LAMMPS data file which has been previously read
in via the :doc:`read_data <read_data>` command.  If the file has a
Triangles or Lines section, then triangles/lines will be read in and
distributed along with any particles the data file includes.

----------

Atom styles for granular surfaces
"""""""""""""""""""""""""""""""""

For all three ways of defining *local* surfaces, the triangles/lines
are stored internally in LAMMPS as triangle-style or line-style
particles.  This means you must use a hybrid atom style which includes
one of these two sub-styles (for 3d or 2d):

* :doc:`atom_style tri <atom_style>` for 3d simulations
* :doc:`atom_style line <atom_style>` for 2d simulations

Typically, the atom_style hybrid command will also define a
:doc:`atom_style sphere <atom_style>` sub-style for the granular
particles which interact with the surfaces.

Note that for molecule or STL file input, the :doc:`fix surface/local
<fix_surface_local>` command reads the file(s) and uses the values for
each surface to creat a single new triangle or line particle.  For
data file input, the triangle/line particles are created when the data
file is read.

For granular simluations with *global* surfaces, a hybrid atom style
which defines triangle-style or line-style particles should NOT be
used.  The triangles/lines are stored by the :doc:`fix surface/global
<fix_surface_global>` command and not as triangle-style or line-style
particles.

----------

Rules for surface topology
""""""""""""""""""""""""""

For both *global* and *local* surfaces, granular particles interact
with both sides of each triangle or line segment.

No check is made to see if two triangles or line segments intersect
each other; this is allowed if it makes sense for the geometry of the
collection of surfaces.

As an example, consider a 2d simulation which mixes a container of
granular particles.  *Global* line segments are used to define both
the box-shaped container and the mixer in the center.  The 4 mixer
blades are in the shape of a large X and are made to rotate using the
:doc:`fix_modify <fix_modify>` command (see below).

The 2 blades could be defined by 2 line segments which cross each
other at their centers.  Or the 2 blades could be defined by 4 line
segments, all of which have a common endpoint at the center of the
mixer.  Or the 2 blades could be defined by 4 non-touching line
segments, all of which have a distinct endpoint near the center of the
mixer, but displaced from it by a distance less than the radius of a
granular particle.  In any of these cases, when a particle gets very
close to the center of the mixer it will interact with both nearby
line segments as expected.

See the next section on connectivity for how two triangles or line
segemnts are treated if they share a common edge (triangle) or point
(triange or line).

----------

Surface connectivity
""""""""""""""""""""

If multiple triangles/lines are used to define a contiguous surface
which is flat or gently curved or has sharp edges or corners, LAMMPS
will detect when two or more line segments (2d) share the same
endpoint.  Or when two or more triangles (3d) share the same edge or
same corner point.

This connectivity is stored internally and is used when appropriate to
calculate accurate forces on particles which simultaneously overlap
with 2 or more connected triangles or line segments.

Consider the simulation model of the previous section for a 2d mixer
now defined by *local* line segments.  The flat surface of each mixer
blade (and container box faces) is defined by multiple small line
segments.  It is imporant that these line segments be "connected" so
that when a particle contacts two adjacent line segments at the same
time, the resulting force on the particle is the same as it would be
if it were contacting the middle of a single long line segment.

Here is how to ensure that LAMMPS detects the appropriate connections.

For either *global* or *local* surfaces, if the triangles/lines are
defined in a molecule or STL file, then 3 corner points (triangle) or
2 end points (line) will be listed for each triangle/line in the file.
LAMMPS will only make a connection between 2 triangles or lines if a
shared point is EXACTLY the same in both.  This is a single point in
both for a corner point or end point connection.  It is two points in
both triangles for an edge connection.

For *local* surfaces, if the triangles/lines are defined in a data
file, then 3 corner points (triangle) or 2 end points (line) will be
listed for each triangle/line in the file.  However in this case,
LAMMPS will allow for an INEXACT match of a shared point to make a
connection between 2 triangles or lines.  Again, this is a single
point in both for a corner point or end point connection.  It is two
points in both triangles for an edge connection.

An INEXACT match means that the two points can be EPSILON apart.
EPSILON is defined as a tiny fraction (1.0e-4) of the size of
the smallest triangle or line in the system.

The reason INEXACT matches are allowed is that data files can be
created in a variety of manners, including by LAMMPS itself as a
simulation runs via the :doc:`write_data <write_data>` command.
Interally, triangle-style and line-style particles do not store their
corner points directly.  Instead, the center point of the
triangle/line is stored, along with an orientation of the
triangle/line and a displacement vector from the center point for each
corner point.  This means that when new corner points values are
written to a data file for two different triangles/line, they may
differ by epsilon due to round-offs in finite-precision arithmetic.

Note that due to how connectivity is defined, two triangles/lines will
not be connected if their corner points are separted by even small
distances (greater than EPSILON).  Likewise they will not be connected
if the corner point of one triangle/line is very close to (or even on)
the surface of another triangle or middle of another line segment.  In
general these kinds of granular surfaces could be problematic and
should be avoided, but LAMMPS does not check for these conditions.

NOTE: maybe add a picture of T-shaped surf with 2 line segments (not
3).  Explain why it could be bad?

Note that if a triangle or line segment has a free edge or free
corner/end point (not connected to any other triangle/line), granular
particles will still interact with the triangle/line if the nearest
contact point to the spherical particle center is on the free edge or
is the free corner/end point.

----------

Surface motion
""""""""""""""

By default, surface triangles/lines are motionless during a
simulation, whether they are *global* or *local*.  Triangles/lines
impart forces and torques to granular particles, but the inverse
forces/torques on the triangles/lines do not cause them to move.

However, triangles/lines can be made to move in a prescribed manner.
E.g. the rotation of 2d mixer blades in the example described above.
These two commands can be used for that purpose:

* :doc:`fix_modify move <fix_modify>` for *global* surfaces
* :doc:`fix move <fix_move>` for *local* surfaces

For *global* surfaces, the :doc:`fix_modify move <fix_modify>` command
can rotate all the surfaces around a specified axis at a specified
rate.

For *local* surfaces, the :doc:`fix move <fix_move>` command can move
a specified subset of the triangles/lines in various ways
(translation, rotation, etc).

More options for moving surfaces in prescribed manners will likely be
added in the future.

----------

Example scripts
"""""""""""""""

The examples/gransurf directory has example input scripts which use
both *global* and *local* surfaces.  Both 2d and 3d models are included.

Each script produces a series of snapshot images using the :doc:`dump
image <dump_image>` command.  The snapshots visualize both the
particles and granular surfaces.  The snaphost can be animated to view
a movie of the simulation.
