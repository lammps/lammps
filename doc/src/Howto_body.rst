Body particles
==============

**Overview:**

In LAMMPS, body particles are generalized finite-size particles.
Individual body particles can represent complex entities, such as
surface meshes of discrete points, collections of sub-particles,
deformable objects, etc.  Note that other kinds of finite-size
spherical and aspherical particles are also supported by LAMMPS, such
as spheres, ellipsoids, line segments, and triangles, but they are
simpler entities that body particles.  See the :doc:`Howto spherical <Howto_spherical>` doc page for a general overview of all
these particle types.

Body particles are used via the :doc:`atom_style body <atom_style>`
command.  It takes a body style as an argument.  The current body
styles supported by LAMMPS are as follows.  The name in the first
column is used as the *bstyle* argument for the :doc:`atom_style body <atom_style>` command.

+----------------------+---------------------------------------------------+
| *nparticle*          | rigid body with N sub-particles                   |
+----------------------+---------------------------------------------------+
| *rounded/polygon*    | 2d polygons with N vertices                       |
+----------------------+---------------------------------------------------+
| *rounded/polyhedron* | 3d polyhedra with N vertices, E edges and F faces |
+----------------------+---------------------------------------------------+

The body style determines what attributes are stored for each body and
thus how they can be used to compute pairwise body/body or
bond/non-body (point particle) interactions.  More details of each
style are described below.

More styles may be added in the future.  See the :doc:`Modify body <Modify_body>` doc page for details on how to add a new body
style to the code.


----------


**When to use body particles:**

You should not use body particles to model a rigid body made of
simpler particles (e.g. point, sphere, ellipsoid, line segment,
triangular particles), if the interaction between pairs of rigid
bodies is just the summation of pairwise interactions between the
simpler particles.  LAMMPS already supports this kind of model via the
:doc:`fix rigid <fix_rigid>` command.  Any of the numerous pair styles
that compute interactions between simpler particles can be used.  The
:doc:`fix rigid <fix_rigid>` command time integrates the motion of the
rigid bodies.  All of the standard LAMMPS commands for thermostatting,
adding constraints, performing output, etc will operate as expected on
the simple particles.

By contrast, when body particles are used, LAMMPS treats an entire
body as a single particle for purposes of computing pairwise
interactions, building neighbor lists, migrating particles between
processors, output of particles to a dump file, etc.  This means that
interactions between pairs of bodies or between a body and non-body
(point) particle need to be encoded in an appropriate pair style.  If
such a pair style were to mimic the :doc:`fix rigid <fix_rigid>` model,
it would need to loop over the entire collection of interactions
between pairs of simple particles within the two bodies, each time a
single body/body interaction was computed.

Thus it only makes sense to use body particles and develop such a pair
style, when particle/particle interactions are more complex than what
the :doc:`fix rigid <fix_rigid>` command can already calculate.  For
example, consider particles with one or more of the following
attributes:

* represented by a surface mesh
* represented by a collection of geometric entities (e.g. planes + spheres)
* deformable
* internal stress that induces fragmentation

For these models, the interaction between pairs of particles is likely
to be more complex than the summation of simple pairwise interactions.
An example is contact or frictional forces between particles with
planar surfaces that inter-penetrate.  Likewise, the body particle may
store internal state, such as a stress tensor used to compute a
fracture criterion.

These are additional LAMMPS commands that can be used with body
particles of different styles

+------------------------------------------------+-----------------------------------------------------+
| :doc:`fix nve/body <fix_nve_body>`             | integrate motion of a body particle in NVE ensemble |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`fix nvt/body <fix_nvt_body>`             | ditto for NVT ensemble                              |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`fix npt/body <fix_npt_body>`             | ditto for NPT ensemble                              |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`fix nph/body <fix_nph_body>`             | ditto for NPH ensemble                              |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`compute body/local <compute_body_local>` | store sub-particle attributes of a body particle    |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`compute temp/body <compute_temp_body>`   | compute temperature of body particles               |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`dump local <dump>`                       | output sub-particle attributes of a body particle   |
+------------------------------------------------+-----------------------------------------------------+
| :doc:`dump image <dump_image>`                 | output body particle attributes as an image         |
+------------------------------------------------+-----------------------------------------------------+

The pair styles defined for use with specific body styles are listed
in the sections below.


----------


**Specifics of body style nparticle:**

The *nparticle* body style represents body particles as a rigid body
with a variable number N of sub-particles.  It is provided as a
vanilla, prototypical example of a body particle, although as
mentioned above, the :doc:`fix rigid <fix_rigid>` command already
duplicates its functionality.

The atom\_style body command for this body style takes two additional
arguments:


.. parsed-literal::

   atom_style body nparticle Nmin Nmax
   Nmin = minimum # of sub-particles in any body in the system
   Nmax = maximum # of sub-particles in any body in the system

The Nmin and Nmax arguments are used to bound the size of data
structures used internally by each particle.

When the :doc:`read_data <read_data>` command reads a data file for this
body style, the following information must be provided for each entry
in the *Bodies* section of the data file:


.. parsed-literal::

   atom-ID 1 M
   N
   ixx iyy izz ixy ixz iyz
   x1 y1 z1
   ...
   xN yN zN

where M = 6 + 3\*N, and N is the number of sub-particles in the body
particle.

The integer line has a single value N.  The floating point line(s)
list 6 moments of inertia followed by the coordinates of the N
sub-particles (x1 to zN) as 3N values.  These values can be listed on
as many lines as you wish; see the :doc:`read_data <read_data>` command
for more details.

The 6 moments of inertia (ixx,iyy,izz,ixy,ixz,iyz) should be the
values consistent with the current orientation of the rigid body
around its center of mass.  The values are with respect to the
simulation box XYZ axes, not with respect to the principal axes of the
rigid body itself.  LAMMPS performs the latter calculation internally.
The coordinates of each sub-particle are specified as its x,y,z
displacement from the center-of-mass of the body particle.  The
center-of-mass position of the particle is specified by the x,y,z
values in the *Atoms* section of the data file, as is the total mass
of the body particle.

The :doc:`pair_style body/nparticle <pair_body_nparticle>` command can be used
with this body style to compute body/body and body/non-body interactions.

For output purposes via the :doc:`compute body/local <compute_body_local>` and :doc:`dump local <dump>`
commands, this body style produces one datum for each of the N
sub-particles in a body particle.  The datum has 3 values:


.. parsed-literal::

   1 = x position of sub-particle
   2 = y position of sub-particle
   3 = z position of sub-particle

These values are the current position of the sub-particle within the
simulation domain, not a displacement from the center-of-mass (COM) of
the body particle itself.  These values are calculated using the
current COM and orientation of the body particle.

For images created by the :doc:`dump image <dump_image>` command, if the
*body* keyword is set, then each body particle is drawn as a
collection of spheres, one for each sub-particle.  The size of each
sphere is determined by the *bflag1* parameter for the *body* keyword.
The *bflag2* argument is ignored.


----------


**Specifics of body style rounded/polygon:**

The *rounded/polygon* body style represents body particles as a 2d
polygon with a variable number of N vertices.  This style can only be
used for 2d models; see the :doc:`boundary <boundary>` command.  See the
"pair\_style body/rounded/polygon" doc page for a diagram of two
squares with rounded circles at the vertices.  Special cases for N = 1
(circle) and N = 2 (rod with rounded ends) can also be specified.

One use of this body style is for 2d discrete element models, as
described in :ref:`Fraige <body-Fraige>`.

Similar to body style *nparticle*\ , the atom\_style body command for
this body style takes two additional arguments:


.. parsed-literal::

   atom_style body rounded/polygon Nmin Nmax
   Nmin = minimum # of vertices in any body in the system
   Nmax = maximum # of vertices in any body in the system

The Nmin and Nmax arguments are used to bound the size of data
structures used internally by each particle.

When the :doc:`read_data <read_data>` command reads a data file for this
body style, the following information must be provided for each entry
in the *Bodies* section of the data file:


.. parsed-literal::

   atom-ID 1 M
   N
   ixx iyy izz ixy ixz iyz
   x1 y1 z1
   ...
   xN yN zN
   i j j k k ...
   diameter

where M = 6 + 3\*N + 2\*N + 1, and N is the number of vertices in the
body particle.

The integer line has a single value N.  The floating point line(s)
list 6 moments of inertia followed by the coordinates of the N
vertices (x1 to zN) as 3N values (with z = 0.0 for each), followed by
2N vertex indices corresponding to the end points of the N edges,
followed by a single diameter value = the rounded diameter of the
circle that surrounds each vertex. The diameter value can be different
for each body particle. These floating-point values can be listed on
as many lines as you wish; see the :doc:`read_data <read_data>` command
for more details.

The 6 moments of inertia (ixx,iyy,izz,ixy,ixz,iyz) should be the
values consistent with the current orientation of the rigid body
around its center of mass.  The values are with respect to the
simulation box XYZ axes, not with respect to the principal axes of the
rigid body itself.  LAMMPS performs the latter calculation internally.
The coordinates of each vertex are specified as its x,y,z displacement
from the center-of-mass of the body particle.  The center-of-mass
position of the particle is specified by the x,y,z values in the
*Atoms* section of the data file.

For example, the following information would specify a square particle
whose edge length is sqrt(2) and rounded diameter is 1.0.  The
orientation of the square is aligned with the xy coordinate axes which
is consistent with the 6 moments of inertia: ixx iyy izz ixy ixz iyz =
1 1 4 0 0 0. Note that only Izz matters in 2D simulations.


.. parsed-literal::

   3 1 27
   4
   1 1 4 0 0 0
   -0.7071 -0.7071 0
   -0.7071 0.7071 0
   0.7071 0.7071 0
   0.7071 -0.7071 0
   0 1
   1 2
   2 3
   3 0
   1.0

A rod in 2D, whose length is 4.0, mass 1.0, rounded at two ends
by circles of diameter 0.5, is specified as follows:


.. parsed-literal::

   1 1 13
   2
   1 1 1.33333 0 0 0
   -2 0 0
   2 0 0
   0.5

A disk, whose diameter is 3.0, mass 1.0, is specified as follows:


.. parsed-literal::

   1 1 10
   1
   1 1 4.5 0 0 0
   0 0 0
   3.0

The :doc:`pair_style body/rounded/polygon <pair_body_rounded_polygon>`
command can be used with this body style to compute body/body
interactions.  The :doc:`fix wall/body/polygon <fix_wall_body_polygon>`
command can be used with this body style to compute the interaction of
body particles with a wall.


----------


**Specifics of body style rounded/polyhedron:**

The *rounded/polyhedron* body style represents body particles as a 3d
polyhedron with a variable number of N vertices, E edges and F faces.
This style can only be used for 3d models; see the
:doc:`boundary <boundary>` command.  See the "pair\_style
body/rounded/polygon" doc page for a diagram of a two 2d squares with
rounded circles at the vertices.  A 3d cube with rounded spheres at
the 8 vertices and 12 rounded edges would be similar.  Special cases
for N = 1 (sphere) and N = 2 (rod with rounded ends) can also be
specified.

This body style is for 3d discrete element models, as described in
:ref:`Wang <body-Wang>`.

Similar to body style *rounded/polygon*\ , the atom\_style body command
for this body style takes two additional arguments:


.. parsed-literal::

   atom_style body rounded/polyhedron Nmin Nmax
   Nmin = minimum # of vertices in any body in the system
   Nmax = maximum # of vertices in any body in the system

The Nmin and Nmax arguments are used to bound the size of data
structures used internally by each particle.

When the :doc:`read_data <read_data>` command reads a data file for this
body style, the following information must be provided for each entry
in the *Bodies* section of the data file:


.. parsed-literal::

   atom-ID 3 M
   N E F
   ixx iyy izz ixy ixz iyz
   x1 y1 z1
   ...
   xN yN zN
   0 1
   1 2
   2 3
   ...
   0 1 2 -1
   0 2 3 -1
   ...
   1 2 3 4
   diameter

where M = 6 + 3\*N + 2\*E + 4\*F + 1, and N is the number of vertices in
the body particle, E = number of edges, F = number of faces.

The integer line has three values: number of vertices (N), number of
edges (E) and number of faces (F). The floating point line(s) list 6
moments of inertia followed by the coordinates of the N vertices (x1
to zN) as 3N values, followed by 2N vertex indices corresponding to
the end points of the E edges, then 4\*F vertex indices defining F
faces.  The last value is the diameter value = the rounded diameter of
the sphere that surrounds each vertex. The diameter value can be
different for each body particle. These floating-point values can be
listed on as many lines as you wish; see the
:doc:`read_data <read_data>` command for more details.  Because the
maximum number of vertices per face is hard-coded to be 4
(i.e. quadrilaterals), faces with more than 4 vertices need to be
split into triangles or quadrilaterals.  For triangular faces, the
last vertex index should be set to -1.

The ordering of the 4 vertices within a face should follow
the right-hand rule so that the normal vector of the face points
outwards from the center of mass.

The 6 moments of inertia (ixx,iyy,izz,ixy,ixz,iyz) should be the
values consistent with the current orientation of the rigid body
around its center of mass.  The values are with respect to the
simulation box XYZ axes, not with respect to the principal axes of the
rigid body itself.  LAMMPS performs the latter calculation internally.
The coordinates of each vertex are specified as its x,y,z displacement
from the center-of-mass of the body particle.  The center-of-mass
position of the particle is specified by the x,y,z values in the
*Atoms* section of the data file.

For example, the following information would specify a cubic particle
whose edge length is 2.0 and rounded diameter is 0.5.
The orientation of the cube is aligned with the xyz coordinate axes
which is consistent with the 6 moments of inertia: ixx iyy izz ixy ixz
iyz = 0.667 0.667 0.667 0 0 0.


.. parsed-literal::

   1 3 79
   8 12 6
   0.667 0.667 0.667 0 0 0
   1 1 1
   1 -1 1
   -1 -1 1
   -1 1 1
   1 1 -1
   1 -1 -1
   -1 -1 -1
   -1 1 -1
   0 1
   1 2
   2 3
   3 0
   4 5
   5 6
   6 7
   7 4
   0 4
   1 5
   2 6
   3 7
   0 1 2 3
   4 5 6 7
   0 1 5 4
   1 2 6 5
   2 3 7 6
   3 0 4 7
   0.5

A rod in 3D, whose length is 4.0, mass 1.0 and rounded at two ends
by circles of diameter 0.5, is specified as follows:


.. parsed-literal::

   1 1 13
   2
   0 1.33333 1.33333 0 0 0
   -2 0 0
   2 0 0
   0.5

A sphere whose diameter is 3.0 and mass 1.0, is specified as follows:


.. parsed-literal::

   1 1 10
   1
   0.9 0.9 0.9 0 0 0
   0 0 0
   3.0

The :doc:`pair_style body/rounded/polhedron <pair_body_rounded_polyhedron>` command can
be used with this body style to compute body/body interactions.  The
:doc:`fix wall/body/polyhedron <fix_wall_body_polygon>` command can be
used with this body style to compute the interaction of body particles
with a wall.


----------


For output purposes via the :doc:`compute body/local <compute_body_local>` and :doc:`dump local <dump>`
commands, this body style produces one datum for each of the N
sub-particles in a body particle.  The datum has 3 values:


.. parsed-literal::

   1 = x position of vertex
   2 = y position of vertex
   3 = z position of vertex

These values are the current position of the vertex within the
simulation domain, not a displacement from the center-of-mass (COM) of
the body particle itself.  These values are calculated using the
current COM and orientation of the body particle.

For images created by the :doc:`dump image <dump_image>` command, if the
*body* keyword is set, then each body particle is drawn as a polygon
consisting of N line segments.  Note that the line segments are drawn
between the N vertices, which does not correspond exactly to the
physical extent of the body (because the :doc:`pair_style rounded/polygon <pair_body_rounded_polygon>` defines finite-size
spheres at those point and the line segments between the spheres are
tangent to the spheres).  The drawn diameter of each line segment is
determined by the *bflag1* parameter for the *body* keyword.  The
*bflag2* argument is ignored.


----------


.. _body-Fraige:



**(Fraige)** F. Y. Fraige, P. A. Langston, A. J. Matchett, J. Dodds,
Particuology, 6, 455 (2008).

.. _body-Wang:



**(Wang)** J. Wang, H. S. Yu, P. A. Langston, F. Y. Fraige, Granular
Matter, 13, 1 (2011).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
