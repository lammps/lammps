.. index:: fix wall/region

fix wall/region command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID wall/region region-ID style args ... cutoff

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/region = style name of this fix command
* region-ID = region whose boundary will act as wall
* style = *lj93* or *lj126* or *lj1043* or *colloid* or *harmonic* or *morse*
* args for styles *lj93* or *lj126* or *lj1043* or *colloid* or *harmonic* =

  .. parsed-literal::

        epsilon = strength factor for wall-particle interaction (energy or energy/distance\^2 units)
        sigma = size factor for wall-particle interaction (distance units)

* args for style *morse* =

  .. parsed-literal::

        D_0 = depth of the potential (energy units)
        alpha = width parameter (1/distance units)
        r_0 = distance of the potential minimum from wall position (distance units)

* cutoff = distance from wall at which wall-particle interaction is cut off (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix wall all wall/region mySphere lj93 1.0 1.0 2.5
   fix wall all wall/region mySphere harmonic 1.0 0.0 2.5
   fix wall all wall/region box_top morse 1.0 1.0 1.5 3.0

Description
"""""""""""

Treat the surface of the geometric region defined by the *region-ID*
as a bounding wall which interacts with nearby particles according to
the specified style.

The distance between a particle and the surface is the distance to the
nearest point on the surface and the force the wall exerts on the
particle is along the direction between that point and the particle,
which is the direction normal to the surface at that point.  Note that
if the region surface is comprised of multiple "faces", then each face
can exert a force on the particle if it is close enough.  E.g. for
:doc:`region_style block <region>`, a particle in the interior, near a
corner of the block, could feel wall forces from 1, 2, or 3 faces of
the block.

Regions are defined using the :doc:`region <region>` command.  Note that
the region volume can be interior or exterior to the bounding surface,
which will determine in which direction the surface interacts with
particles, i.e. the direction of the surface normal.  The surface of
the region only exerts forces on particles "inside" the region; if a
particle is "outside" the region it will generate an error, because it
has moved through the wall.

Regions can either be primitive shapes (block, sphere, cylinder, etc)
or combinations of primitive shapes specified via the *union* or
*intersect* region styles.  These latter styles can be used to
construct particle containers with complex shapes.  Regions can also
change over time via the :doc:`region <region>` command keywords (move)
and *rotate*\ .  If such a region is used with this fix, then the of
region surface will move over time in the corresponding manner.

.. note::

   As discussed on the :doc:`region <region>` command doc page,
   regions in LAMMPS do not get wrapped across periodic boundaries.  It
   is up to you to insure that periodic or non-periodic boundaries are
   specified appropriately via the :doc:`boundary <boundary>` command when
   using a region as a wall that bounds particle motion.  This also means
   that if you embed a region in your simulation box and want it to
   repulse particles from its surface (using the "side out" option in the
   :doc:`region <region>` command), that its repulsive force will not be
   felt across a periodic boundary.

.. note::

   For primitive regions with sharp corners and/or edges (e.g. a
   block or cylinder), wall/particle forces are computed accurately for
   both interior and exterior regions.  For *union* and *intersect*
   regions, additional sharp corners and edges may be present due to the
   intersection of the surfaces of 2 or more primitive volumes.  These
   corners and edges can be of two types: concave or convex.  Concave
   points/edges are like the corners of a cube as seen by particles in
   the interior of a cube.  Wall/particle forces around these features
   are computed correctly.  Convex points/edges are like the corners of a
   cube as seen by particles exterior to the cube, i.e. the points jut
   into the volume where particles are present.  LAMMPS does NOT compute
   the location of these convex points directly, and hence wall/particle
   forces in the cutoff volume around these points suffer from
   inaccuracies.  The basic problem is that the outward normal of the
   surface is not continuous at these points.  This can cause particles
   to feel no force (they don't "see" the wall) when in one location,
   then move a distance epsilon, and suddenly feel a large force because
   they now "see" the wall.  In a worst-case scenario, this can blow
   particles out of the simulation box.  Thus, as a general rule you
   should not use the fix wall/gran/region command with *union* or
   *interesect* regions that have convex points or edges resulting from
   the union/intersection (convex points/edges in the union/intersection
   due to a single sub-region are still OK).

.. note::

   Similarly, you should not define *union* or *intersert* regions
   for use with this command that share an overlapping common face that
   is part of the overall outer boundary (interior boundary is OK), even
   if the face is smooth.  E.g. two regions of style block in a *union*
   region, where the two blocks overlap on one or more of their faces.
   This is because LAMMPS discards points that are part of multiple
   sub-regions when calculating wall/particle interactions, to avoid
   double-counting the interaction.  Having two coincident faces could
   cause the face to become invisible to the particles.  The solution is
   to make the two faces differ by epsilon in their position.

The energy of wall-particle interactions depends on the specified
style.

For style *lj93*\ , the energy E is given by the 9/3 potential:

.. math::

 E = \epsilon \left[ \frac{2}{15} \left(\frac{\sigma}{r}\right)^{9} -
                       \left(\frac{\sigma}{r}\right)^3 \right]
                       \qquad r < r_c

For style *lj126*\ , the energy E is given by the 12/6 potential:

.. math::

 E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

For style *wall/lj1043*\ , the energy E is given by the 10/4/3 potential:

.. math::

 E = 2 \pi \epsilon \left[ \frac{2}{5} \left(\frac{\sigma}{r}\right)^{10} -
                       \left(\frac{\sigma}{r}\right)^4 -
                       \frac{\sqrt(2)\sigma^3}{3\left(r+\left(0.61/\sqrt(2)\right)\sigma\right)^3}\right]
                       \qquad r < r_c

For style *colloid*\ , the energy E is given by an integrated form of
the :doc:`pair_style colloid <pair_colloid>` potential:

.. math::

   E = & \epsilon \left[ \frac{\sigma^{6}}{7560}
   \left(\frac{6R-D}{D^{7}} + \frac{D+8R}{(D+2R)^{7}} \right) \right. \\
    & \left. - \frac{1}{6} \left(\frac{2R(D+R) + D(D+2R)
    \left[ \ln D - \ln (D+2R) \right]}{D(D+2R)} \right) \right] \qquad r < r_c

For style *wall/harmonic*\ , the energy E is given by a harmonic spring
potential (the distance parameter is ignored):

.. math::

   E = \epsilon \quad (r - r_c)^2 \qquad r < r_c

For style *wall/morse*\ , the energy E is given by the Morse potential:

.. math::

   E = D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
       \qquad r < r_c

Unlike other styles, this requires three parameters (:math:`D_0`,
:math:`\alpha`, and :math:`r_0` in this order) instead of two like
for the other wall styles.

In all cases, *r* is the distance from the particle to the region
surface, and Rc is the *cutoff* distance at which the particle and
surface no longer interact.  The cutoff is always the last argument.
The energy of the wall potential is shifted so that the wall-particle
interaction energy is 0.0 at the cutoff distance.

For a full description of these wall styles, see fix\_style
:doc:`wall <fix_wall>`

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy of interaction between atoms and the wall to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the interaction between
atoms and each wall to the system's virial as part of :doc:`thermodynamic output <thermo_style>`. The default is *virial no*

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar energy and a global 3-length vector
of forces, which can be accessed by various :doc:`output commands <Howto_output>`.  The scalar energy is the sum of energy
interactions for all particles interacting with the wall represented
by the region surface.  The 3 vector quantities are the x,y,z
components of the total force acting on the wall due to the particles.
The scalar and vector values calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the atom/wall interaction energy to be included in
   the total potential energy of the system (the quantity being
   minimized), you MUST enable the :doc:`fix_modify <fix_modify>` *energy*
   option for this fix.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix wall/lj93 <fix_wall>`,
:doc:`fix wall/lj126 <fix_wall>`,
:doc:`fix wall/lj1043 <fix_wall>`,
:doc:`fix wall/colloid <fix_wall>`,
:doc:`fix wall/harmonic <fix_wall>`,
:doc:`fix wall/gran <fix_wall_gran>`

**Default:** none
