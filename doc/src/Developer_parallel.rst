Parallel algorithms
-------------------

LAMMPS is from ground up designed to be running in parallel using the
MPI standard with distributed data via domain decomposition.  The
parallelization has to be efficient to enable good strong scaling (=
good speedup for the same system) and good weak scaling (= the
computational cost of enlarging the system is linear with the system
size).  Additional parallelization using GPUs or OpenMP can then be
applied within the sub-domain assigned to an MPI process.

Partitioning
^^^^^^^^^^^^

The underlying spatial decomposition strategy used by LAMMPS for
distributed-memory parallelism is set with the :doc:`comm_style command
<comm_style>` and can be either "brick" (a regular grid) or "tiled".

.. _domain-decomposition:
.. figure:: img/domain-decomp.png
   :align: center

   domain decomposition

   This figure shows the different kinds of domain decomposition used
   for MPI parallelization: "brick" on the left with an orthogonal
   (left) and a triclinic (middle) simulation domain, and a "tiled"
   decomposition (right).  The black lines show the division into
   sub-domains and the contained atoms are "owned" by the corresponding
   MPI process. The green dashed lines indicate how sub-domains are
   extended with "ghost" atoms up to the communication cutoff distance.

The LAMMPS simulation box is a 3d or 2d volume, which can be orthogonal
or triclinic in shape, as illustrated in the :ref:`domain-decomposition`
figure for the 2d case.  Orthogonal means the box edges are aligned with
the *x*, *y*, *z* Cartesian axes, and the box faces are thus all
rectangular.  Triclinic allows for a more general parallelepiped shape
in which edges are aligned with three arbitrary vectors and the box
faces are parallelograms.  In each dimension box faces can be periodic,
or non-periodic with fixed or shrink-wrapped boundaries.  In the fixed
case, atoms which move outside the face are deleted; shrink-wrapped
means the position of the box face adjusts continuously to enclose all
the atoms.

For distributed-memory MPI parallelism, the simulation box is spatially
decomposed (partitioned) into non-overlapping sub-domains which fill the
box. The default partitioning, "brick", is most suitable when atom
density is roughly uniform, as shown in the left-side images of the
:ref:`domain-decomposition` figure.  The sub-domains comprise a regular
grid and all sub-domains are identical in size and shape.  Both the
orthogonal and triclinic boxes can deform continuously during a
simulation, e.g. to compress a solid or shear a liquid, in which case
the processor sub-domains likewise deform.


For models with non-uniform density, the number of particles per
processor can be load-imbalanced with the default partitioning.  This
reduces parallel efficiency, as the overall simulation rate is limited
by the slowest processor, i.e. the one with the largest computational
load.  For such models, LAMMPS supports multiple strategies to reduce
the load imbalance:

- The processor grid decomposition is by default based on the simulation
  cell volume and tries to optimize the volume to surface ratio for the sub-domains.
  This can be changed with the :doc:`processors command <processors>`.
- The parallel planes defining the size of the sub-domains can be shifted
  with the :doc:`balance command <balance>`. Which can be done in addition
  to choosing a more optimal processor grid.
- The recursive bisectioning algorithm in combination with the "tiled"
  communication style can produce a partitioning with equal numbers of
  particles in each sub-domain.


.. |decomp1| image:: img/decomp-regular.png
   :width: 24%

.. |decomp2| image:: img/decomp-processors.png
   :width: 24%

.. |decomp3| image:: img/decomp-balance.png
   :width: 24%

.. |decomp4| image:: img/decomp-rcb.png
   :width: 24%

|decomp1|  |decomp2|  |decomp3|  |decomp4|

The pictures above demonstrate different decompositions for a 2d system
with 12 MPI ranks.  The atom colors indicate the load imbalance with
green being optimal and red the least optimal.  Due to the vacuum in the system, the default
decomposition is unbalanced with several MPI ranks without atoms
(left). By forcing a 1x12x1 processor grid, every MPI rank does
computations now, but number of atoms per sub-domain is still uneven and
the thin slice shape increases the amount of communication between sub-domains
(center left). With a 2x6x1 processor grid and shifting the
sub-domain divisions, the load imbalance is further reduced and the amount
of communication required between sub-domains is less (center right).
And using the recursive bisectioning leads to further improved
decomposition (right).


Communication
^^^^^^^^^^^^^

Following the partitioning scheme the data of the system is distributed
and each MPI process stores information (positions, velocities, etc.)
for the subset of atoms within its sub-domain, called "owned" atoms.  It
also stores copies of some of that information for "ghost" atoms within
the communication cutoff distance of its sub-domain, which are owned by
nearby MPI processes. This enables calculating all short-range
interactions which involve atoms the MPI process "owns" in parallel.
The dashed-line boxes in the :ref:`domain-decomposition` figure
illustrate the extended ghost-atom sub-domain for one processor.

This approach is also used to implement periodic boundary conditions:
atoms that lie within the cutoff distance across a periodic boundary are
also stored as ghost atoms and taken from the periodic replication of
the sub-domain, which may be the same sub-domain, e.g. if running in
serial.  As a consequence of this, force computation in LAMMPS is not
subject to minimum image conventions and thus cutoffs may be larger than
half the simulation domain.

.. _ghost-atom-comm:
.. figure:: img/ghost-comm.png
   :align: center

   ghost atom communication

   This figure shows the ghost atom communication patterns between
   sub-domains for "brick" (left) and "tiled" communication styles for
   2d simulations.  The numbers indicate MPI process ranks.  Here the
   sub-domains are drawn spatially separated for clarity.  The
   dashed-line box is the extended sub-domain of processor 0 which
   includes its ghost atoms.  The red- and blue-shaded boxes are the
   regions of communicated ghost atoms.

The diagrams of the `ghost-atom-comm` figure illustrate how ghost atom
communication is performed in two stages for a 2d simulation (three in
3d) for both a regular and irregular partitioning of the simulation box.
For the regular case (left) atoms are exchanged first in the
*x*-direction, then in *y*, with four neighbors in the grid of processor
sub-domains.

In the *x* stage, processor ranks 1 and 2 send owned atoms in their
red-shaded regions to rank 0 (and vice versa).  Then in the *y* stage,
ranks 3 and 4 send atoms in their blue-shaded regions to rank 0, which
includes ghost atoms they received in the *x* stage.  Rank 0 thus
acquires all its ghost atoms; atoms in the solid blue corner regions
are communicated twice before rank 0 receives them.

For the irregular case (right) the two stages are similar, but a
processor can have more than one neighbor in each direction.  In the
*x* stage, MPI ranks 1,2,3 send owned atoms in their red-shaded regions to
rank 0 (and vice versa).  These include only atoms between the lower
and upper *y*-boundary of rank 0's sub-domain.  In the *y* stage, ranks
4,5,6 send atoms in their blue-shaded regions to rank 0.  This may
include ghost atoms they received in the *x* stage, but only if they
are needed by rank 0 to fill its extended ghost atom regions in the
+/-*y* directions (blue rectangles).  Thus in this case, ranks 5 and
6 do not include ghost atoms they received from each other (in the *x*
stage) in the atoms they send to rank 0.  The key point is that while
the pattern of communication is more complex in the irregular
partitioning case, it can still proceed in two stages (three in 3d)
via atom exchanges with only neighboring processors.

When attributes of owned atoms are sent to neighboring processors to
become attributes of their ghost atoms, LAMMPS calls this a "forward"
communication.  On timesteps when atoms migrate to new owning processors
and neighbor lists are rebuilt, each processor creates a list of its
owned atoms which are ghost atoms in each of its neighbor processors.
These lists are used to pack per-atom coordinates (for example) into
message buffers in subsequent steps until the next reneighboring.

A "reverse" communication is when computed ghost atom attributes are
sent back to the processor who owns the atom.  This is used (for
example) to sum partial forces on ghost atoms to the complete force on
owned atoms.  The order of the two stages described in the
:ref:`ghost-atom-comm` figure is inverted and the same lists of atoms
are used to pack and unpack message buffers with per-atom forces.  When
a received buffer is unpacked, the ghost forces are summed to owned atom
forces.  As in forward communication, forces on atoms in the four blue
corners of the diagrams are sent, received, and summed twice (once at
each stage) before owning processors have the full force.

These two operations are used many places within LAMMPS aside from
exchange of coordinates and forces, for example by manybody potentials
to share intermediate per-atom values, or by rigid-body integrators to
enable each atom in a body to access body properties.  Here are
additional details about how these communication operations are
performed in LAMMPS:

- When exchanging data with different processors, forward and reverse
  communication is done using ``MPI_Send()`` and ``MPI_IRecv()`` calls.
  If a processor is "exchanging" atoms with itself, only the pack and
  unpack operations are performed, e.g. to create ghost atoms across
  periodic boundaries when running on a single processor.

- For forward communication of owned atom coordinates, periodic box
  lengths are added and subtracted when the receiving processor is
  across a periodic boundary from the sender.  There is then no need to
  apply a minimum image convention when calculating distances between
  atom pairs when building neighbor lists or computing forces.

- The cutoff distance for exchanging ghost atoms is typically equal to
  the neighbor cutoff.  But it can also chosen to be longer if needed,
  e.g. half the diameter of a rigid body composed of multiple atoms or
  over 3x the length of a stretched bond for dihedral interactions.  It
  can also exceed the periodic box size.  For the regular communication
  pattern (left), if the cutoff distance extends beyond a neighbor
  processor's sub-domain, then multiple exchanges are performed in the
  same direction.  Each exchange is with the same neighbor processor,
  but buffers are packed/unpacked using a different list of atoms. For
  forward communication, in the first exchange a processor sends only
  owned atoms.  In subsequent exchanges, it sends ghost atoms received
  in previous exchanges.  For the irregular pattern (right) overlaps of
  a processor's extended ghost-atom sub-domain with all other processors
  in each dimension are detected.

Neighbor lists
^^^^^^^^^^^^^^

Long-range interactions
^^^^^^^^^^^^^^^^^^^^^^^
