Parallel algorithms
-------------------

LAMMPS is from ground up designed to be running in parallel using the
MPI standard with distributed data via domain decomposition.  The
parallelization has to be efficient to enable good strong scaling (=
good speedup for the same system) and good weak scaling (= the
computational cost of enlarging the system is linear with the system
size).  Additional parallelization using GPUs or OpenMP can then be
applied within the sub-domain assigned to an MPI process.  For clarity,
most of the following illustrations show the 2d simulation case.

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
with 12 MPI ranks.  The atom colors indicate the load imbalance of each
sub-domain with green being optimal and red the least optimal.

Due to the vacuum in the system, the default decomposition is unbalanced
with several MPI ranks without atoms (left). By forcing a 1x12x1
processor grid, every MPI rank does computations now, but number of
atoms per sub-domain is still uneven and the thin slice shape increases
the amount of communication between sub-domains (center left). With a
2x6x1 processor grid and shifting the sub-domain divisions, the load
imbalance is further reduced and the amount of communication required
between sub-domains is less (center right).  And using the recursive
bisectioning leads to further improved decomposition (right).


Communication
^^^^^^^^^^^^^

Following the partitioning scheme in use all per-atom data is
distributed across the MPI processes, which allows LAMMPS to handle very
large systems provided it uses a correspondingly large number of MPI
processes.  Since The per-atom data (atom IDs, positions, velocities,
types, etc.)  To be able to compute the short-range interactions MPI
processes need not only access to data of atoms they "own" but also
information about atoms from neighboring sub-domains, in LAMMPS referred
to as "ghost" atoms.  These are copies of atoms storing required
per-atom data for up to the communication cutoff distance. The green
dashed-line boxes in the :ref:`domain-decomposition` figure illustrate
the extended ghost-atom sub-domain for one processor.

This approach is also used to implement periodic boundary
conditions: atoms that lie within the cutoff distance across a periodic
boundary are also stored as ghost atoms and taken from the periodic
replication of the sub-domain, which may be the same sub-domain, e.g. if
running in serial.  As a consequence of this, force computation in
LAMMPS is not subject to minimum image conventions and thus cutoffs may
be larger than half the simulation domain.

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

Efficient communication patterns are needed to update the "ghost" atom
data, since that needs to be done at every MD time step or minimization
step.  The diagrams of the `ghost-atom-comm` figure illustrate how ghost
atom communication is performed in two stages for a 2d simulation (three
in 3d) for both a regular and irregular partitioning of the simulation
box.  For the regular case (left) atoms are exchanged first in the
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

To compute forces efficiently, each processor creates a Verlet-style
neighbor list which enumerates all pairs of atoms *i,j* (*i* = owned,
*j* = owned or ghost) with separation less than the applicable
neighbor list cutoff distance.  In LAMMPS the neighbor lists are stored
in a multiple-page data structure; each page is a contiguous chunk of
memory which stores vectors of neighbor atoms *j* for many *i* atoms.
This allows pages to be incrementally allocated or deallocated in blocks
as needed.  Neighbor lists typically consume the most memory of any data
structure in LAMMPS.  The neighbor list is rebuilt (from scratch) once
every few timesteps, then used repeatedly each step for force or other
computations.  The neighbor cutoff distance is :math:`R_n = R_f +
\Delta_s`, where :math:`R_f` is the (largest) force cutoff defined by
the interatomic potential for computing short-range pairwise or manybody
forces and :math:`\Delta_s` is a "skin" distance that allows the list to
be used for multiple steps assuming that atoms do not move very far
between consecutive time steps.  Typically the code triggers
reneighboring when any atom has moved half the skin distance since the
last reneighboring; this and other options of the neighbor list rebuild
can be adjusted with the :doc:`neigh_modify <neigh_modify>` command.

On steps when reneighboring is performed, atoms which have moved outside
their owning processor's sub-domain are first migrated to new processors
via communication.  Periodic boundary conditions are also (only)
enforced on these steps to ensure each atom is re-assigned to the
correct processor.  After migration, the atoms owned by each processor
are stored in a contiguous vector.  Periodically each processor
spatially sorts owned atoms within its vector to reorder it for improved
cache efficiency in force computations and neighbor list building.  For
that atoms are spatially binned and then reordered so that atoms in the
same bin are adjacent in the vector.  Atom sorting can be disabled or
its settings modified with the :doc:`atom_modify <atom_modify>` command.

.. _neighbor-stencil:
.. figure:: img/neigh-stencil.png
   :align: center

   neighbor list stencils

   A 2d simulation sub-domain (thick black line) and the corresponding
   ghost atom cutoff region (dashed blue line) for both orthogonal
   (left) and triclinic (right) domains.  A regular grid of neighbor
   bins (thin lines) overlays the entire simulation domain and need not
   align with sub-domain boundaries; only the portion overlapping the
   augmented sub-domain is shown.  In the triclinic case it overlaps the
   bounding box of the tilted rectangle.  The blue- and red-shaded bins
   represent a stencil of bins searched to find neighbors of a particular
   atom (black dot).

To build a local neighbor list in linear time, the simulation domain is
overlaid (conceptually) with a regular 3d (or 2d) grid of neighbor bins,
as shown in the :ref:`neighbor-stencil` figure for 2d models and a
single MPI processor's sub-domain.  Each processor stores a set of
neighbor bins which overlap its sub-domain extended by the neighbor
cutoff distance :math:`R_n`.  As illustrated, the bins need not align
with processor boundaries; an integer number in each dimension is fit to
the size of the entire simulation box.

Most often LAMMPS builds what it calls a "half" neighbor list where
each *i,j* neighbor pair is stored only once, with either atom *i* or
*j* as the central atom.  The build can be done efficiently by using a
pre-computed "stencil" of bins around a central origin bin which
contains the atom whose neighbors are being searched for.  A stencil
is simply a list of integer offsets in *x,y,z* of nearby bins
surrounding the origin bin which are close enough to contain any
neighbor atom *j* within a distance :math:`R_n` from any atom *i* in the
origin bin.  Note that for a half neighbor list, the stencil can be
asymmetric since each atom only need store half its nearby neighbors.

These stencils are illustrated in the figure for a half list and a bin
size of :math:`\frac{1}{2} R_n`.  There are 13 red+blue stencil bins in
2d (for the orthogonal case, 15 for triclinic).  In 3d there would be
63, 13 in the plane of bins that contain the origin bin and 25 in each
of the two planes above it in the *z* direction (75 for triclinic).  The
reason the triclinic stencil has extra bins is because the bins tile the
bounding box of the entire triclinic domain and thus are not periodic
with respect to the simulation box itself.  The stencil and logic for
determining which *i,j* pairs to include in the neighbor list are
altered slightly to account for this.

To build a neighbor list, a processor first loops over its "owned" plus
"ghost" atoms and assigns each to a neighbor bin.  This uses an integer
vector to create a linked list of atom indices within each bin.  It then
performs a triply-nested loop over its owned atoms *i*, the stencil of
bins surrounding atom *i*'s bin, and the *j* atoms in each stencil bin
(including ghost atoms).  If the distance :math:`r_{ij} < R_n`, then
atom *j* is added to the vector of atom *i*'s neighbors.

Here are additional details about neighbor list build options LAMMPS
supports:

- The choice of bin size is an option; a size half of :math:`R_n` has
  been found to be optimal for many typical cases.  Smaller bins incur
  additional overhead to loop over; larger bins require more distance
  calculations.  Note that for smaller bin sizes, the 2d stencil in the
  figure would be more semi-circular in shape (hemispherical in 3d),
  with bins near the corners of the square eliminated due to their
  distance from the origin bin.

- Depending on the interatomic potential(s) and other commands used in
  an input script, multiple neighbor lists and stencils with different
  attributes may be needed.  This includes lists with different cutoff
  distances, e.g. for force computation versus occasional diagnostic
  computations such as a radial distribution function, or for the
  r-RESPA time integrator which can partition pairwise forces by
  distance into subsets computed at different time intervals.  It
  includes "full" lists (as opposed to half lists) where each *i,j* pair
  appears twice, stored once with *i* and *j*, and which use a larger
  symmetric stencil.  It also includes lists with partial enumeration of
  ghost atom neighbors.  The full and ghost-atom lists are used by
  various manybody interatomic potentials.  Lists may also use different
  criteria for inclusion of a pair interaction.  Typically this simply
  depends only on the distance between two atoms and the cutoff
  distance.  But for finite-size coarse-grained particles with
  individual diameters (e.g. polydisperse granular particles), it can
  also depend on the diameters of the two particles.

- When using :doc:`pair style hybrid <pair_hybrid>` multiple sub-lists
  of the master neighbor list for the full system need to be generated,
  one for each sub-style, which contains only the *i,j* pairs needed to
  compute interactions between subsets of atoms for the corresponding
  potential.  This means not all *i* or *j* atoms owned by a processor
  are included in a particular sub-list.

- Some models use different cutoff lengths for pairwise interactions
  between different kinds of particles which are stored in a single
  neighbor list.  One example is a solvated colloidal system with large
  colloidal particles where colloid/colloid, colloid/solvent, and
  solvent/solvent interaction cutoffs can be dramatically different.
  Another is a model of polydisperse finite-size granular particles;
  pairs of particles interact only when they are in contact with each
  other.  Mixtures with particle size ratios as high as 10-100x may be
  used to model realistic systems.  Efficient neighbor list building
  algorithms for these kinds of systems are available in LAMMPS.  They
  include a method which uses different stencils for different cutoff
  lengths and trims the stencil to only include bins that straddle the
  cutoff sphere surface.  More recently a method which uses both
  multiple stencils and multiple bin sizes was developed; it builds
  neighbor lists efficiently for systems with particles of any size
  ratio, though other considerations (timestep size, force computations)
  may limit the ability to model systems with huge polydispersity.

- For small and sparse systems and as a fallback method, LAMMPS also
  supports neighbor list construction without binning by using a full
  :math:`O(N^2)` loop over all *i,j* atom pairs in a sub-domain when
  using the :doc:`neighbor nsq <neighbor>` command.


Long-range interactions
^^^^^^^^^^^^^^^^^^^^^^^
