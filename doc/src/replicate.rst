.. index:: replicate

replicate command
=================

Syntax
""""""

.. code-block:: LAMMPS

   replicate nx ny nz keyword ...

nx,ny,nz = replication factors in each dimension

* zero or more keywords may be appended
* keyword = *bbox* or *bond/periodic*

  .. parsed-literal::

       *bbox* = use a bounding-box algorithm which is faster for large proc counts
       *bond/periodic* = use an algorithm that correctly replicates periodic bond loops

Examples
""""""""

For examples of replicating simple linear polymer chains (periodic or
non-periodic) or periodic carbon nanotubes, see examples/replicate.

.. code-block:: LAMMPS

   replicate 2 3 2
   replicate 2 3 2 bbox
   replicate 2 3 2 bond/periodic

Description
"""""""""""

Replicate the current system one or more times in each dimension.  For
example, replication factors of 2,2,2 will create a simulation with 8x
as many atoms by doubling the size of the simulation box in each
dimension.  A replication factor of 1 leaves the simulation domain
unchanged in that dimension.

When the new simulation box is created it is partitioned into a
regular 3d grid of rectangular bricks, one per processor, based on the
number of processors being used and the settings of the
:doc:`processors <processors>` command.  The partitioning can be
changed by subsequent :doc:`balance <balance>` or :doc:`fix balance
<fix_balance>` commands.

All properties of each atom are replicated (except per-atom fix data,
see the Restrictions section below).  This includes their velocities,
which may or may not be desirable.  New atom IDs are assigned to new
atoms, as are new molecule IDs.  Bonds and other topology interactions
are created between pairs of new atoms as well as between old and new
atoms.

.. note::

   The bond discussion which follows only refers to models with
   permanent covalent bonds typically defined in LAMMPS via a data
   file.  It is not relevant to systems modeled with many-body
   potentials which can define bonds on-the-fly, based on the current
   positions of nearby atoms, e.g. models using the :doc:`AIREBO
   <pair_airebo>` or :doc:`ReaxFF <pair_reaxff>` potentials.

If the *bond/periodic* keyword is not specified, bond replication is
done by using the image flag for each atom to "unwrap" it out of the
periodic box before replicating it.  After replication is performed,
atoms outside the new periodic box are wrapped back into it.  This
assigns correct images flags to all atoms in the system.  For this to
work, all original atoms in the original simulation box must have
consistent image flags.  This means that if two atoms have a bond
between them which crosses a periodic boundary, their respective image
flags will differ by 1 in that dimension.

Image flag consistency is not possible if a system has a periodic bond
loop, meaning there is a chain of bonds which crosses an entire
dimension and re-connects to itself across a periodic boundary.  In
this case you MUST use the *bond/periodic* keyword to correctly
replicate the system.  This option zeroes the image flags for all
atoms and uses a different algorithm to find new (nearby) bond
neighbors in the replicated system.  In the final replicated system
all image flags are zero (in each dimension).

.. note::

   LAMMPS does not check for image flag consistency before performing
   the replication (it does issue a warning about this before a
   simulation is run).  If the original image flags are inconsistent,
   the replicated system will also have inconsistent image flags, but
   will otherwise be correctly replicated.  This is NOT the case if
   there is a periodic bond loop.  See the next note.

.. note::

   LAMMPS does not check for periodic bond loops.  If you use the
   *bond/periodic* keyword for a system without periodic bond loops,
   the system will be correctly replicated, but image flag information
   will be lost (which may or may not be important to your model).  If
   you do not use the *bond/periodic* keyword for a system with
   periodic bond loops, the replicated system will have invalid bonds
   (typically very long), resulting in bad dynamics.

If possible, the *bbox* keyword should be used when running on a large
number of processors, as it can result in a substantial speed-up for
the replication operation.  It uses a bounding box to only check atoms
in replicas that overlap with each processor's new subdomain when
assigning atoms to processors.  It also preserves image flag
information.  The only drawback to the *bbox* option is that it
requires a temporary use of more memory.  Each processor must be able
to store all atoms (and their per-atom data) in the original system,
before it is replicated.

.. note::

  The algorithm used by the *bond/periodic* keyword builds on the
  algorithm used by the *bbox* keyword and thus has the same memory
  requirements.  If you specify only the *bond/peridoic* keyword it
  will internally set the *bbox* keyword as well.

----------

Restrictions
""""""""""""

A 2d simulation cannot be replicated in the z dimension.

If a simulation is non-periodic in a dimension, care should be used
when replicating it in that dimension, as it may generate atoms nearly
on top of each other.

If the current simulation was read in from a restart file (before a
run is performed), there must not be any fix information stored in the
file for individual atoms.  Similarly, no fixes can be defined at the
time the replicate command is used that require vectors of atom
information to be stored.  This is because the replicate command does
not know how to replicate that information for new atoms it creates.

To work around this restriction two options are possible.  (1) Fixes
which use the stored data in the restart file can be defined before
replication and then deleted via the :doc:`unfix <unfix>` command and
re-defined after it.  Or (2) the restart file can be converted to a
data file (which deletes the stored fix information) and fixes defined
after the replicate command.  In both these scenarios, the per-atom
fix information in the restart file is lost.

Related commands
""""""""""""""""

none

Default
"""""""

No settings for using the *bbox* or *bond/periodic* algorithms.
