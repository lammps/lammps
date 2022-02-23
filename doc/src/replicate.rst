.. index:: replicate

replicate command
=================

Syntax
""""""

.. parsed-literal::

   replicate nx ny nz *keyword*

nx,ny,nz = replication factors in each dimension

* optional *keyword* = *bbox* or *bondlist*

  .. parsed-literal::

       *bbox* = only check atoms in replicas that overlap with a processor's sub-domain
       *bondlist* = use a generalized algorithm that correctly replicates periodic loops

Examples
""""""""

For an example of replicating periodically looped carbon nanotubes, see examples/replicate.

.. code-block:: LAMMPS

   replicate 2 3 2

Description
"""""""""""

Replicate the current simulation one or more times in each dimension.
For example, replication factors of 2,2,2 will create a simulation
with 8x as many atoms by doubling the simulation domain in each
dimension.  A replication factor of 1 in a dimension leaves the
simulation domain unchanged.  When the new simulation box is created
it is also partitioned into a regular 3d grid of rectangular bricks,
one per processor, based on the number of processors being used and
the settings of the :doc:`processors <processors>` command.  The
partitioning can later be changed by the :doc:`balance <balance>` or
:doc:`fix balance <fix_balance>` commands.

All properties of the atoms are replicated, including their
velocities, which may or may not be desirable.  New atom IDs are
assigned to new atoms, as are molecule IDs.  Bonds and other topology
interactions are created between pairs of new atoms as well as between
old and new atoms.  This is done by using the image flag for each atom
to "unwrap" it out of the periodic box before replicating it.

This means that any molecular bond you specify in the original data
file that crosses a periodic boundary should be between two atoms with
image flags that differ by 1.  This will allow the bond to be
unwrapped appropriately.

The optional keyword *bbox* uses a bounding box to only check atoms in
replicas that overlap with a processor's sub-domain when assigning
atoms to processors.  It typically results in a substantial speedup
when using the replicate command on a large number of processors.  It
does require temporary use of more memory, specifically that each
processor can store all atoms in the entire system before it is
replicated.

The optional keyword *bondlist* correctly treats molecules that span
the box and are bonded to themselves across a periodic boundary, by
relying on self-consistent nearest-image assumptions (rather than
using image flags).  The *bondlist* option resets all image flags to
zero.  Therefore, the *bondlist* keyword can also be used in general
for systems that may not have consistent image flags.  The *bondlist*
algorithm builds off the *bbox* algorithm, so it is fast when using a
large number of processors, but does require temporary use of more
memory.  Specifically, each processor must be able to store arrays for
all atoms in the entire system before it is replicated.

.. note::

   For systems that contain a molecule that spans the box and is
   bonded to itself across a periodic boundary (so that the molecule
   is effectively a loop), the *bondlist* keyword must be used.  A
   simple example would be a linear polymer chain that spans the
   simulation box and bonds back to itself across the periodic
   boundary.  More realistic examples would be a CNT (meant to be an
   infinitely long CNT) or a graphene sheet or a bulk periodic crystal
   where there are explicit bonds specified between near neighbors.
   (Note that this restriction only applies to systems that have
   permanent bonds as specified in the data file.  A CNT that is just
   atoms modeled with the :doc:`AIREBO potential <pair_airebo>` has no
   such permanent bonds, so it can be replicated without the
   *bondlist* keyword.)

Restrictions
""""""""""""

A 2d simulation cannot be replicated in the z dimension.

If a simulation is non-periodic in a dimension, care should be used
when replicating it in that dimension, as it may put atoms nearly on
top of each other.

If the current simulation was read in from a restart file (before a
run is performed), there must not be any fix information stored in
the file for individual atoms.  Similarly, no fixes can be defined at
the time the replicate command is used that require vectors of atom
information to be stored.  This is because the replicate command does
not know how to replicate that information for new atoms it creates.
To work around this restriction, restart files may be converted into
data files and fixes may be undefined via the :doc:`unfix <unfix>`
command before and redefined after the replicate command.

Related commands
""""""""""""""""

none


Default
"""""""

none
