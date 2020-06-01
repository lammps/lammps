.. index:: compute msd/chunk

compute msd/chunk command
=========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID msd/chunk chunkID

* ID, group-ID are documented in :doc:`compute <compute>` command
* msd/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all msd/chunk molchunk

Description
"""""""""""

Define a computation that calculates the mean-squared displacement
(MSD) for multiple chunks of atoms.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>`
doc pages for details of how chunks can be defined and examples of how
they can be used to measure properties of a system.

Four quantities are calculated by this compute for each chunk.  The
first 3 quantities are the squared dx,dy,dz displacements of the
center-of-mass.  The 4th component is the total squared displacement,
i.e. (dx\*dx + dy\*dy + dz\*dz) of the center-of-mass.  These
calculations include all effects due to atoms passing through periodic
boundaries.

Note that only atoms in the specified group contribute to the
calculation.  The :doc:`compute chunk/atom <compute_chunk_atom>` command
defines its own group; atoms will have a chunk ID = 0 if they are not
in that group, signifying they are not assigned to a chunk, and will
thus also not contribute to this calculation.  You can specify the
"all" group for this command if you simply want to include atoms with
non-zero chunk IDs.

The slope of the mean-squared displacement (MSD) versus time is
proportional to the diffusion coefficient of the diffusing chunks.

The displacement of the center-of-mass of the chunk is from its
original center-of-mass position, calculated on the timestep this
compute command was first invoked.

.. note::

   The number of chunks *Nchunk* calculated by the :doc:`compute chunk/atom <compute_chunk_atom>` command must remain constant each
   time this compute is invoked, so that the displacement for each chunk
   from its original position can be computed consistently.  If *Nchunk*
   does not remain constant, an error will be generated.  If needed, you
   can enforce a constant *Nchunk* by using the *nchunk once* or *ids
   once* options when specifying the :doc:`compute chunk/atom <compute_chunk_atom>` command.

.. note::

   This compute stores the original position (of the
   center-of-mass) of each chunk.  When a displacement is calculated on a
   later timestep, it is assumed that the same atoms are assigned to the
   same chunk ID.  However LAMMPS has no simple way to insure this is the
   case, though you can use the *ids once* option when specifying the
   :doc:`compute chunk/atom <compute_chunk_atom>` command.  Note that if
   this is not the case, the MSD calculation does not have a sensible
   meaning.

.. note::

   The initial coordinates of the atoms in each chunk are stored in
   "unwrapped" form, by using the image flags associated with each atom.
   See the :doc:`dump custom <dump>` command for a discussion of
   "unwrapped" coordinates.  See the Atoms section of the
   :doc:`read_data <read_data>` command for a discussion of image flags and
   how they are set for each atom.  You can reset the image flags
   (e.g. to 0) before invoking this compute by using the :doc:`set image <set>` command.

.. note::

   If you want the quantities calculated by this compute to be
   continuous when running from a :doc:`restart file <read_restart>`, then
   you should use the same ID for this compute, as in the original run.
   This is so that the fix this compute creates to store per-chunk
   quantities will also have the same ID, and thus be initialized
   correctly with chunk reference positions from the restart file.

The simplest way to output the results of the compute msd/chunk
calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>`
command, for example:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   compute myChunk all msd/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk[*] file tmp.out mode vector

**Output info:**

This compute calculates a global array where the number of rows = the
number of chunks *Nchunk* as calculated by the specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  The number of columns =
4 for dx,dy,dz and the total displacement.  These values can be
accessed by any command that uses global array values from a compute
as input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

The array values are "intensive".  The array values will be in
distance\^2 :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute msd <compute_msd>`

**Default:** none
