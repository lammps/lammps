.. index:: compute vcm/chunk

compute vcm/chunk command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID vcm/chunk chunkID

* ID, group-ID are documented in :doc:`compute <compute>` command
* vcm/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 fluid vcm/chunk molchunk

Description
"""""""""""

Define a computation that calculates the center-of-mass velocity for
multiple chunks of atoms.

In LAMMPS, chunks are collections of atoms defined by a
:doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified as
chunkID.  For example, a single chunk could be the atoms in a molecule or atoms
in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and
:doc:`Howto chunk <Howto_chunk>` doc pages for details of how chunks can be
defined and examples of how they can be used to measure properties of a system.

This compute calculates the :math:`(x,y,z)` components of the center-of-mass
velocity for each chunk.  This is done by summing mass\*velocity for
each atom in the chunk and dividing the sum by the total mass of the
chunk.

Note that only atoms in the specified group contribute to the
calculation.  The :doc:`compute chunk/atom <compute_chunk_atom>` command
defines its own group; atoms will have a chunk ID = 0 if they are not
in that group, signifying they are not assigned to a chunk, and will
thus also not contribute to this calculation.  You can specify the
"all" group for this command if you simply want to include atoms with
non-zero chunk IDs.

The simplest way to output the results of the compute vcm/chunk
calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>`
command, for example:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   compute myChunk all vcm/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk[*] file tmp.out mode vector

Output info
"""""""""""

This compute calculates a global array where the number of rows is the
number of chunks *Nchunk* as calculated by the specified
:doc:`compute chunk/atom <compute_chunk_atom>` command.  The number of
columns is 3 for the :math:`(x,y,z)` center-of-mass velocity coordinates of
each chunk. These values can be accessed by any command that uses global array
values from a compute as input.  See the :doc:`Howto output <Howto_output>`
page for an overview of LAMMPS output options.

The array values are "intensive".  The array values will be in
velocity :doc:`units <units>`.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

none


Default
"""""""

none
