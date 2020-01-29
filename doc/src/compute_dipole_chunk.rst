.. index:: compute dipole/chunk

compute dipole/chunk command
============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID dipole/chunk chunkID charge-correction

* ID, group-ID are documented in :doc:`compute <compute>` command
* dipole/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* charge-correction = *mass* or *geometry*\ , use COM or geometric center for charged chunk correction (optional)

Examples
""""""""


.. parsed-literal::

   compute 1 fluid dipole/chunk molchunk
   compute dw water dipole/chunk 1 geometry

Description
"""""""""""

Define a computation that calculates the dipole vector and total dipole
for multiple chunks of atoms.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>`
doc pages for details of how chunks can be defined and examples of how
they can be used to measure properties of a system.

This compute calculates the x,y,z coordinates of the dipole vector
and the total dipole moment for each chunk, which includes all effects
due to atoms passing through periodic boundaries. For chunks with a net
charge the resulting dipole is made position independent by subtracting
the position vector of the center of mass or geometric center times the
net charge from the computed dipole vector.

Note that only atoms in the specified group contribute to the
calculation.  The :doc:`compute chunk/atom <compute_chunk_atom>` command
defines its own group; atoms will have a chunk ID = 0 if they are not
in that group, signifying they are not assigned to a chunk, and will
thus also not contribute to this calculation.  You can specify the
"all" group for this command if you simply want to include atoms with
non-zero chunk IDs.

.. note::

   The coordinates of an atom contribute to the chunk's
   dipole in "unwrapped" form, by using the image flags
   associated with each atom.  See the :doc:`dump custom <dump>` command
   for a discussion of "unwrapped" coordinates.  See the Atoms section of
   the :doc:`read_data <read_data>` command for a discussion of image flags
   and how they are set for each atom.  You can reset the image flags
   (e.g. to 0) before invoking this compute by using the :doc:`set image <set>` command.

The simplest way to output the results of the compute com/chunk
calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>`
command, for example:


.. parsed-literal::

   compute cc1 all chunk/atom molecule
   compute myChunk all dipole/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk[\*] file tmp.out mode vector

**Output info:**

This compute calculates a global array where the number of rows = the
number of chunks *Nchunk* as calculated by the specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  The number of columns =
4 for the x,y,z dipole vector components and the total dipole of each
chunk. These values can be accessed by any command that uses global
array values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The array values are "intensive".  The array values will be in
dipole units, i.e. charge units times distance :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute com/chunk <compute_com_chunk>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
