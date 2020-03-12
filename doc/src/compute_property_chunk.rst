.. index:: compute property/chunk

compute property/chunk command
==============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID property/chunk chunkID input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* property/chunk = style name of this compute command
* input = one or more attributes

  .. parsed-literal::

       attributes = count, id, coord1, coord2, coord3
         count = # of atoms in chunk
         id = original chunk IDs before compression by :doc:`compute chunk/atom <compute_chunk_atom>`
         coord123 = coordinates for spatial bins calculated by :doc:`compute chunk/atom <compute_chunk_atom>`

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all property/chunk count
   compute 1 all property/chunk ID coord1

Description
"""""""""""

Define a computation that stores the specified attributes of chunks of
atoms.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>`
doc pages for details of how chunks can be defined and examples of how
they can be used to measure properties of a system.

This compute calculates and stores the specified attributes of chunks
as global data so they can be accessed by other :doc:`output commands <Howto_output>` and used in conjunction with other
commands that generate per-chunk data, such as :doc:`compute com/chunk <compute_com_chunk>` or :doc:`compute msd/chunk <compute_msd_chunk>`.

Note that only atoms in the specified group contribute to the
calculation of the *count* attribute.  The :doc:`compute chunk/atom <compute_chunk_atom>` command defines its own group;
atoms will have a chunk ID = 0 if they are not in that group,
signifying they are not assigned to a chunk, and will thus also not
contribute to this calculation.  You can specify the "all" group for
this command if you simply want to include atoms with non-zero chunk
IDs.

The *count* attribute is the number of atoms in the chunk.

The *id* attribute stores the original chunk ID for each chunk.  It
can only be used if the *compress* keyword was set to *yes* for the
:doc:`compute chunk/atom <compute_chunk_atom>` command referenced by
chunkID.  This means that the original chunk IDs (e.g. molecule IDs)
will have been compressed to remove chunk IDs with no atoms assigned
to them.  Thus a compressed chunk ID of 3 may correspond to an original
chunk ID (molecule ID in this case) of 415.  The *id* attribute will
then be 415 for the 3rd chunk.

The *coordN* attributes can only be used if a *binning* style was used
in the :doc:`compute chunk/atom <compute_chunk_atom>` command referenced
by chunkID.  For *bin/1d*\ , *bin/2d*\ , and *bin/3d* styles the attribute
is the center point of the bin in the corresponding dimension.  Style
*bin/1d* only defines a *coord1* attribute.  Style *bin/2d* adds a
*coord2* attribute.  Style *bin/3d* adds a *coord3* attribute.

Note that if the value of the *units* keyword used in the :doc:`compute chunk/atom command <compute_chunk_atom>` is *box* or *lattice*\ , the
*coordN* attributes will be in distance :doc:`units <units>`.  If the
value of the *units* keyword is *reduced*\ , the *coordN* attributes
will be in unitless reduced units (0-1).

The simplest way to output the results of the compute property/chunk
calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>`
command, for example:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   compute myChunk1 all property/chunk cc1 count
   compute myChunk2 all com/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk1 c_myChunk2[*] file tmp.out mode vector

**Output info:**

This compute calculates a global vector or global array depending on
the number of input values.  The length of the vector or number of
rows in the array is the number of chunks.

This compute calculates a global vector or global array where the
number of rows = the number of chunks *Nchunk* as calculated by the
specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  If a
single input is specified, a global vector is produced.  If two or
more inputs are specified, a global array is produced where the number
of columns = the number of inputs.  The vector or array can be
accessed by any command that uses global values from a compute as
input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

The vector or array values are "intensive".  The values will be
unitless or in the units discussed above.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix ave/chunk <fix_ave_chunk>`

**Default:** none
