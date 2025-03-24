.. index:: compute vacf/chunk

compute vacf/chunk command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID vacf/chunk chunkID

* ID, group-ID are documented in :doc:`compute <compute>` command
* vacf/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all vacf/chunk molchunk

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the velocity auto-correlation
function (VACF) for multiple chunks of atoms.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute
chunk/atom <compute_chunk_atom>` command, which assigns each atom to a
single chunk (or no chunk).  The ID for this command is specified as
chunkID.  For example, a single chunk could be the atoms in a molecule
or atoms in a spatial bin.  See the :doc:`compute chunk/atom
<compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>` doc pages for
details of how chunks can be defined and examples of how they can be
used to measure properties of a system.

Four quantities are calculated by this compute for each chunk.  The
first 3 quantities are the product of the initial center of mass
velocity (VCM) for each chunk in *x*, *y*, and *z* direction with the
current center of mass velocity in the same direction.  The fourth
component is the total VACF, i.e. the sum of the three components.

Note that only atoms in the specified group contribute to the
calculation.  The :doc:`compute chunk/atom <compute_chunk_atom>` command
defines its own group; atoms will have a chunk ID = 0 if they are not in
that group, signifying they are not assigned to a chunk, and will thus
also not contribute to this calculation.  You can specify the "all"
group for this command if you simply want to include atoms with non-zero
chunk IDs.

The integral of the VACF versus time is proportional to the diffusion
coefficient of the diffusing chunks.

.. note::

   The number of chunks *Nchunk* calculated by the
   :doc:`compute chunk/atom <compute_chunk_atom>` command must remain constant
   each time this compute is invoked, so that the dot product for each chunk
   from its original position can be computed consistently.  If *Nchunk*
   does not remain constant, an error will be generated.  If needed, you
   can enforce a constant *Nchunk* by using the *nchunk once* or *ids once*
   options when specifying the :doc:`compute chunk/atom <compute_chunk_atom>`
   command.

.. note::

   This compute stores the original center-of-mass velocities of each
   chunk.  When a VACF is calculated on a later timestep, it is assumed
   that the same atoms are assigned to the same chunk ID.  However
   LAMMPS has no simple way to ensure this is the case, though you can
   use the *ids once* option when specifying the :doc:`compute
   chunk/atom <compute_chunk_atom>` command.  Note that if this is not
   the case, the VACF calculation does not have a sensible meaning.

.. note::

   If you want the quantities calculated by this compute to be
   continuous when running from a :doc:`restart file <read_restart>`, then
   you should use the same ID for this compute, as in the original run.
   This is so that the fix this compute creates to store per-chunk
   quantities will also have the same ID, and thus be initialized
   correctly with chunk reference positions from the restart file.

The simplest way to output the results of the compute vacf/chunk
calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>`
command, for example:

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   compute myChunk all vacf/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk[*] file tmp.out mode vector

Output info
"""""""""""

This compute calculates a global array where the number of rows = the
number of chunks *Nchunk* as calculated by the specified :doc:`compute
chunk/atom <compute_chunk_atom>` command.  The number of columns = 4 for
the *x*, *y*, *z*, component and the total VACF.  These values can be
accessed by any command that uses global array values from a compute as
input.  See the :doc:`Howto output <Howto_output>` page for an overview
of LAMMPS output options.

The array values are "intensive".  The array values will be in
distance\ :math:`^2` divided by time\ :math:`^2` :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute vacf <compute_vacf>`, :doc:`compute msd/chunk <compute_msd_chunk>`

Default
"""""""

none
