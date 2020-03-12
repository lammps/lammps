.. index:: fix spring/chunk

fix spring/chunk command
========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID spring/chunk K chunkID comID

* ID, group-ID are documented in :doc:`fix <fix>` command
* spring/chunk = style name of this fix command
* K = spring constant for each chunk (force/distance units)
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* comID = ID of :doc:`compute com/chunk <compute_com_chunk>` command

Examples
""""""""

.. code-block:: LAMMPS

   fix restrain all spring/chunk 100 chunkID comID

Description
"""""""""""

Apply a spring force to the center-of-mass (COM) of chunks of atoms as
defined by the :doc:`compute chunk/atom <compute_chunk_atom>` command.
Chunks can be molecules or spatial bins or other groupings of atoms.
This is a way of tethering each chunk to its initial COM coordinates.

The *chunkID* is the ID of a compute chunk/atom command defined in the
input script.  It is used to define the chunks.  The *comID* is the ID
of a compute com/chunk command defined in the input script.  It is
used to compute the COMs of each chunk.

At the beginning of the first :doc:`run <run>` or
:doc:`minimize <minimize>` command after this fix is defined, the
initial COM of each chunk is calculated and stored as R0m, where M is
the chunk number.  Thereafter, at every timestep (or minimization
iteration), the current COM of each chunk is calculated as Rm.  A
restoring force of magnitude K (Rm - R0m) Mi / Mm is applied to each
atom in each chunk where *K* is the specified spring constant, Mi is
the mass of the atom, and Mm is the total mass of all atoms in the
chunk.  Note that *K* thus represents the spring constant for the
total force on each chunk of atoms, not for a spring applied to each
atom.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the energy stored in all the springs to the system's potential
energy as part of :doc:`thermodynamic output <thermo_style>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the energy of all
the springs, i.e. 0.5 \* K \* r\^2 per-spring.

The scalar value calculated by this fix is "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the spring energies to be included in the total
   potential energy of the system (the quantity being minimized), you
   MUST enable the :doc:`fix_modify <fix_modify>` *energy* option for this
   fix.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix spring <fix_spring>`, :doc:`fix spring/self <fix_spring_self>`,
:doc:`fix spring/rg <fix_spring_rg>`

**Default:** none
