.. index:: compute displace/atom

compute displace/atom command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID displace/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* displace/atom = style name of this compute command
* zero or more keyword/arg pairs may be appended
* keyword = *refresh*

  .. parsed-literal::

       *refresh* arg = name of per-atom variable

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all displace/atom
   compute 1 all displace/atom refresh myVar

Description
"""""""""""

Define a computation that calculates the current displacement of each
atom in the group from its original (reference) coordinates, including
all effects due to atoms passing through periodic boundaries.

A vector of four quantities per atom is calculated by this compute.
The first three elements of the vector are the :math:`(dx,dy,dz)`
displacements.  The fourth component is the total displacement
(i.e., :math:`\sqrt{dx^2 + dy^2 + dz^2}`).

The displacement of an atom is from its original position at the time
the compute command was issued.  The value of the displacement will be
0.0 for atoms not in the specified compute group.

.. note::

   Initial coordinates are stored in "unwrapped" form, by using the
   image flags associated with each atom.  See the :doc:`dump custom
   <dump>` command for a discussion of "unwrapped" coordinates.  See
   the Atoms section of the :doc:`read_data <read_data>` command for a
   discussion of image flags and how they are set for each atom.  You
   can reset the image flags (e.g., to 0) before invoking this compute
   by using the :doc:`set image <set>` command.

.. note::

   If you want the quantities calculated by this compute to be
   continuous when running from a :doc:`restart file <read_restart>`, then
   you should use the same ID for this compute, as in the original run.
   This is so that the fix this compute creates to store per-atom
   quantities will also have the same ID, and thus be initialized
   correctly with time = 0 atom coordinates from the restart file.

----------

The *refresh* option can be used in conjunction with the "dump_modify
refresh" command to generate incremental dump files.

The definition and motivation of an incremental dump file is as
follows.  Instead of outputting all atoms at each snapshot (with some
associated values), you may only wish to output the subset of atoms
with a value that has changed in some way compared to the value the
last time that atom was output.  In some scenarios this can result in
a dramatically smaller dump file.  If desired, by post-processing the
sequence of snapshots, the values for all atoms at all timesteps can
be inferred.

A concrete example using this compute, is a simulation of atom
diffusion in a solid, represented as atoms on a lattice.  Diffusive
hops are rare.  Imagine that when a hop occurs an atom moves more than
a distance *Dhop*\ .  For any snapshot we only want to output atoms that
have hopped since the last snapshot.  This can be accomplished with
something like the following commands:

.. code-block:: LAMMPS

   write_dump      all custom tmp.dump id type x y z    # see comment below

   variable        Dhop equal 0.6
   variable        check atom "c_dsp[4] > v_Dhop"
   compute         dsp all displace/atom refresh check
   dump            1 all custom 100 tmp.dump id type x y z
   dump_modify     1 append yes thresh c_dsp[4] > ${Dhop} &
                   refresh c_dsp delay 100

The :doc:`dump_modify thresh <dump_modify>` command will only output
atoms that have displaced more than :math:`0.6~\AA` on each
snapshot (assuming metal units).  The dump_modify *refresh* option triggers a
call to this compute at the end of every dump.

The *refresh* argument for this compute is the ID of an
:doc:`atom-style variable <variable>` which calculates a Boolean value (0 or 1)
based on the same criterion used by dump_modify thresh.  This compute
evaluates the atom-style variable.  For each atom that returns 1 (true),
the original (reference) coordinates of the atom (stored by
this compute) are updated.

The effect of these commands is that a particular atom will only be
output in the dump file on the snapshot after it makes a diffusive
hop.  It will not be output again until it makes another hop.

Note that in the first snapshot of a subsequent run, no atoms will be
typically be output.  That is because the initial displacement for all
atoms is 0.0.  If an initial dump snapshot is desired, containing the
initial reference positions of all atoms, one way to do this is
illustrated above.  An initial write_dump command can be used before
the first run.  It will contain the positions of all the atoms,
Options in the :doc:`dump_modify <dump_modify>` command above will
append new output to that same file and delay the output until a later
timestep.  The *delay* setting avoids a second time = 0 snapshot which
would be empty.

----------

Output info
"""""""""""

This compute calculates a per-atom array with four columns, which can be
accessed by indices 1--4 by any command that uses per-atom values from
a compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

The per-atom array values will be in distance :doc:`units <units>`.

This compute supports the *refresh* option as explained above, for use
in conjunction with :doc:`dump_modify refresh <dump_modify>` to generate
incremental dump files.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute msd <compute_msd>`, :doc:`dump custom <dump>`, :doc:`fix store/state <fix_store_state>`

Default
"""""""

none
