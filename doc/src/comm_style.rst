.. index:: comm_style

comm_style command
==================

Syntax
""""""

.. code-block:: LAMMPS

   comm_style style

* style = *brick* or *brick/direct* or *tiled*

Examples
""""""""

.. code-block:: LAMMPS

   comm_style brick
   comm_style brick/direct
   comm_style tiled

Description
"""""""""""

This command sets the style of inter-processor communication of atom
information that occurs each timestep as coordinates and other
properties are exchanged between neighboring processors and stored as
properties of ghost atoms.

For the default *brick* style, the domain decomposition used by LAMMPS
to partition the simulation box must be a regular 3d grid of bricks,
one per processor.  Each processor communicates with its 6 Cartesian
neighbors in the grid to acquire information for nearby atoms.

.. versionadded:: TBD

The *brick/direct* style uses the same regular 3d grid of bricks as the
*brick* style, but acquires ghost atoms differently.  Instead of the
6-way exchange above, in which atoms are relayed through intermediate
processors and each stage must finish before the next begins, every
processor exchanges atoms directly with each of the nearby processors
that owns atoms within its ghost cutoff.  All of those exchanges are
posted at once, so no processor waits on a relay of messages it is not
part of.

This is most useful when the ghost cutoff is large compared to the size
of a subdomain, which happens when a simulation is run on many
processors, or with a long cutoff, or both.  In that regime the *brick*
style relays atoms through several stages in each dimension, and the
number of stages grows with the cutoff, while *brick/direct* always
communicates in a single stage.  For a small number of processors, where
each subdomain is large compared to the cutoff, *brick* exchanges fewer
messages and is usually faster, so *brick/direct* is not a good default.

Because the direct exchange is built from the same regular grid, the
*brick/direct* style requires a uniform processor grid and cannot be
combined with the :doc:`balance <balance>` or
:doc:`fix balance <fix_balance>` commands, which make the decomposition
non-uniform.  It also does not support the *multi* mode or the *group*
keyword of :doc:`comm_modify <comm_modify>`, and it does not yet
implement the communication some bond styles require.  LAMMPS will stop
with an error in each of those cases.

For the *tiled* style, a more general domain decomposition can be
used, as triggered by the :doc:`balance <balance>` or :doc:`fix balance <fix_balance>` commands.  The simulation box can be
partitioned into non-overlapping rectangular-shaped "tiles" or varying
sizes and shapes.  Again there is one tile per processor.  To acquire
information for nearby atoms, communication must now be done with a
more complex pattern of neighboring processors.

Note that this command does not actually define a partitioning of the
simulation box (a domain decomposition), rather it determines what
kinds of decompositions are allowed and the pattern of communication
used to enable the decomposition.  A decomposition is created when the
simulation box is first created, via the :doc:`create_box <create_box>`
or :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.  For the *brick*, *brick/direct* and *tiled*
styles, the initial decomposition will be the same, as described by
:doc:`create_box <create_box>` and :doc:`processors <processors>`
commands.  The decomposition can be changed via the
:doc:`balance <balance>` or :doc:`fix balance <fix_balance>` commands.

Restrictions
""""""""""""

The *brick/direct* style requires a uniform processor grid, and does not
support :doc:`comm_modify <comm_modify>` *mode multi* or the *group*
keyword, variable-size reverse communication from a fix, or the
communication used by some bond styles.

Related commands
""""""""""""""""

:doc:`comm_modify <comm_modify>`, :doc:`processors <processors>`,
:doc:`balance <balance>`, :doc:`fix balance <fix_balance>`

Default
"""""""

The default style is *brick*.
