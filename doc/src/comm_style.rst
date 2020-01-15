.. index:: comm\_style

comm\_style command
===================

Syntax
""""""


.. parsed-literal::

   comm_style style

* style = *brick* or *tiled*

Examples
""""""""


.. parsed-literal::

   comm_style brick
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
simulation box is first created, via the :doc:`create\_box <create_box>`
or :doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands.  For both the *brick* and *tiled* styles, the initial
decomposition will be the same, as described by
:doc:`create\_box <create_box>` and :doc:`processors <processors>`
commands.  The decomposition can be changed via the
:doc:`balance <balance>` or :doc:`fix balance <fix_balance>` commands.

Restrictions
""""""""""""


Communication style *tiled* cannot be used with *triclinic* simulation
cells.

Related commands
""""""""""""""""

:doc:`comm\_modify <comm_modify>`, :doc:`processors <processors>`,
:doc:`balance <balance>`, :doc:`fix balance <fix_balance>`

Default
"""""""

The default style is brick.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
