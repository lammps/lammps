.. index:: dimension

dimension command
=================

Syntax
""""""


.. parsed-literal::

   dimension N

* N = 2 or 3

Examples
""""""""


.. parsed-literal::

   dimension 2

Description
"""""""""""

Set the dimensionality of the simulation.  By default LAMMPS runs 3d
simulations.  To run a 2d simulation, this command should be used
prior to setting up a simulation box via the
:doc:`create_box <create_box>` or :doc:`read_data <read_data>` commands.
Restart files also store this setting.

See the discussion on the :doc:`Howto 2d <Howto_2d>` doc page for
additional instructions on how to run 2d simulations.

.. note::

   Some models in LAMMPS treat particles as finite-size spheres or
   ellipsoids, as opposed to point particles.  In 2d, the particles will
   still be spheres or ellipsoids, not circular disks or ellipses,
   meaning their moment of inertia will be the same as in 3d.

Restrictions
""""""""""""


This command must be used before the simulation box is defined by a
:doc:`read_data <read_data>` or :doc:`create_box <create_box>` command.

Related commands
""""""""""""""""

:doc:`fix enforce2d <fix_enforce2d>`

Default
"""""""


.. parsed-literal::

   dimension 3
