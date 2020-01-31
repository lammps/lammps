.. index:: compute erotate/sphere/atom

compute erotate/sphere/atom command
===================================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID erotate/sphere/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* erotate/sphere/atom = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all erotate/sphere/atom

Description
"""""""""""

Define a computation that calculates the rotational kinetic energy for
each particle in a group.

The rotational energy is computed as 1/2 I w\^2, where I is the moment
of inertia for a sphere and w is the particle's angular velocity.

.. note::

   For :doc:`2d models <dimension>`, particles are treated as
   spheres, not disks, meaning their moment of inertia will be the same
   as in 3d.

The value of the rotational kinetic energy will be 0.0 for atoms not
in the specified compute group or for point particles with a radius =
0.0.

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

**Default:** none


