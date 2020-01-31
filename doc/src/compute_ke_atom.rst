.. index:: compute ke/atom

compute ke/atom command
=======================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID ke/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* ke/atom = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all ke/atom

Description
"""""""""""

Define a computation that calculates the per-atom translational
kinetic energy for each atom in a group.

The kinetic energy is simply 1/2 m v\^2, where m is the mass and v is
the velocity of each atom.

The value of the kinetic energy will be 0.0 for atoms not in the
specified compute group.

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


