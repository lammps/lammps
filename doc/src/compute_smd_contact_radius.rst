.. index:: compute smd/contact/radius

compute smd/contact/radius command
==================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/contact/radius

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/contact/radius = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/contact/radius

Description
"""""""""""

Define a computation which outputs the contact radius, i.e., the
radius used to prevent particles from penetrating each other.  The
contact radius is used only to prevent particles belonging to
different physical bodies from penetrating each other. It is used by
the contact pair styles, e.g., smd/hertz and smd/tri_surface.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

The value of the contact radius will be 0.0 for particles not in the
specified compute group.

**Output info:**

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-particle vector values will be in distance :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`dump custom <dump>` smd/hertz smd/tri_surface

**Default:** none
