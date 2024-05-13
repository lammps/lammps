.. index:: compute smd/vol

compute smd/vol command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID smd/vol

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/vol = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/vol

Description
"""""""""""

Define a computation that provides the per-particle volume and the sum
of the per-particle volumes of the group for which the fix is defined.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-particle vector values will be given in :doc:`units <units>` of
volume.

Additionally, the compute returns a scalar, which is the sum of the
per-particle volumes of the group for which the fix is defined.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute smd/rho <compute_smd_rho>`

Default
"""""""

none
