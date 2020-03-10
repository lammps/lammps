.. index:: compute smd/internal/energy

compute smd/internal/energy command
===================================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/internal/energy

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/smd/internal/energy = style name of this compute command

Examples
""""""""

.. parsed-literal::

   compute 1 all smd/internal/energy

Description
"""""""""""

Define a computation which outputs the per-particle enthalpy, i.e.,
the sum of potential energy and heat.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to use Smooth
Mach Dynamics in LAMMPS.

**Output Info:**

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-particle vector values will be given in :doc:`units <units>` of energy.

Restrictions
""""""""""""

This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info. This compute can
only be used for particles which interact via the updated Lagrangian
or total Lagrangian SPH pair styles.

**Related Commands:**

Default
"""""""
