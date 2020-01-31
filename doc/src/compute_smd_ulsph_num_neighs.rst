.. index:: compute smd/ulsph/num/neighs

compute smd/ulsph/num/neighs command
====================================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID smd/ulsph/num/neighs

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/ulsph/num/neighs = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all smd/ulsph/num/neighs

Description
"""""""""""

Define a computation that returns the number of neighbor particles
inside of the smoothing kernel radius for particles interacting via
the updated Lagrangian SPH pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to using Smooth
Mach Dynamics in LAMMPS.

**Output info:**

This compute returns a per-particle vector, which can be accessed by
any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-particle values will be given dimensionless, see :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.  This compute can
only be used for particles which interact with the updated Lagrangian
SPH pair style.

Related commands
""""""""""""""""

:doc:`compute smd/tlsph/num/neighs <compute_smd_tlsph_num_neighs>`

**Default:** none
