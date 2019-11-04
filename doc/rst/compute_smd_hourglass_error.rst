.. index:: compute smd/hourglass/error

compute smd/hourglass/error command
===================================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID smd/hourglass/error

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/hourglass/error = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all smd/hourglass/error

Description
"""""""""""

Define a computation which outputs the error of the approximated
relative separation with respect to the actual relative separation of
the particles i and j. Ideally, if the deformation gradient is exact,
and there exists a unique mapping between all particles' positions
within the neighborhood of the central node and the deformation
gradient, the approximated relative separation will coincide with the
actual relative separation of the particles i and j in the deformed
configuration.  This compute is only really useful for debugging the
hourglass control mechanism which is part of the Total-Lagrangian SPH
pair style.

See `this PDF guide <PDF/SMD_LAMMPS_userguide.pdf>`_ to use Smooth
Mach Dynamics in LAMMPS.

**Output Info:**

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-particle vector values will are dimensionless. See
:doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This quantity will be computed only for particles which interact with
tlsph pair style.

**Related Commands:**

:doc:`smd/tlsph\_defgrad <compute_smd_tlsph_defgrad>`

Default
"""""""


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
