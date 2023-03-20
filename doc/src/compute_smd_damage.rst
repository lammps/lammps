.. index:: compute smd/damage

compute smd/damage command
==========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID smd/damage

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/damage = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/damage

Description
"""""""""""

Define a computation that calculates the damage status of SPH particles
according to the damage model which is defined via the SMD SPH pair styles, e.g., the maximum plastic strain failure criterion.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to use Smooth Mach Dynamics in LAMMPS.

**Output Info:**

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-particle values are dimensionless an in the range of zero to one.

Restrictions
""""""""""""

This compute is part of the MACHDYN package.  It is only enabled if
LAMMPS was built with that package.  See the "Build

Related commands
""""""""""""""""

:doc:`smd/plastic_strain <compute_smd_plastic_strain>`, :doc:`smd/tlsph_stress <compute_smd_tlsph_stress>`

Default
"""""""

none
