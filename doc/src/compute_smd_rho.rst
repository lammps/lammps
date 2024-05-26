.. index:: compute smd/rho

compute smd/rho command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID smd/rho

* ID, group-ID are documented in :doc:`compute <compute>` command
* smd/rho = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all smd/rho

Description
"""""""""""

Define a computation that calculates the per-particle mass density.
The mass density is the mass of a particle which is constant during
the course of a simulation, divided by its volume, which can change
due to mechanical deformation.

See `this PDF guide <PDF/MACHDYN_LAMMPS_userguide.pdf>`_ to use Smooth
Mach Dynamics in LAMMPS.

Output info
"""""""""""

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.
See the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-particle values will be in :doc:`units <units>` of mass over volume.

Restrictions
""""""""""""

This compute is part of the MACHDYN package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute smd/vol <compute_smd_vol>`

Default
"""""""

none
