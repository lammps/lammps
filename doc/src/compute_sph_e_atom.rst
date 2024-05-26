.. index:: compute sph/e/atom

compute sph/e/atom command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID sph/e/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* sph/e/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all sph/e/atom

Description
"""""""""""

Define a computation that calculates the per-atom internal energy
for each atom in a group.

The internal energy is the energy associated with the internal degrees
of freedom of an SPH particle, i.e. a Smooth-Particle Hydrodynamics
particle.

See `this PDF guide <PDF/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The value of the internal energy will be 0.0 for atoms not in the
specified compute group.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the SPH package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

Default
"""""""

none
