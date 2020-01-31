.. index:: compute meso/e/atom

compute meso/e/atom command
===========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID meso/e/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* meso/e/atom = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all meso/e/atom

Description
"""""""""""

Define a computation that calculates the per-atom internal energy
for each atom in a group.

The internal energy is the energy associated with the internal degrees
of freedom of a mesoscopic particles, e.g. a Smooth-Particle
Hydrodynamics particle.

See `this PDF guide <USER/sph/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The value of the internal energy will be 0.0 for atoms not in the
specified compute group.

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""


This compute is part of the USER-SPH package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

**Default:** none


