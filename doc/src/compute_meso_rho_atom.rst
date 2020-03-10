.. index:: compute meso/rho/atom

compute meso/rho/atom command
=============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID meso/rho/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* meso/rho/atom = style name of this compute command

Examples
""""""""

.. parsed-literal::

   compute 1 all meso/rho/atom

Description
"""""""""""

Define a computation that calculates the per-atom mesoscopic density
for each atom in a group.

The mesoscopic density is the mass density of a mesoscopic particle,
calculated by kernel function interpolation using "pair style
sph/rhosum".

See `this PDF guide <USER/sph/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The value of the mesoscopic density will be 0.0 for atoms not in the
specified compute group.

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values will be in mass/volume :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the USER-SPH package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`dump custom <dump>`

**Default:** none
