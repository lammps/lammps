.. index:: compute sph/rho/atom

compute sph/rho/atom command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID sph/rho/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* sph/rho/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all sph/rho/atom

Description
"""""""""""

Define a computation that calculates the per-atom SPH density for each
atom in a group, i.e. a Smooth-Particle Hydrodynamics density.

The SPH density is the mass density of an SPH particle, calculated by
kernel function interpolation using "pair style sph/rhosum".

See `this PDF guide <PDF/SPH_LAMMPS_userguide.pdf>`_ to using SPH in
LAMMPS.

The value of the SPH density will be 0.0 for atoms not in the
specified compute group.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-atom vector values will be in mass/volume :doc:`units <units>`.

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
