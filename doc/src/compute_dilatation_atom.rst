.. index:: compute dilatation/atom

compute dilatation/atom command
===============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID dilatation/atom

* ID, group-ID are documented in compute command
* dilation/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all dilatation/atom

Description
"""""""""""

Define a computation that calculates the per-atom dilatation for each
atom in a group.  This is a quantity relevant for :doc:`Peridynamics models <pair_peri>`.  See `this document <PDF/PDLammps_overview.pdf>`_
for an overview of LAMMPS commands for Peridynamics modeling.

For small deformation, dilatation of is the measure of the volumetric
strain.

The dilatation "theta" for each peridynamic particle I is calculated
as a sum over its neighbors with unbroken bonds, where the
contribution of the IJ pair is a function of the change in bond length
(versus the initial length in the reference state), the volume
fraction of the particles and an influence function.  See the
`PDLAMMPS user guide <http://www.sandia.gov/~mlparks/papers/PDLAMMPS.pdf>`_ for a formal
definition of dilatation.

This command can only be used with a subset of the Peridynamic :doc:`pair styles <pair_peri>`: peri/lps, peri/ves and peri/eps.

The dilatation value will be 0.0 for atoms not in the specified
compute group.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The per-atom vector values are unitless numbers (theta) >= 0.0.

Restrictions
""""""""""""

This compute is part of the PERI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute damage/atom <compute_damage_atom>`,
:doc:`compute plasticity/atom <compute_plasticity_atom>`

Default
"""""""

none
