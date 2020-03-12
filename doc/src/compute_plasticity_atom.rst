.. index:: compute plasticity/atom

compute plasticity/atom command
===============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID plasticity/atom

* ID, group-ID are documented in compute command
* plasticity/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all plasticity/atom

Description
"""""""""""

Define a computation that calculates the per-atom plasticity for each
atom in a group.  This is a quantity relevant for :doc:`Peridynamics models <pair_peri>`.  See `this document <PDF/PDLammps_overview.pdf>`_
for an overview of LAMMPS commands for Peridynamics modeling.

The plasticity for a Peridynamic particle is the so-called consistency
parameter (lambda).  For elastic deformation lambda = 0, otherwise
lambda > 0 for plastic deformation.  For details, see
:ref:`(Mitchell) <Mitchell>` and the PDF doc included in the LAMMPS
distribution in `doc/PDF/PDLammps\_EPS.pdf <PDF/PDLammps_EPS.pdf>`_.

This command can be invoked for one of the Peridynamic :doc:`pair styles <pair_peri>`: peri/eps.

The plasticity value will be 0.0 for atoms not in the specified
compute group.

**Output info:**

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` doc page for an overview of
LAMMPS output options.

The per-atom vector values are unitless numbers (lambda) >= 0.0.

Restrictions
""""""""""""

This compute is part of the PERI package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`compute damage/atom <compute_damage_atom>`,
:doc:`compute dilatation/atom <compute_dilatation_atom>`

**Default:** none

----------

.. _Mitchell:

**(Mitchell)** Mitchell, "A non-local, ordinary-state-based
viscoelasticity model for peridynamics", Sandia National Lab Report,
8064:1-28 (2011).
