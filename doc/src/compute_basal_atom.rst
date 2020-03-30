.. index:: compute basal/atom

compute basal/atom command
==========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID basal/atom

* ID, group-ID are documented in :doc:`compute <compute>` command
* basal/atom = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all basal/atom

Description
"""""""""""

Defines a computation that calculates the hexagonal close-packed "c"
lattice vector for each atom in the group.  It does this by
calculating the normal unit vector to the basal plane for each atom.
The results enable efficient identification and characterization of
twins and grains in hexagonal close-packed structures.

The output of the compute is thus the 3 components of a unit vector
associated with each atom.  The components are set to 0.0 for
atoms not in the group.

Details of the calculation are given in :ref:`(Barrett) <Barrett>`.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently or to have multiple compute/dump commands, each of
which computes this quantity.

An example input script that uses this compute is provided
in examples/USER/misc/basal.

**Output info:**

This compute calculates a per-atom array with 3 columns, which can be
accessed by indices 1-3 by any command that uses per-atom values from
a compute as input.  See the :doc:`Howto output <Howto_output>` doc page
for an overview of LAMMPS output options.

The per-atom vector values are unitless since the 3 columns represent
components of a unit vector.

Restrictions
""""""""""""

This compute is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

The output of this compute will be meaningless unless the atoms are on
(or near) hcp lattice sites, since the calculation assumes a
well-defined basal plane.

Related commands
""""""""""""""""

:doc:`compute centro/atom <compute_centro_atom>`, :doc:`compute ackland/atom <compute_ackland_atom>`

**Default:** none

----------

.. _Barrett:

**(Barrett)** Barrett, Tschopp, El Kadiri, Scripta Mat. 66, p.666 (2012).
