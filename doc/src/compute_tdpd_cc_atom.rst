.. index:: compute tdpd/cc/atom

compute tdpd/cc/atom command
============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID tdpd/cc/atom index

* ID, group-ID are documented in :doc:`compute <compute>` command
* tdpd/cc/atom = style name of this compute command
* index = index of chemical species (1 to Nspecies)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all tdpd/cc/atom 2

Description
"""""""""""

Define a computation that calculates the per-atom chemical
concentration of a specified species for each tDPD particle in a
group.

The chemical concentration of each species is defined as the number of
molecules carried by a tDPD particle for dilute solution.  For more
details see :ref:`(Li2015) <Li2015a>`.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input. See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS
output options.

The per-atom vector values will be in the units of chemical species
per unit mass.

Restrictions
""""""""""""

This compute is part of the DPD-MESO package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style tdpd <pair_mesodpd>`

Default
"""""""

none

----------

.. _Li2015a:

**(Li2015)** Li, Yazdani, Tartakovsky, Karniadakis, J Chem Phys, 143:
014101 (2015).  DOI: 10.1063/1.4923254
