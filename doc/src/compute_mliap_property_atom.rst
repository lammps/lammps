.. index:: compute mliap/property/atom
.. index:: compute mliap/property/atom/kk

compute mliap/property/atom command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID mliap/property/atom name property_name dim property_dim hybrid_index index

* ID, group-ID are documented in :doc:`compute <compute>` command
* mliap/property/atom = style name of this compute command
* name property_name = specifies the name of the extra property to get from the unified MLIAP
* dim property_dim = specifies the dimension of the per-atom vector to grab from the unified MLIAP (must be >= 1)
* hybrid_index index = specifies which mliap pair style (see :doc:`pair_style mliap <pair_mliap>` for more information on mliap pair style) in a hybrid pair style to grab (see :doc:`pair_style hybrid/scaled <pair_hybrid>` for more information about the hybrid pair style)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all mliap/property/atom name energy_std dim 1
   compute 1 all mliap/property/atom name pair_force_std dim 3

   pair_style hybrid/scaled &
                         0.3 mliap unified mliap_unified_model1.pt &
                         0.7 mliap unified mliap_unified_model2.pt &
   compute 1 all mliap/property/atom name energy_std_2 dim 1 hybrid_index 2


Description
"""""""""""

Define a computation that computes a per-atom vector of some named extra property
that the unified MLIAP produces. This compute will register a request for an extra property
with the given name and dimension. This request occurs every time the `pair mliap unified`
pair style (see :doc:`pair_style mliap unified <pair_mliap>` for more information about the MLIAP pair style).

.. note::

    This compute requires a :doc:`pair_style mliap unified <pair_mliap>` to function.
   
Upon request, a call is made to the MLIAPUnifiedInterface that acts as a request to the MLIAP to
update the extra property with the given name (as a Python string). After the MLIAPUnifiedInterface
returns from the request, the data is copied into the ExtraProperties struct for access when needed
from the compute.

As an example of use of this compute, we demonstrate how this compute can be used to dump a per-atom
potential energy standard deviation that is returned by an ensemble model. Suppose that the
MLIAP has an extra property named 'energy_std'.

.. code-block:: LAMMPS

   compute        peAtom all pe/atom
   compute        peStdAtom all mliap/property/atom name energy_std dim 1
   dump 1 all custom 1 eatoms_std.xyz id c_peAtom c_peStdAtom[1]

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

Restrictions
""""""""""""

The hybrid_index feature is currently only allowed in non-Kokkos LAMMPS.

Related commands
""""""""""""""""

:doc:`pair_style mliap <pair_mliap>`
:doc:`pair_style hybrid <pair_hybrid>`

Default
"""""""

none
