.. index:: compute force/tally
.. index:: compute heat/flux/tally
.. index:: compute heat/flux/virial/tally
.. index:: compute pe/tally
.. index:: compute pe/mol/tally
.. index:: compute stress/tally

compute force/tally command
===========================

compute heat/flux/tally command
===============================

compute heat/flux/virial/tally command
======================================

compute pe/tally command
========================

compute pe/mol/tally command
============================

compute stress/tally command
============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID style group2-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = *force/tally* or *heat/flux/tally* or *heat/flux/virial/tally* or *pe/tally* or *pe/mol/tally* or *stress/tally*
* group2-ID = group ID of second (or same) group

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 lower force/tally upper
   compute 1 left pe/tally right
   compute 1 lower stress/tally lower
   compute 1 subregion heat/flux/tally all
   compute 1 liquid heat/flux/virial/tally solid

Description
"""""""""""

Define a computation that calculates properties between two groups of
atoms by accumulating them from pairwise non-bonded computations.
Except for *heat/flux/virial/tally*, the two groups can be the same.
This is similar to :doc:`compute group/group <compute_group_group>`
only that the data is
accumulated directly during the non-bonded force computation. The
computes *force/tally*, *pe/tally*, *stress/tally*, and
*heat/flux/tally* are primarily provided as example how to program
additional, more sophisticated computes using the tally callback
mechanism. Compute *pe/mol/tally* is one such style, that can---through using
this mechanism---separately tally intermolecular
and intramolecular energies. Something that would otherwise be
impossible without integrating this as a core functionality into
the base classes of LAMMPS.

----------

Compute *heat/flux/tally* obtains the heat flux
(strictly speaking, heat flow) inside the first group,
which is the sum of the convective contribution
due to atoms in the first group and the virial contribution
due to interaction between the first and second groups:

.. math::

   \mathbf{Q}=  \sum_{i \in \text{group 1}} e_i \mathbf{v}_i + \frac{1}{2} \sum_{i \in \text{group 1}} \sum_{\substack{j \in \text{group 2} \\ j \neq i } } \left( \mathbf{F}_{ij} \cdot \mathbf{v}_j \right) \mathbf{r}_{ij}

When the second group in *heat/flux/tally* is set to "all",
the resulting values will be identical
to that obtained by :doc:`compute heat/flux <compute_heat_flux>`,
provided only pairwise interactions exist.

Compute *heat/flux/virial/tally* obtains the total virial heat flux
(strictly speaking, heat flow) into the first group due to interaction
with the second group, and is defined as:

.. math::

   Q = \frac{1}{2} \sum_{i \in \text{group 1}} \sum_{j \in \text{group 2}} \mathbf{F}_{ij} \cdot \left(\mathbf{v}_i + \mathbf{v}_j \right)

Although, the *heat/flux/virial/tally* compute
does not include the convective term,
it can be used to obtain the total heat flux over control surfaces,
when there are no particles crossing over,
such as is often in solid--solid and solid--liquid interfaces.
This would be identical to the method of planes method.
Note that the *heat/flux/virial/tally* compute is distinctly different
from the *heat/flux* and *heat/flux/tally* computes,
that are essentially volume averaging methods.
The following example demonstrates the difference:

.. code-block:: LAMMPS

   # System with only pairwise interactions.
   # Non-periodic boundaries in the x direction.
   # Has LeftLiquid and RightWall groups along x direction.

   # Heat flux over the solid-liquid interface
   compute hflow_hfvt RightWall heat/flux/virial/tally LeftLiquid
   variable hflux_hfvt equal c_hflow_hfvt/(ly*lz)

   # x component of approximate heat flux vector inside the liquid region,
   # two approaches.
   #
   compute myKE all ke/atom
   compute myPE all pe/atom
   compute myStress all stress/atom NULL virial
   compute hflow_hf LeftLiquid heat/flux myKE myPE myStress
   variable hflux_hf equal c_hflow_hf[1]/${volLiq}
   #
   compute hflow_hft LeftLiquid heat/flux/tally all
   variable hflux_hft equal c_hflow_hft[1]/${volLiq}

   # Pressure over the solid-liquid interface, three approaches.
   #
   compute force_gg RightWall group/group LeftLiquid
   variable press_gg equal c_force_gg[1]/(ly*lz)
   #
   compute force_ft RightWall force/tally LeftLiquid
   compute rforce_ft RightWall reduce sum c_force_ft[1]
   variable press_ft equal c_rforce_ft/(ly*lz)
   #
   compute rforce_hfvt all reduce sum c_hflow_hfvt[1]
   variable press_hfvt equal c_rforce_hfvt/(ly*lz)

----------

The pairwise contributions are computing via a callback that the
compute registers with the non-bonded pairwise force computation.
This limits the use to systems that have no bonds, no Kspace, and no
many-body interactions. On the other hand, the computation does not
have to compute forces or energies a second time and thus can be much
more efficient. The callback mechanism allows to write more complex
pairwise property computations.

----------

Output info
"""""""""""

- Compute *pe/tally* calculates a global scalar (the energy) and a per
  atom scalar (the contributions of the single atom to the global
  scalar).

- Compute *pe/mol/tally* calculates a global four-element vector containing
  (in this order): *evdwl* and *ecoul* for intramolecular pairs and
  *evdwl* and *ecoul* for intermolecular pairs. Since molecules are
  identified by their molecule IDs, the partitioning does not have to be
  related to molecules, but the energies are tallied into the respective
  slots depending on whether the molecule IDs of a pair are the same or
  different.

- Compute *force/tally* calculates a global scalar (the force magnitude)
  and a per atom 3-element vector (force contribution from each atom).

- Compute *stress/tally* calculates a global scalar
  (average of the diagonal elements of the stress tensor) and a per atom
  vector (the six elements of stress tensor contributions from the
  individual atom).

- As in :doc:`compute heat/flux <compute_heat_flux>`,
  compute *heat/flux/tally* calculates a global vector of length 6,
  where the first three components are the :math:`x`, :math:`y`, :math:`z`
  components of the full heat flow vector,
  and the next three components are the corresponding components
  of just the convective portion of the flow (i.e., the
  first term in the equation for :math:`\mathbf{Q}`).

- Compute *heat/flux/virial/tally* calculates a global scalar (heat flow)
  and a per atom three-element vector
  (contribution to the force acting over atoms in the first group
  from individual atoms in both groups).

Both the scalar and vector values calculated by this compute are
"extensive".

Restrictions
""""""""""""

This compute is part of the TALLY package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Not all pair styles can be evaluated in a pairwise mode as required by
this compute.  For example, 3-body and other many-body potentials,
such as :doc:`Tersoff <pair_tersoff>` and
:doc:`Stillinger-Weber <pair_sw>` cannot be used.  :doc:`EAM <pair_eam>`
potentials only include the pair potential portion of the EAM
interaction when used by this compute, not the embedding term.  Also
bonded or Kspace interactions do not contribute to this compute.

When used with dynamic groups, a :doc:`run 0 <run>` command needs to
be inserted in order to initialize the dynamic groups before accessing
the computes.

Related commands
""""""""""""""""

* :doc:`compute group/group <compute_group_group>`
* :doc:`compute heat/flux <compute_heat_flux>`

Default
"""""""

none
