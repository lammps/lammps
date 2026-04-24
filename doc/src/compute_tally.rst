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

   compute ID group-ID style group2-ID keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* style = *force/tally* or *heat/flux/tally* or *heat/flux/virial/tally* or *pe/tally* or *pe/mol/tally* or *stress/tally*
* group2-ID = group ID of second (or same) group
* zero or more keywords may be appended
* keyword = :code:`two_body` and/or :code:`three_body` and/or :code:`four_body`

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 lower force/tally upper
   compute 1 left pe/tally right
   compute 1 lower stress/tally lower
   compute 1 subregion heat/flux/tally all
   compute 1 liquid heat/flux/virial/tally solid
   compute 1 left heat/flux/virial/tally right two_body four_body

Description
"""""""""""

Define a computation that calculates properties between two groups of
atoms by accumulating them from force and potential energy computations.
This is similar to :doc:`compute group/group <compute_group_group>`
only that the data is
accumulated directly during the force computation by a callback mechanism.
By default, only non-bonded pairwise interactions are supported, however the compute *heat/flux/virial/tally* has been extended to calculate contributions from many-body interactions. Also, for this compute the two groups cannot be the same. The
computes *force/tally*, *pe/tally*, *stress/tally*, and
*heat/flux/tally* are primarily provided as example how to program
additional, more sophisticated computes using the tally callback
mechanism. Compute *pe/mol/tally* is one such style, that can---through using
this mechanism---separately tally intermolecular
and intramolecular energies. Something that would otherwise be
impossible without integrating this as a core functionality into
the base classes of LAMMPS.

----------

For pairwise interactions, compute *heat/flux/tally* obtains the heat flux
(strictly speaking, heat flow) inside the first group,
which is the sum of the convective contribution
due to atoms in the first group and the virial contribution, strictly due to pairwise non-bonded interactions between the first and second groups:

.. math::

   \mathbf{Q}_{2 \rightarrow 1}=  \sum_{i \in \text{group 1}} e_i \mathbf{v}_i + \frac{1}{2} \sum_{i \in \text{group 1}} \sum_{\substack{j \in \text{group 2} \\ j \neq i } } \left( \mathbf{F}_{ij} \cdot \mathbf{v}_j \right) \mathbf{r}_{ij}

When the second group in *heat/flux/tally* is set to "all",
the resulting values will be identical
to that obtained by :doc:`compute heat/flux <compute_heat_flux>`,
provided only pairwise interactions exist.

Compute *heat/flux/virial/tally* obtains the total virial heat flux
(strictly speaking, heat flow) into the first group due to interaction
with the second group. For pairwise interactions, it is given by:

.. math::

   Q^{\text{virial}}_{2 \rightarrow 1} = \frac{1}{2} \sum_{i \in \text{group 1}} \sum_{j \in \text{group 2}} \mathbf{F}_{ij} \cdot \left(\mathbf{v}_i + \mathbf{v}_j \right)

.. versionchanged:: TBD

For many-body-interactions, the above expression needs to be modified. Details of this derivation can be found in :ref:`(Poulos2026)  <Poulos1>`. The resulting expression for the virial heat flow is given by:

.. math::

  Q^{\text{virial}}_{2 \rightarrow 1} = \sum_i \mathbf{F}^{2 \rightarrow 1}_i \cdot \mathbf{v}_i


where now :math:`i` runs over all atoms in both groups. The many-body force term :math:`\mathbf{F}^{2 \rightarrow 1}_i` is defined as

.. math::
    \mathbf{F}^{2 \rightarrow 1}_i &=   \sum_{k}^{K_N} \left(H_i^k - H_0^k \right) \mathbf{F}_i^k \\
    H_0^k &= \sum_{j \in k} p_j^k H_j \\
    H_i &= \begin{cases}
    1 & \text{if } i \in 1 \\
    0 & \text{otherwise}
    \end{cases}

where :math:`k` runs over all :math:`K_N` many-body potential energy terms that include atom :math:`i` and :math:`p_j^k` is the percentage of the potential energy term :math:`U^k` attributed to atom :math:`j`. For each atom :math:`i` and each interaction term :math:`k`, the factor :math:`( H_i^k - H_0^k )` is equal to the total percentage (in absolute value) of the potential energy term :math:`U^k` attributed to all the atoms contained in the interaction term :math:`k` that belong to the opposite group than :math:`i`. For pairwise interactions, the above expressions reduce to the standard formula provided above (:math:`\mathbf{F}^{2 \rightarrow 1}_i = \sum_{j} \mathbf{F}_{ij}`).

The :code:`two_body`, :code:`three_body` and :code:`four_body` keywords are only available for the *heat/flux/virial/tally* compute and function as flags controlling the inclusion of 2-body, 3-body and 4-body interactions terms in the calculation of :math:`\mathbf{F}^{2 \rightarrow 1}_i`.

The *heat/flux/virial/tally* compute can also be used to easily obtain the spectral decomposition of the heat current with many-body interactions, as described in :ref:`(Poulos2026)  <Poulos1_tally>`.

.. math::

  Q^{\text{virial}}_{2 \rightarrow 1}(\omega) = \frac{2}{\tau_{total}}\sum_i Re\{\langle \hat{\mathbf{F}}^{2 \rightarrow 1}_i(\omega) \cdot \hat{\mathbf{v}}^*_i(\omega) \rangle \}

where :math:`\hat{A}(\omega)` is the Fourier transforms of the time-dependent property :math:`A(t)`, :math:`*` denotes complex conjugate and :math:`\tau_{total}` is the total simulation time. Below is an example of how to obtain :math:`Q^{\text{virial}}_{2 \rightarrow 1}(\omega)` through LAMMPS using the *heat/flux/virial/tally* compute:

.. code-block:: LAMMPS

   # Calculate the heat flow from group G_B to G_A,
   # and its spectral decomposition.
   # Group 'Interface' contains the interacting atoms from both groups.

   # i. Variables
   variable thermo_steps equal 1000
   variable N_repeat_TP  equal 100
   variable N_every_TP   equal 10
   variable N_force      equal 10

   # ii. Calculations
   compute  force_mb   G_A heat/flux/virial/tally G_B
   fix      ave_heat   all ave/time ${N_every_TP} ${N_repeat_TP} ${thermo_steps} c_force_mb

   # iii. Outputting
   thermo         ${thermo_steps}
   thermo_style   custom step temp etotal f_ave_heat
   dump     data  Interface custom ${N_force} Virial.dat id c_force_mb[1] c_force_mb[2] c_force_mb[3] vx vy vz

The :code:`f_ave_heat` variable then gives the instantaneous heat current, while the Fourier transforms of the :code:`c_force_mb` vector and the atomic velocities :code:`vx`, :code:`vy`, :code:`vz` can be used to compute the spectral decomposition :math:`Q(\omega)` as detailed above.

The *heat/flux/virial/tally* compute enables the calculation of both the instantaneous heat flow as well as its spectral decomposition across any arbitrary control surface defined as the physical boundary between the two groups. This thus enables the study of heat flow across interfaces of any arbitrary geometry and not necessarily planar ones.

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

The force contributions are computed via a callback that the
compute registers with the force computation (strictly non-bonded pairwise, except for *heat/flux/virial/tally*).
This limits the use to systems that have no bonds, no Kspace, and no
many-body interactions (except for *heat/flux/virial/tally*). On the other hand, the computation does not
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
  (the force :math:`\mathbf{F}^{2 \rightarrow 1}_i`, which is the contribution to the force acting over each atom in the first group
  from all atoms in the second group).

Both the scalar and vector values calculated by this compute are
"extensive".

Restrictions
""""""""""""

This compute is part of the TALLY package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Not all pair styles can be evaluated in a pairwise mode as required by
this compute.  Currently, only the *heat/flux/virial/tally* compute supports many-body interactions, and to date only the :doc:`Tersoff <pair_tersoff>` and :doc:`Stillinger-Weber <pair_sw>` many-body potentials have been implemented. For all other computes, pair styles with many-body interactions
cannot be used. :doc:`EAM <pair_eam>` potentials only include the pair
potential portion of the EAM interaction when used by this compute, not
the embedding term.  Also bonded or Kspace interactions do not
contribute to this compute.

The :code:`two_body`, :code:`three_body`, and :code:`four_body` keywords are only available for the *heat/flux/virial/tally* compute.

These computes are not compatible with accelerated pair styles from the
GPU, INTEL, KOKKOS, or OPENMP packages. They will either create an error
or print a warning when required data was not tallied in the required way
and thus the data acquisition functions from these computes not called.

When used with dynamic groups, a :doc:`run 0 <run>` command needs to
be inserted in order to initialize the dynamic groups before accessing
the computes.

Related commands
""""""""""""""""

* :doc:`compute group/group <compute_group_group>`
* :doc:`compute heat/flux <compute_heat_flux>`

Default
"""""""

By default, the compute includes contributions from all many-body interactions, that is, the keywords :code:`two_body`, :code:`three_body` and :code:`four_body` are all activated by default.

----------

.. _Poulos1_tally:

**(Poulos2026)** Poulos, Surblys, Termentzidis, Phys Rev B, 113, 045414 (2026).
