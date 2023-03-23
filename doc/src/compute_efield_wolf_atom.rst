.. index:: compute efield/wolf/atom

compute efield/wolf/atom command
================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID efield/wolf/atom alpha keyword val

* ID, group-ID are documented in :doc:`compute <compute>` command
* efield/atom/wolf = style name of this compute command
* alpha = damping parameter (inverse distance units)
* zero or more keyword/value pairs may be appended
* keyword = *limit* or *cutoff*

  .. parsed-literal::

     *limit* group2-ID = limit computing the electric field contributions to a group (default: all)
     *cutoff* value = set cutoff for computing contributions to this value (default: maximum cutoff of pair style)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all efield/wolf/atom 0.2
   compute 1 mols efield/wolf/atom 0.25 limit water cutoff 10.0

Description
"""""""""""

.. versionadded:: 8Feb2023

Define a computation that approximates the electric field at each atom in a group.

.. math::

   \vec{E}_i = \frac{\vec{F}coul_i}{q_i} = \sum_{j \neq i} \frac{q_j}{r_{ij}^2} \qquad r < r_c

The electric field at the position of the atom *i* is the coulomb force
on a unit charge at that point, which is equivalent to dividing the
Coulomb force by the charge of the individual atom.

In this compute the electric field is approximated as the derivative of
the potential energy using the Wolf summation method, described in
:ref:`Wolf <Wolf4>`, given by:

.. math::
   E_i = \frac{1}{2} \sum_{j \neq i}
   \frac{q_i q_j {\rm erfc}(\alpha r_{ij})}{r_{ij}} +
   \frac{1}{2} \sum_{j \neq i}
   \frac{q_i q_j {\rm erf}(\alpha r_{ij})}{r_{ij}} \qquad r < r_c

where :math:`\alpha` is the damping parameter, and *erf()* and *erfc()*
are error-function and complementary error-function terms.  This
potential is essentially a short-range, spherically-truncated,
charge-neutralized, force-shifted, pairwise *1/r* summation.  With a
manipulation of adding and subtracting a self term (for i = j) to the
first and second term on the right-hand-side, respectively, and a small
enough :math:`\alpha` damping parameter, the second term shrinks and the
potential becomes a rapidly-converging real-space summation.  With a
long enough cutoff and small enough :math:`\alpha` parameter, the
electric field calculated by the Wolf summation method approaches that
computed using the Ewald sum.

The value of the electric field components will be 0.0 for atoms not in
the specified compute group.

When the *limit* keyword is used, only contributions from atoms in the
selected group will be considered, otherwise contributions from all
atoms within the cutoff are included.

When the *cutoff* keyword is used, the cutoff used for the electric
field approximation can be set explicitly.  By default it is the largest
cutoff of any pair style force computation.

.. admonition:: Computational Efficiency
   :class: note

   This compute will loop over a full neighbor list just like a pair
   style does when computing forces, thus it can be quite time-consuming
   and slow down a calculation significantly when its data is used in
   every time step.  The :doc:`compute efield/atom
   <compute_efield_atom>` command of the DIELECTRIC package is more
   efficient in comparison, since the electric field data is collected
   and stored as part of the force computation at next to no extra
   computational cost.

Output info
"""""""""""

This compute calculates a per-atom vector, which can be accessed by
any command that uses per-atom values from a compute as input.  See
the :doc:`Howto output <Howto_output>` page for an overview of
LAMMPS output options.

The vector contains 3 values per atom which are the x-, y-, and
z-direction electric field components in force units.

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package. It is only enabled if
LAMMPS was built with that package.

Related commands
""""""""""""""""

:doc:`pair_style coul/wolf <pair_coul>`,
:doc:`compute efield/atom <compute_efield_atom>`

Default
"""""""

The option defaults are *limit* = all and *cutoff* = largest cutoff
for pair styles.

----------

.. _Wolf4:

**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J Chem
Phys, 110, 8254 (1999).
