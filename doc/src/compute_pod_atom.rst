.. index:: compute pod/atom
.. index:: compute podd/atom
.. index:: compute pod/local
.. index:: compute pod/global

compute pod/atom command
========================

compute podd/atom command
=========================

compute pod/local command
=========================

compute pod/global command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID pod/atom param.pod coefficients.pod
   compute ID group-ID podd/atom param.pod coefficients.pod
   compute ID group-ID pod/local param.pod coefficients.pod
   compute ID group-ID pod/global param.pod coefficients.pod

* ID, group-ID are documented in :doc:`compute <compute>` command
* pod/atom = style name of this compute command
* param.pod = the parameter file specifies parameters of the POD descriptors
* coefficients.pod = the coefficient file specifies coefficients of the POD potential

Examples
""""""""

.. code-block:: LAMMPS

   compute d all pod/atom Ta_param.pod
   compute dd all podd/atom Ta_param.pod
   compute ldd all pod/local Ta_param.pod
   compute gdd all podd/global Ta_param.pod
   compute d all pod/atom Ta_param.pod Ta_coefficients.pod
   compute dd all podd/atom Ta_param.pod Ta_coefficients.pod
   compute ldd all pod/local Ta_param.pod Ta_coefficients.pod
   compute gdd all podd/global Ta_param.pod Ta_coefficients.pod

Description
"""""""""""

.. versionadded:: 27June2024

Define a computation that calculates a set of quantities related to the
POD descriptors of the atoms in a group. These computes are used
primarily for calculating the dependence of energy and force components
on the linear coefficients in the :doc:`pod pair_style <pair_pod>`,
which is useful when training a POD potential to match target data. POD
descriptors of an atom are characterized by the radial and angular
distribution of neighbor atoms. The detailed mathematical definition is
given in the papers by :ref:`(Nguyen and Rohskopf) <Nguyen20222c>`,
:ref:`(Nguyen2023) <Nguyen20232c>`, :ref:`(Nguyen2024) <Nguyen20242c>`,
and :ref:`(Nguyen and Sema) <Nguyen20243c>`.

Compute *pod/atom* calculates the per-atom POD descriptors.

Compute *podd/atom* calculates derivatives of the per-atom POD
descriptors with respect to atom positions.

Compute *pod/local* calculates the per-atom POD descriptors and their
derivatives with respect to atom positions.

Compute *pod/global* calculates the global POD descriptors and their
derivatives with respect to atom positions.

Examples how to use Compute POD commands are found in the directory
``examples/PACKAGES/pod``.


.. warning::

   All of these compute styles produce *very* large per-atom output
   arrays that scale with the total number of atoms in the system.
   This will result in *very* large memory consumption for systems
   with a large number of atoms.

----------

Output info
"""""""""""

Compute *pod/atom* produces an 2D array of size :math:`N \times M`,
where :math:`N` is the number of atoms and :math:`M` is the number of
descriptors. Each column corresponds to a particular POD descriptor.

Compute *podd/atom* produces an 2D array of size :math:`N \times (M * 3
N)`. Each column corresponds to a particular derivative of a POD
descriptor.

Compute *pod/local* produces an 2D array of size :math:`(1 + 3N) \times
(M * N)`.  The first row contains the per-atom descriptors, and the last
3N rows contain the derivatives of the per-atom descriptors with respect
to atom positions.

Compute *pod/global* produces an 2D array of size :math:`(1 + 3N) \times
(M)`.  The first row contains the global descriptors, and the last 3N
rows contain the derivatives of the global descriptors with respect to
atom positions.

Restrictions
""""""""""""

These computes are part of the ML-POD package.  They are only enabled
if LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fitpod <fitpod_command>`,
:doc:`pair_style pod <pair_pod>`


Default
"""""""

none

----------

.. _Nguyen20222c:

**(Nguyen and Rohskopf)** Nguyen and Rohskopf,  Journal of Computational Physics, 480, 112030, (2023).

.. _Nguyen20232c:

**(Nguyen2023)** Nguyen, Physical Review B, 107(14), 144103, (2023).

.. _Nguyen20242c:

**(Nguyen2024)** Nguyen, Journal of Computational Physics, 113102, (2024).

.. _Nguyen20243c:

**(Nguyen and Sema)** Nguyen and Sema, https://arxiv.org/abs/2405.00306, (2024).


