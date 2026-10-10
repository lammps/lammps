.. index:: pair_style lj/cut/coul/long/soft/functor
.. index:: pair_style lj/cut/coul/long/soft/functor/omp

pair_style lj/cut/coul/long/soft/functor command
================================================

Accelerator Variants: *lj/cut/coul/long/soft/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/soft/functor n alphalj alphac cutoff (cutoff2)

* n = soft-core coupling exponent for the lambda activation
* alphalj = soft-core parameter for the Lennard-Jones term
* alphac = soft-core parameter for the Coulomb term
* cutoff = global cutoff for LJ (and Coulomb if only one cutoff) (distance units)
* cutoff2 = global cutoff for Coulomb (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/soft/functor 2.0 0.5 10.0 9.5
   pair_coeff * *  1.0 3.0 1.0
   pair_coeff 1 1  1.0 3.0 0.4 9.5

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/coul/long/soft/functor* computes the soft-core (free-energy
perturbation) Lennard-Jones potential together with a soft-core long-range
(Ewald/PPPM) Coulomb term, and is numerically equivalent to :doc:`pair_style
lj/cut/coul/long/soft <pair_fep_soft>`.  It is provided by the
:ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using the
template/functor framework described in :doc:`Developer_write_pair_functor`.

The arguments and ``pair_coeff`` coefficients (:math:`\epsilon`, :math:`\sigma`,
and the per-type-pair coupling parameter :math:`\lambda`) are as in
:doc:`pair_style lj/cut/coul/long/soft <pair_fep_soft>`.  The same :math:`\lambda`
scales both the van der Waals and the Coulomb interaction, which is why this is a
fused style rather than a Coulomb policy applied to an unmodified evaluator.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.  It additionally requires a :doc:`KSpace style <kspace_style>`
(from the KSPACE package) and that atoms store a charge.

Related commands
""""""""""""""""

:doc:`pair_style lj/cut/coul/long/soft <pair_fep_soft>`,
:doc:`pair_style lj/cut/coul/long/functor <pair_lj_cut_coul_long_functor>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
