.. index:: pair_style coul/long/soft/functor
.. index:: pair_style coul/long/soft/functor/omp

pair_style coul/long/soft/functor command
=========================================

Accelerator Variants: *coul/long/soft/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style coul/long/soft/functor n alphac cutoff

* n = soft-core coupling exponent for the lambda activation
* alphac = soft-core parameter for the Coulomb term
* cutoff = global cutoff for Coulomb interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style coul/long/soft/functor 2.0 10.0 9.5
   pair_coeff * * 1.0
   pair_coeff 1 1 0.4

Description
"""""""""""

.. versionadded:: TBD

Style *coul/long/soft/functor* computes the soft-core (free-energy perturbation)
long-range (Ewald/PPPM) Coulomb interaction and is numerically equivalent to
:doc:`pair_style coul/long/soft <pair_fep_soft>`.  It is provided by the
:ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using the
template/functor framework described in :doc:`Developer_write_pair_functor` (a
zero van der Waals evaluator combined with a soft-core long-range-Coulomb policy).

As for :doc:`pair_style coul/long/soft <pair_fep_soft>`, the only ``pair_coeff``
argument is the per-type-pair coupling parameter :math:`\lambda`.  It must be used
with a long-range solver such as :doc:`kspace_style ewald or pppm <kspace_style>`.

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

:doc:`pair_style coul/long/soft <pair_fep_soft>`,
:doc:`pair_style coul/long/functor <pair_coul_long_functor>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
