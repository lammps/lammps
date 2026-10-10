.. index:: pair_style coul/cut/soft/functor
.. index:: pair_style coul/cut/soft/functor/omp

pair_style coul/cut/soft/functor command
========================================

Accelerator Variants: *coul/cut/soft/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style coul/cut/soft/functor n alphac cutoff

* n = soft-core coupling exponent for the lambda activation
* alphac = soft-core parameter for the Coulomb term
* cutoff = global cutoff for Coulomb interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style coul/cut/soft/functor 2.0 10.0 9.5
   pair_coeff * * 1.0
   pair_coeff 1 1 0.4 8.0

Description
"""""""""""

.. versionadded:: TBD

Style *coul/cut/soft/functor* computes the soft-core (free-energy perturbation)
cutoff Coulomb interaction and is numerically equivalent to :doc:`pair_style
coul/cut/soft <pair_fep_soft>`.  It is provided by the :ref:`FUNCTOR
<PKG-FUNCTOR>` package as a reimplementation using the template/functor framework
described in :doc:`Developer_write_pair_functor` (a zero van der Waals evaluator
combined with a soft-core cutoff-Coulomb policy).

As for :doc:`pair_style coul/cut/soft <pair_fep_soft>`, the ``pair_coeff``
arguments are the per-type-pair coupling parameter :math:`\lambda` and an optional
per-type-pair cutoff; the Coulomb cutoff of unset pairs is obtained by mixing.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.

This style requires that atoms store a charge, e.g. via :doc:`atom_style full
<atom_style>` or :doc:`atom_style charge <atom_style>`.

Related commands
""""""""""""""""

:doc:`pair_style coul/cut/soft <pair_fep_soft>`,
:doc:`pair_style coul/cut/functor <pair_coul_cut_functor>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`Developer_write_pair_functor`

Default
"""""""

none
