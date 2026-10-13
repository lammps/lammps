.. index:: pair_style coul/cut/functor
.. index:: pair_style coul/cut/functor/omp

pair_style coul/cut/functor command
===================================

Accelerator Variants: *coul/cut/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style coul/cut/functor cutoff

* cutoff = global cutoff for Coulomb interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style coul/cut/functor 10.0
   pair_coeff * *
   pair_coeff 1 1 8.0

Description
"""""""""""

.. versionadded:: TBD

Style *coul/cut/functor* computes the cutoff (unscreened, :math:`1/r`) Coulomb
interaction and is numerically equivalent to :doc:`pair_style coul/cut
<pair_coul>`.  It is provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a
reimplementation using the template/functor framework described in
:doc:`Developer_write_pair_functor` (a zero van der Waals evaluator combined with
a cutoff-Coulomb policy).

As for :doc:`pair_style coul/cut <pair_coul>`, the only ``pair_coeff`` argument
is an optional per-type-pair cutoff; the Coulomb cutoff of unset pairs is
obtained by mixing.

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

:doc:`pair_style coul/cut <pair_coul>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
