.. index:: pair_style coul/long/functor

pair_style coul/long/functor command
====================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style coul/long/functor cutoff

* cutoff = global cutoff for Coulomb interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style coul/long/functor 10.0
   pair_coeff * *

Description
"""""""""""

.. versionadded:: TBD

Style *coul/long/functor* computes the long-range (Ewald/PPPM) Coulomb
interaction and is numerically equivalent to :doc:`pair_style coul/long
<pair_coul>`.  It is provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a
reimplementation using the template/functor framework described in
:doc:`Developer_write_pair_functor` (a zero van der Waals evaluator combined
with a long-range-Coulomb policy).

It must be used with a long-range solver such as :doc:`kspace_style ewald or
pppm <kspace_style>`, and honors the :doc:`pair_modify table <pair_modify>`
option.  As for :doc:`pair_style coul/long <pair_coul>`, ``pair_coeff`` takes no
numeric coefficients.

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.  It additionally requires a :doc:`KSpace style <kspace_style>`
(from the KSPACE package) and that atoms store a charge.

Related commands
""""""""""""""""

:doc:`pair_style coul/long <pair_coul>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`kspace_style <kspace_style>`, :doc:`Developer_write_pair_functor`

Default
"""""""

none
