.. index:: pair_style coul/long/cs/functor
.. index:: pair_style coul/long/cs/functor/omp

pair_style coul/long/cs/functor command
=======================================

Accelerator Variants: *coul/long/cs/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style coul/long/cs/functor cutoff

* cutoff = global cutoff for Coulomb interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style coul/long/cs/functor 10.0
   pair_coeff * *

Description
"""""""""""

.. versionadded:: TBD

Style *coul/long/cs/functor* is the core-shell (adiabatic core-shell / Drude)
variant of :doc:`pair_style coul/long/functor <pair_coul_long_functor>` and is
numerically equivalent to :doc:`pair_style coul/long/cs <pair_cs>`.  It is
provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using
the template/functor framework described in :doc:`Developer_write_pair_functor`.

As for :doc:`pair_style coul/long/cs <pair_cs>`, the long-range Coulomb kernel is
regularized so that two oppositely charged core/shell particles that approach
each other very closely (:math:`r \rightarrow 0`) do not cause a division by
zero: a tiny offset is added to :math:`r^2`, a higher-accuracy series expansion
of the Ewald correction is used, and bonded pairs are treated with a small
displacement of the singular :math:`1/r` term.  This is intended for use with
the core-shell model; see the :doc:`core-shell Howto <Howto_coreshell>`.

It must be used with a long-range solver such as :doc:`kspace_style ewald or
pppm <kspace_style>`, and honors the :doc:`pair_modify table <pair_modify>`
option.  As for :doc:`pair_style coul/long/cs <pair_cs>`, ``pair_coeff`` takes no
numeric coefficients, and the :doc:`single command <Howto_coreshell>` is not
supported (the regularized kernel does not match the energy decomposition used by
:doc:`compute group/group <compute_group_group>`).

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

:doc:`pair_style coul/long/cs <pair_cs>`,
:doc:`pair_style coul/long/functor <pair_coul_long_functor>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
