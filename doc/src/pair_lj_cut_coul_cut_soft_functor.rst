.. index:: pair_style lj/cut/coul/cut/soft/functor
.. index:: pair_style lj/cut/coul/cut/soft/functor/omp

pair_style lj/cut/coul/cut/soft/functor command
===============================================

Accelerator Variants: *lj/cut/coul/cut/soft/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/cut/soft/functor n alphalj alphac cutoff (cutoff2)

* n = soft-core coupling exponent for the lambda activation
* alphalj = soft-core parameter for the Lennard-Jones term
* alphac = soft-core parameter for the Coulomb term
* cutoff = global cutoff for LJ (and Coulomb if only one cutoff) (distance units)
* cutoff2 = global cutoff for Coulomb (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/cut/soft/functor 2.0 0.5 10.0 9.5
   pair_coeff * *  1.0 3.0 1.0
   pair_coeff 1 1  1.0 3.0 0.4 9.5

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/coul/cut/soft/functor* computes the soft-core (free-energy
perturbation) Lennard-Jones potential together with a soft-core cutoff Coulomb
term, and is numerically equivalent to :doc:`pair_style lj/cut/coul/cut/soft
<pair_fep_soft>`.  It is provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a
reimplementation using the template/functor framework described in
:doc:`Developer_write_pair_functor`.

The arguments and ``pair_coeff`` coefficients (:math:`\epsilon`, :math:`\sigma`,
and the per-type-pair coupling parameter :math:`\lambda`) are as in
:doc:`pair_style lj/cut/coul/cut/soft <pair_fep_soft>`.  The same :math:`\lambda`
scales both the van der Waals and the Coulomb interaction, which is why this is a
fused style rather than a Coulomb policy applied to an unmodified evaluator.

A per-pair Coulomb cutoff is supported only through the seven-argument
``pair_coeff`` form (``I J epsilon sigma lambda cut_lj cut_coul``).  Unlike
:doc:`pair_style lj/cut/coul/cut/soft <pair_fep_soft>`, the six-argument form
(``I J epsilon sigma lambda cut_lj``) does *not* also set the Coulomb cutoff to
``cut_lj``; the Coulomb cutoff then remains the global value from the
``pair_style`` command.

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

:doc:`pair_style lj/cut/coul/cut/soft <pair_fep_soft>`,
:doc:`pair_style lj/cut/coul/cut/functor <pair_lj_cut_coul_cut_functor>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`Developer_write_pair_functor`

Default
"""""""

none
