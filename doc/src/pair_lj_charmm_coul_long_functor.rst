.. index:: pair_style lj/charmm/coul/long/functor
.. index:: pair_style lj/charmm/coul/long/functor/omp

pair_style lj/charmm/coul/long/functor command
==============================================

Accelerator Variants: *lj/charmm/coul/long/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/charmm/coul/long/functor inner outer (inner2) (outer2)

* inner,outer = inner and outer cutoff for the LJ switching function (distance units)
* inner2,outer2 = optional inner/outer Coulomb cutoff (unused here; the Coulomb
  cutoff is *outer2* if a fourth value is given, else *outer*)

(as for :doc:`pair_style lj/charmm/coul/long <pair_charmm>`, the common forms are
``inner outer`` and ``inner outer cut_coul``)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/charmm/coul/long/functor 8.0 10.0
   pair_style lj/charmm/coul/long/functor 8.0 10.0 9.0
   pair_coeff * * 0.07 3.5
   pair_coeff 1 1 0.07 3.5 0.04 3.0

Description
"""""""""""

.. versionadded:: TBD

Style *lj/charmm/coul/long/functor* computes the CHARMM Lennard-Jones potential
(with the force/energy switching function between the inner and outer cutoff)
together with a long-range (Ewald/PPPM) Coulomb term, and is numerically
equivalent to :doc:`pair_style lj/charmm/coul/long <pair_charmm>`.  It is
provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using
the template/functor framework described in :doc:`Developer_write_pair_functor`.

The arguments and ``pair_coeff`` coefficients (:math:`\epsilon`, :math:`\sigma`,
and optionally the 1-4 :math:`\epsilon_{14}`, :math:`\sigma_{14}`) are as in
:doc:`pair_style lj/charmm/coul/long <pair_charmm>`.  As for that style, the
default mixing rule is arithmetic, the 1-4 LJ coefficients are made available to
:doc:`dihedral_style charmm <dihedral_charmm>`, and a long-range
:doc:`kspace_style <kspace_style>` is required.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.  It additionally requires a :doc:`KSpace style <kspace_style>`
(from the KSPACE package) and a charge per atom.

Related commands
""""""""""""""""

:doc:`pair_style lj/charmm/coul/long <pair_charmm>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
