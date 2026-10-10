.. index:: pair_style lj/cut/coul/long/cs/functor
.. index:: pair_style lj/cut/coul/long/cs/functor/omp

pair_style lj/cut/coul/long/cs/functor command
==============================================

Accelerator Variants: *lj/cut/coul/long/cs/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/cs/functor cutoff (cutoff2)

* cutoff = global cutoff for LJ (and Coulomb if only one cutoff) (distance units)
* cutoff2 = global cutoff for Coulomb (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/cs/functor 10.0
   pair_style lj/cut/coul/long/cs/functor 10.0 8.0
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.0 9.0

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/coul/long/cs/functor* is the core-shell (adiabatic core-shell /
Drude) variant of :doc:`pair_style lj/cut/coul/long/functor
<pair_lj_cut_coul_long_functor>` and is numerically equivalent to
:doc:`pair_style lj/cut/coul/long/cs <pair_cs>`.  It is provided by the
:ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using the
template/functor framework described in :doc:`Developer_write_pair_functor`.

It computes the cutoff Lennard-Jones potential together with the regularized
long-range (Ewald/PPPM) Coulomb kernel of the core-shell model (see
:doc:`pair_style coul/long/cs <pair_cs>` and the :doc:`core-shell Howto
<Howto_coreshell>`).  The arguments and ``pair_coeff`` coefficients are as in
:doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`.

As for :doc:`pair_style lj/cut/coul/long/cs <pair_cs>`, the :doc:`single command
<Howto_coreshell>` is not supported, because the regularized core-shell Coulomb
kernel does not match the energy decomposition used by :doc:`compute group/group
<compute_group_group>`.

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

:doc:`pair_style lj/cut/coul/long/cs <pair_cs>`,
:doc:`pair_style lj/cut/coul/long/functor <pair_lj_cut_coul_long_functor>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
