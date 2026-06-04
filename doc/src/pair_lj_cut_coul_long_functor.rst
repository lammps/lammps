.. index:: pair_style lj/cut/coul/long/functor

pair_style lj/cut/coul/long/functor command
===========================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/functor cutoff (cutoff2)

* cutoff = global cutoff for LJ (and Coulomb if only one cutoff) (distance units)
* cutoff2 = global cutoff for Coulomb (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/long/functor 10.0
   pair_style lj/cut/coul/long/functor 10.0 12.0
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.0 9.0

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/coul/long/functor* computes the cutoff Lennard-Jones potential
together with a long-range (Ewald/PPPM) Coulomb term, and is numerically
equivalent to :doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`.  It is
provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation using
the template/functor framework described in :doc:`Developer_write_pair_functor`
(the Lennard-Jones evaluator combined with a long-range-Coulomb policy).

Like :doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`, it must be used with
a long-range solver such as :doc:`kspace_style ewald or pppm <kspace_style>`, and
it honors the :doc:`pair_modify table <pair_modify>` option (direct ``erfc``
evaluation or the bitmapped Coulomb interpolation tables).

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.  It additionally requires a :doc:`KSpace style <kspace_style>`
(from the KSPACE package) and that atoms store a charge.

The Coulomb cutoff is global (as for the original style).

Related commands
""""""""""""""""

:doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`kspace_style <kspace_style>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
