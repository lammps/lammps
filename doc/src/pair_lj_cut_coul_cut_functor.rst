.. index:: pair_style lj/cut/coul/cut/functor

pair_style lj/cut/coul/cut/functor command
==========================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/cut/functor cutoff (cutoff2)

* cutoff = global cutoff for LJ (and Coulomb if only one cutoff) (distance units)
* cutoff2 = global cutoff for Coulomb (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/cut/functor 10.0
   pair_style lj/cut/coul/cut/functor 10.0 8.0
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.0 9.0

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/coul/cut/functor* computes the cutoff Lennard-Jones potential
together with a cutoff (unscreened, :math:`1/r`) Coulomb term, and is
numerically equivalent to :doc:`pair_style lj/cut/coul/cut <pair_lj_cut_coul>`.
It is provided by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation
using the template/functor framework described in
:doc:`Developer_write_pair_functor` (the Lennard-Jones evaluator combined with a
cutoff-Coulomb policy).

The ``pair_style`` arguments and the LJ ``pair_coeff`` coefficients are as in
:doc:`pair_style lj/cut/coul/cut <pair_lj_cut_coul>`.

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.

The Coulomb cutoff is global.  Unlike :doc:`pair_style lj/cut/coul/cut
<pair_lj_cut_coul>`, a per-pair Coulomb cutoff (a sixth ``pair_coeff`` argument)
is not yet supported; a per-pair Coulomb cutoff is a planned extension.

Related commands
""""""""""""""""

:doc:`pair_style lj/cut/coul/cut <pair_lj_cut_coul>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`Developer_write_pair_functor`

Default
"""""""

none
