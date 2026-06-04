.. index:: pair_style lj/cut/functor
.. index:: pair_style lj/cut/functor/omp

pair_style lj/cut/functor command
=================================

Accelerator Variants: *lj/cut/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style lj/cut/functor cutoff

* cutoff = global cutoff for Lennard-Jones interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/functor 2.5
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.0 3.0

Description
"""""""""""

.. versionadded:: TBD

Style *lj/cut/functor* computes the standard cutoff Lennard-Jones potential and
is numerically equivalent to :doc:`pair_style lj/cut <pair_lj>`.  It is provided
by the :ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation of *lj/cut*
using the template/functor framework described in
:doc:`Developer_write_pair_functor`, which applies several inner-loop
optimizations.  The accepted arguments, ``pair_coeff`` coefficients, mixing
rules, and :doc:`pair_modify <pair_modify>` options (``shift``, ``tail``,
``mix``) are identical to :doc:`pair_style lj/cut <pair_lj>`.

The *lj/cut/functor* style can be used by itself or, for type-pair combinations
that have no dedicated fused style, together with other ``/functor`` styles
through :doc:`pair_style hybrid/overlay <pair_hybrid>`.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style supports the :doc:`pair_modify <pair_modify>` *mix*, *shift*, and
*tail* options in the same way as :doc:`pair_style lj/cut <pair_lj>`.  It writes
its information to :doc:`binary restart files <restart>`, so a pair_style and
pair_coeff command do not need to be specified in an input script that reads a
restart file.

This style does not support the :doc:`run_style respa <run_style>` *inner*,
*middle*, *outer* keywords.

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.

Related commands
""""""""""""""""

:doc:`pair_style lj/cut <pair_lj>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
