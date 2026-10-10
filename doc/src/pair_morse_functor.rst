.. index:: pair_style morse/functor
.. index:: pair_style morse/functor/omp

pair_style morse/functor command
================================

Accelerator Variants: *morse/functor/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style morse/functor cutoff

* cutoff = global cutoff for Morse interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style morse/functor 2.5
   pair_coeff * * 100.0 2.0 1.5
   pair_coeff 1 1 100.0 2.0 1.5 3.0

Description
"""""""""""

.. versionadded:: TBD

Style *morse/functor* computes the Morse potential and is numerically equivalent
to :doc:`pair_style morse <pair_morse>`.  It is provided by the
:ref:`FUNCTOR <PKG-FUNCTOR>` package as a reimplementation of *morse* using the
template/functor framework described in :doc:`Developer_write_pair_functor`.  The
accepted arguments and ``pair_coeff`` coefficients (:math:`D_0`, :math:`\alpha`,
:math:`r_0`, and an optional per-pair cutoff) are identical to
:doc:`pair_style morse <pair_morse>`.  Like *morse*, this style has no mixing
rule, so coefficients for all pairs of atom types must be specified explicitly.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style does not support mixing.  It supports the
:doc:`pair_modify <pair_modify>` *shift* option, writes its information to
:doc:`binary restart files <restart>`, and does not support the
:doc:`pair_modify <pair_modify>` *tail* option or the
:doc:`run_style respa <run_style>` *inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

This pair style is part of the FUNCTOR package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` page
for more info.

Related commands
""""""""""""""""

:doc:`pair_style morse <pair_morse>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`Developer_write_pair_functor`

Default
"""""""

none
