.. index:: pair_style yukawa
.. index:: pair_style yukawa/gpu
.. index:: pair_style yukawa/omp
.. index:: pair_style yukawa/kk

pair_style yukawa command
=========================

Accelerator Variants: *yukawa/gpu*, *yukawa/omp*, *yukawa/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style yukawa kappa cutoff

* kappa = screening length (inverse distance units)
* cutoff = global cutoff for Yukawa interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style yukawa 2.0 2.5
   pair_coeff 1 1 100.0 2.3
   pair_coeff * * 100.0

Description
"""""""""""

Style *yukawa* computes pairwise interactions with the formula

.. math::

   E = A \frac{e^{- \kappa r}}{r} \qquad r < r_c

:math:`r_c` is the cutoff.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* A (energy\*distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global yukawa
cutoff is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the A coefficient and cutoff
distance for this pair style can be mixed.  A is an energy value mixed
like a LJ epsilon.  The default mix value is *geometric*\ .  See the
"pair_modify" command for details.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
