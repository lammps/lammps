.. index:: pair_style beck
.. index:: pair_style beck/gpu
.. index:: pair_style beck/omp

pair_style beck command
=======================

Accelerator Variants: *beck/gpu*, *beck/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style beck Rc

* Rc = cutoff for interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style beck 8.0
   pair_coeff * * 399.671876712 0.0000867636112694 0.675 4.390 0.0003746
   pair_coeff 1 1 399.671876712 0.0000867636112694 0.675 4.390 0.0003746 6.0

Description
"""""""""""

Style *beck* computes interactions based on the potential by
:ref:`(Beck) <Beck>`, originally designed for simulation of Helium.  It
includes truncation at a cutoff distance Rc.

.. math::

   E(r) &= A \exp\left[-\alpha r - \beta r^6\right] - \frac{B}{\left(r^2+a^2\right)^3} \left(1+\frac{2.709+3a^2}{r^2+a^2}\right) \qquad r < R_c \\

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.

* :math:`A` (energy units)
* :math:`B` (energy-distance\^6 units)
* :math:`a` (distance units)
* :math:`\alpha` (1/distance units)
* :math:`\beta`  (1/distance\^6 units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global cutoff
:math:`R_c` is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, coefficients must be specified.
No default mixing rules are used.

This pair style does not support the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the EXTRA-PAIR package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Beck:

**(Beck)** Beck, Molecular Physics, 14, 311 (1968).
