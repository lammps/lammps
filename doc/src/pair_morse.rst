.. index:: pair_style morse
.. index:: pair_style morse/gpu
.. index:: pair_style morse/omp
.. index:: pair_style morse/opt
.. index:: pair_style morse/smooth/linear
.. index:: pair_style morse/smooth/linear/omp
.. index:: pair_style morse/kk

pair_style morse command
========================

Accelerator Variants: *morse/gpu*, *morse/omp*, *morse/opt*, *morse/kk*

pair_style morse/smooth/linear command
======================================

Accelerator Variants: *morse/smooth/linear/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *morse* or *morse/smooth/linear* or *morse/soft*
* args = list of arguments for a particular style

.. parsed-literal::

    *morse* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)
    *morse/smooth/linear* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style morse 2.5
   pair_style morse/smooth/linear 2.5
   pair_coeff * * 100.0 2.0 1.5
   pair_coeff 1 1 100.0 2.0 1.5 3.0

Description
"""""""""""

Style *morse* computes pairwise interactions with the formula

.. math::

   E = D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
       \qquad r < r_c

Rc is the cutoff.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`D_0` (energy units)
* :math:`\alpha` (1/distance units)
* :math:`r_0` (distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global morse
cutoff is used.

----------

The *morse/smooth/linear* variant is similar to the lj/smooth/linear
variant in that it adds to the potential a shift and a linear term
so that both, potential energy and force, go to zero at the cut-off:

.. math::

   \phi\left(r\right) & =  D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right] \qquad r < r_c \\
   E\left(r\right) & =  \phi\left(r\right)  - \phi\left(R_c\right) - \left(r - R_c\right) \left.\frac{d\phi}{d r} \right|_{r=R_c}       \qquad r < R_c

The syntax of the pair_style and pair_coeff commands are the same for
the *morse* and *morse/smooth/linear* styles.

----------

A version of the *morse* style with a soft core, *morse/soft*,
suitable for use in free energy calculations, is part of the FEP
package and is documented with the :doc:`pair_style */soft
<pair_fep_soft>` styles. The version with soft core is only available if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

None of these pair styles support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

All of these pair styles support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table options is not relevant for
the Morse pair styles.

None of these pair styles support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

All of these pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *morse/smooth/linear* pair style is only enabled if LAMMPS was
built with the EXTRA-PAIR package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style */soft <pair_fep_soft>`

Default
"""""""

none
