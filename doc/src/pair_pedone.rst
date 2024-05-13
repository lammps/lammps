.. index:: pair_style pedone
.. index:: pair_style pedone/omp

pair_style pedone command
=========================

Accelerator Variants: *pedone/omp*


Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = pedone*
* args = list of arguments for a particular style

.. parsed-literal::

    *pedone* args = cutoff
      cutoff = global cutoff for Pedone interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

pair_style hybrid/overlay pedone 15.0 coul/long 15.0
kspace_style pppm 1.0e-5

pair_coeff * * coul/long
pair_coeff 1 2 pedone 0.030211 2.241334 2.923245 5.0
pair_coeff 2 2 pedone 0.042395 1.379316 3.618701 22.0

Used in input scripts:

   .. parsed-literal::

      examples/PACKAGES/pedone/in.pedone.relax
      examples/PACKAGES/pedone/in.pedone.melt



Description
"""""""""""

.. versionadded:: 17Apr2024

Pair style *pedone* computes the **non-Coulomb** interactions of the Pedone
(or PMMCS) potential :ref:`(Pedone) <Pedone>` which combines Coulomb
interactions, Morse potential, and repulsive :math:`r^{-12}`
Lennard-Jones terms (see below).  The *pedone* pair style is meant
to be used in addition to a :doc:`Coulomb pair style <pair_coul>` via
pair style :doc:`hybrid/overlay <pair_hybrid>` (see example above).
Using *coul/long* or *could/dsf* (for solids) is recommended.

The full Pedone potential function from :ref:`(Pedone) <Pedone>` for each
pair of atoms is:

.. math::

   E =  \frac{C q_i q_j}{\epsilon  r}
       + D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
       + \frac{B_0}{r^{12}} \qquad r < r_c

:math:`r_c` is the cutoff and :math:`C` is a conversion factor that is
specific to the choice of :doc:`units <units>` so that the entire
Coulomb term is in energy units with :math:`q_i` and :math:`q_j` as the
assigned charges in multiples of the elementary charge.

The following coefficients must be defined for the selected pairs of
atom types via the :doc:`pair_coeff <pair_coeff>` command as in the
example above:

* :math:`D_0` (energy units)
* :math:`\alpha` (1/distance units)
* :math:`r_0` (distance units)
* :math:`C_0` (energy units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global *pedone*
cutoff is used.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.

This pair style support the :doc:`pair_modify <pair_modify>` shift
option for the energy of the pair interaction.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`,
so pair_style and pair_coeff commands does not need to be specified in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, or *outer* keywords.

----------

Restrictions
""""""""""""

The *pedone* pair style is only enabled if LAMMPS was built with the
EXTRA-PAIR package.  See the :doc:`Build package <Build_package>` page
for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style <pair_style>`,
:doc:`pair style coul/long and coul/dsf <pair_coul>`,
:doc:`pair style morse <pair_morse>`

Default
"""""""

none

-------------

.. _Pedone:

**(Pedone)** A. Pedone, G. Malavasi, M. C. Menziani, A. N. Cormack, and U. Segre, J. Phys. Chem. B, 110, 11780 (2006)
