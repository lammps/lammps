.. index:: pair_style pedone
.. index:: pair_style pedone/omp
.. index:: pair_style pedone/coul/long
.. index:: pair_style pedone/coul/long/omp

pair_style pedone command
=========================

Accelerator Variants: *pedone/omp*

pair_style pedone/coul/long command
===================================

Accelerator Variants: *pedone/coul/long*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = pedone* or *pedone/coul/long*
* args = list of arguments for a particular style

.. parsed-literal::

    *pedone* args = cutoff
      cutoff = global cutoff for Pedone interactions (distance units)
    *pedone/coul/long* args = cutoff (cutoff2)
      cutoff = global cutoff for Pedone (and Coulombic if only one arg) (distance units)
      cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style morse 2.5
   pair_style morse/smooth/linear 2.5
   pair_coeff * * 100.0 2.0 1.5
   pair_coeff 1 1 100.0 2.0 1.5 3.0

Description
"""""""""""

Pair style *pedone* computes the non-Coulomb interactions of the Pedone
(or PMMCS) potential :ref:`Pedone <Pedone>` which combines Coulomb
interactions, a Morse potential, and a repulsive :math:`r^{-12}`
Lennard-Jones term (see below).  The plain *pedone* pair style is meant
to be used in addition to a :doc:`Coulomb pair style <pair_coul>` via
pair style :doc:`hybrid/overlay <pair_hybrid>` and thus allows to be
combined with different Coulomb variants available in LAMMPS.

Pair style *pedone/coul/long* includes the Coulomb part with a damping
function applied so it can be used in conjunction with the
:doc:`kspace_style <kspace_style>` command and its *ewald* or *pppm*
option.  The Coulombic cutoff specified for this style means that
pairwise interactions within this distance are computed directly;
interactions outside that distance are computed in reciprocal space.
This combination is the preferred way to compute the Pedone potential
and should be simpler to use and faster than adding :doc:`pair style
coul/long <pair_coul>` to pair style *pedone* via :doc:`pair style
hybrid/overlay <pair_hybrid>`.

.. math::

   E =  \frac{C q_i q_j}{\epsilon  r}
       + D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
       + \frac{B_0}{r^{12}} \qquad r < r_c

:math:`r_c` is the cutoff and :math:`C` is a conversion factor so that
the entire Coulomb term is in energy units.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

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

None of these pair styles support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

All of these pair styles support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table options are only relevant for
pair style *pedone*

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

-------------

.. _Pedone:

**(Pedone)** A. Pedone, G. Malavasi, M. C. Menziani, A. N. Cormack, and U. Segre, J. Phys. Chem. B, 110, 11780 (2006)
