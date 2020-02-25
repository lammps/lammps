.. index:: pair_style born

pair_style born command
=======================

pair_style born/omp command
===========================

pair_style born/gpu command
===========================

pair_style born/coul/long command
=================================

pair_style born/coul/long/gpu command
=====================================

pair_style born/coul/long/omp command
=====================================

pair_style born/coul/msm command
================================

pair_style born/coul/msm/omp command
====================================

pair_style born/coul/wolf command
=================================

pair_style born/coul/wolf/gpu command
=====================================

pair_style born/coul/wolf/omp command
=====================================

pair_style born/coul/dsf command
================================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style style args

* style = *born* or *born/coul/long* or *born/coul/msm* or *born/coul/wolf*
* args = list of arguments for a particular style


.. parsed-literal::

     *born* args = cutoff
       cutoff = global cutoff for non-Coulombic interactions (distance units)
     *born/coul/long* args = cutoff (cutoff2)
       cutoff = global cutoff for non-Coulombic (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *born/coul/msm* args = cutoff (cutoff2)
       cutoff = global cutoff for non-Coulombic (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *born/coul/wolf* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for non-Coulombic (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *born/coul/dsf* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for non-Coulombic (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (distance units)

Examples
""""""""


.. code-block:: LAMMPS

   pair_style born 10.0
   pair_coeff * * 6.08 0.317 2.340 24.18 11.51
   pair_coeff 1 1 6.08 0.317 2.340 24.18 11.51

   pair_style born/coul/long 10.0
   pair_style born/coul/long 10.0 8.
   pair_coeff * * 6.08 0.317 2.340 24.18 11.51
   pair_coeff 1 1 6.08 0.317 2.340 24.18 11.51

   pair_style born/coul/msm 10.0
   pair_style born/coul/msm 10.0 8.0
   pair_coeff * * 6.08 0.317 2.340 24.18 11.51
   pair_coeff 1 1 6.08 0.317 2.340 24.18 11.51

   pair_style born/coul/wolf 0.25 10.0
   pair_style born/coul/wolf 0.25 10.0 9.0
   pair_coeff * * 6.08 0.317 2.340 24.18 11.51
   pair_coeff 1 1 6.08 0.317 2.340 24.18 11.51

   pair_style born/coul/dsf 0.1 10.0 12.0
   pair_coeff * *   0.0 1.00 0.00 0.00 0.00
   pair_coeff 1 1 480.0 0.25 0.00 1.05 0.50

Description
"""""""""""

The *born* style computes the Born-Mayer-Huggins or Tosi/Fumi
potential described in :ref:`(Fumi and Tosi) <FumiTosi>`, given by

.. math::

   E = A \exp \left(\frac{\sigma - r}{\rho} \right) - 
   \frac{C}{r^6} + \frac{D}{r^8} \qquad r < r_c


where :math:`\sigma` is an interaction-dependent length parameter,
:math:`\rho` is an ionic-pair dependent length parameter, and
:math:`r_c` is the cutoff.

The styles with *coul/long* or *coul/msm* add a Coulombic term as
described for the :doc:`lj/cut <pair_lj>` pair styles.  An additional
damping factor is applied to the Coulombic term so it can be used in
conjunction with the :doc:`kspace_style <kspace_style>` command and its
*ewald* or *pppm* of *msm* option.  The Coulombic cutoff specified for
this style means that pairwise interactions within this distance are
computed directly; interactions outside that distance are computed in
reciprocal space.

If one cutoff is specified for the *born/coul/long* and
*born/coul/msm* style, it is used for both the A,C,D and Coulombic
terms.  If two cutoffs are specified, the first is used as the cutoff
for the A,C,D terms, and the second is the cutoff for the Coulombic
term.

The *born/coul/wolf* style adds a Coulombic term as described for the
Wolf potential in the :doc:`coul/wolf <pair_coul>` pair style.

The *born/coul/dsf* style computes the Coulomb contribution with the
damped shifted force model as in the :doc:`coul/dsf <pair_coul>` style.

Note that these potentials are related to the :doc:`Buckingham potential <pair_buck>`.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* A (energy units)
* :math:`\rho` (distance units)
* :math:`\sigma` (distance units)
* C (energy units \* distance units\^6)
* D (energy units \* distance units\^8)
* cutoff (distance units)

The second coefficient, rho, must be greater than zero.

The last coefficient is optional.  If not specified, the global A,C,D
cutoff specified in the pair\_style command is used.

For *born/coul/long*\ , *born/coul/wolf* and *born/coul/dsf* no
Coulombic cutoff can be specified for an individual I,J type pair.
All type pairs use the same global Coulombic cutoff specified in the
pair\_style command.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

These pair styles do not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

These styles support the :doc:`pair_modify <pair_modify>` shift option
for the energy of the exp(), 1/r\^6, and 1/r\^8 portion of the pair
interaction.

The *born/coul/long* pair style supports the
:doc:`pair_modify <pair_modify>` table option to tabulate the
short-range portion of the long-range Coulombic interaction.

These styles support the pair\_modify tail option for adding long-range
tail corrections to energy and pressure.

Thess styles writes thei information to binary :doc:`restart <restart>`
files, so pair\_style and pair\_coeff commands do not need to be
specified in an input script that reads a restart file.

These styles can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command.  They do not support the *inner*\ ,
*middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


The *born/coul/long* style is part of the KSPACE package.  It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style buck <pair_buck>`

**Default:** none


----------


.. _FumiTosi:



Fumi and Tosi, J Phys Chem Solids, 25, 31 (1964),
Fumi and Tosi, J Phys Chem Solids, 25, 45 (1964).
