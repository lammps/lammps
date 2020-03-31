.. index:: pair_style buck/long/coul/long

pair_style buck/long/coul/long command
=======================================

pair_style buck/long/coul/long/omp command
===========================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style buck/long/coul/long flag_buck flag_coul cutoff (cutoff2)

* flag_buck = *long* or *cut*

  .. parsed-literal::

       *long* = use Kspace long-range summation for the dispersion term 1/r\^6
       *cut* = use a cutoff

* flag_coul = *long* or *off*

  .. parsed-literal::

       *long* = use Kspace long-range summation for the Coulombic term 1/r
       *off* = omit the Coulombic term

* cutoff = global cutoff for Buckingham (and Coulombic if only 1 cutoff) (distance units)
* cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style buck/long/coul/long cut off 2.5
   pair_style buck/long/coul/long cut long 2.5 4.0
   pair_style buck/long/coul/long long long 4.0
   pair_coeff * * 1 1
   pair_coeff 1 1 1 3 4

Description
"""""""""""

The *buck/long/coul/long* style computes a Buckingham potential (exp/6
instead of Lennard-Jones 12/6) and Coulombic potential, given by

.. math::

   E = & A e^{-r / \rho} - \frac{C}{r^6} \qquad r < r_c \\
   E = & \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c

:math:`r_c` is the cutoff.  If one cutoff is specified in the pair_style
command, it is used for both the Buckingham and Coulombic terms.  If
two cutoffs are specified, they are used as cutoffs for the Buckingham
and Coulombic terms respectively.

The purpose of this pair style is to capture long-range interactions
resulting from both attractive 1/r\^6 Buckingham and Coulombic 1/r
interactions.  This is done by use of the *flag_buck* and *flag_coul*
settings.  The :ref:`Ismail <Ismail>` paper has more details on when it is
appropriate to include long-range 1/r\^6 interactions, using this
potential.

If *flag_buck* is set to *long*\ , no cutoff is used on the Buckingham
1/r\^6 dispersion term.  The long-range portion can be calculated by
using the :doc:`kspace_style ewald/disp or pppm/disp <kspace_style>`
commands.  The specified Buckingham cutoff then determines which
portion of the Buckingham interactions are computed directly by the
pair potential versus which part is computed in reciprocal space via
the Kspace style.  If *flag_buck* is set to *cut*\ , the Buckingham
interactions are simply cutoff, as with :doc:`pair_style buck <pair_buck>`.

If *flag_coul* is set to *long*\ , no cutoff is used on the Coulombic
interactions.  The long-range portion can calculated by using any of
several :doc:`kspace_style <kspace_style>` command options such as
*pppm* or *ewald*\ .  Note that if *flag_buck* is also set to long, then
the *ewald/disp* or *pppm/disp* Kspace style needs to be used to
perform the long-range calculations for both the Buckingham and
Coulombic interactions.  If *flag_coul* is set to *off*\ , Coulombic
interactions are not computed.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* A (energy units)
* rho (distance units)
* C (energy-distance\^6 units)
* cutoff (distance units)
* cutoff2 (distance units)

The second coefficient, rho, must be greater than zero.

The latter 2 coefficients are optional.  If not specified, the global
Buckingham and Coulombic cutoffs specified in the pair_style command
are used.  If only one cutoff is specified, it is used as the cutoff
for both Buckingham and Coulombic interactions for this type pair.  If
both coefficients are specified, they are used as the Buckingham and
Coulombic cutoffs for this type pair.  Note that if you are using
*flag_buck* set to *long*\ , you cannot specify a Buckingham cutoff for
an atom type pair, since only one global Buckingham cutoff is allowed.
Similarly, if you are using *flag_coul* set to *long*\ , you cannot
specify a Coulombic cutoff for an atom type pair, since only one
global Coulombic cutoff is allowed.

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

This pair styles does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the exp() and 1/r\^6 portion of the pair
interaction, assuming *flag_buck* is *cut*\ .

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the Buckingham portion of the pair
interaction.

This pair style supports the :doc:`pair_modify <pair_modify>` table and
table/disp options since they can tabulate the short-range portion of
the long-range Coulombic and dispersion interactions.

This pair style write its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style supports the use of the *inner*\ , *middle*\ , and *outer*
keywords of the :doc:`run_style respa <run_style>` command, meaning the
pairwise forces can be partitioned by distance at different levels of
the rRESPA hierarchy.  See the :doc:`run_style <run_style>` command for
details.

----------

Restrictions
""""""""""""

This style is part of the KSPACE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none

----------

.. _Ismail:

**(Ismail)** Ismail, Tsige, In 't Veld, Grest, Molecular Physics
(accepted) (2007).
