.. index:: pair_style nm/cut

pair_style nm/cut command
=========================

pair_style nm/cut/coul/cut command
==================================

pair_style nm/cut/coul/long command
===================================

pair_style nm/cut/omp command
=============================

pair_style nm/cut/coul/cut/omp command
======================================

pair_style nm/cut/coul/long/omp command
=======================================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style style args

* style = *nm/cut* or *nm/cut/coul/cut* or *nm/cut/coul/long*
* args = list of arguments for a particular style
  
  .. parsed-literal::
  
       *nm/cut* args = cutoff
         cutoff = global cutoff for Pair interactions (distance units)
       *nm/cut/coul/cut* args = cutoff (cutoff2)
         cutoff = global cutoff for Pair (and Coulombic if only 1 arg) (distance units)
         cutoff2 = global cutoff for Coulombic (optional) (distance units)
       *nm/cut/coul/long* args = cutoff (cutoff2)
         cutoff = global cutoff for Pair (and Coulombic if only 1 arg) (distance units)
         cutoff2 = global cutoff for Coulombic (optional) (distance units)



Examples
""""""""


.. code-block:: LAMMPS

   pair_style nm/cut 12.0
   pair_coeff * * 0.01 5.4 8.0 7.0
   pair_coeff 1 1 0.01 4.4 7.0 6.0

   pair_style nm/cut/coul/cut 12.0 15.0
   pair_coeff * * 0.01 5.4 8.0 7.0
   pair_coeff 1 1 0.01 4.4 7.0 6.0

   pair_style nm/cut/coul/long 12.0 15.0
   pair_coeff * * 0.01 5.4 8.0 7.0
   pair_coeff 1 1 0.01 4.4 7.0 6.0

Description
"""""""""""

Style *nm* computes site-site interactions based on the N-M potential
by :ref:`Clarke <Clarke>`, mainly used for ionic liquids.  A site can
represent a single atom or a united-atom site.  The energy of an
interaction has the following form:

.. math::

   E = \frac{E_0}{(n-m)} \left[ m \left(\frac{r_0}{r}\right)^n - n
   \left(\frac{r_0}{r}\right)^m \right] \qquad r < r_c

where :math:`r_c` is the cutoff.

Style *nm/cut/coul/cut* adds a Coulombic pairwise interaction given by

.. math::

   E = \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c

where :math:`C` is an energy-conversion constant, :math:`q_i` and :math:`q_j`
are the charges on the 2 atoms, and epsilon is the dielectric constant which can
be set by the :doc:`dielectric <dielectric>` command.  If one cutoff is
specified in the pair\_style command, it is used for both the N-M and Coulombic
terms.  If two cutoffs are specified, they are used as cutoffs for the N-M and
Coulombic terms respectively.

Styles *nm/cut/coul/long* compute the same
Coulombic interactions as style *nm/cut/coul/cut* except that an
additional damping factor is applied to the Coulombic term so it can
be used in conjunction with the :doc:`kspace_style <kspace_style>`
command and its *ewald* or *pppm* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within this
distance are computed directly; interactions outside that distance are
computed in reciprocal space.

For all of the *nm* pair styles, the following coefficients must
be defined for each pair of atoms types
via the :doc:`pair_coeff <pair_coeff>` command as in the
examples above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.

* :math:`E_0` (energy units)
* :math:`r_0` (distance units)
* :math:`n` (unitless)
* :math:`m` (unitless)
* cutoff1 (distance units)
* cutoff2 (distance units)

The latter 2 coefficients are optional.  If not specified, the global
N-M and Coulombic cutoffs specified in the pair\_style command are used.
If only one cutoff is specified, it is used as the cutoff for both N-M
and Coulombic interactions for this type pair.  If both coefficients
are specified, they are used as the N-M and Coulombic cutoffs for this
type pair.  You cannot specify 2 cutoffs for style *nm*\ , since it
has no Coulombic terms.

For *nm/cut/coul/long* only the N-M cutoff can be specified since a
Coulombic cutoff cannot be specified for an individual I,J type pair.
All type pairs use the same global Coulombic cutoff specified in the
pair\_style command.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

These pair styles do not support mixing. Thus, coefficients for all
I,J pairs must be specified explicitly.

All of the *nm* pair styles supports the
:doc:`pair_modify <pair_modify>` shift option for the energy of the pair
interaction.

The *nm/cut/coul/long* pair styles support the
:doc:`pair_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the *nm* pair styles support the :doc:`pair_modify <pair_modify>`
tail option for adding a long-range tail correction to the energy and
pressure for the N-M portion of the pair interaction.

All of the *nm* pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

All of the *nm* pair styles can only be used via the *pair* keyword of
the :doc:`run_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.


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

Restrictions
""""""""""""


These pair styles are part of the MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


----------


.. _Clarke:



**(Clarke)** Clarke and Smith, J Chem Phys, 84, 2290 (1986).
