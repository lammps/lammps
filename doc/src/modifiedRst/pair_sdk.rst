.. index:: pair\_style lj/sdk

pair\_style lj/sdk command
==========================

pair\_style lj/sdk/gpu command
==============================

pair\_style lj/sdk/kk command
=============================

pair\_style lj/sdk/omp command
==============================

pair\_style lj/sdk/coul/long command
====================================

pair\_style lj/sdk/coul/long/gpu command
========================================

pair\_style lj/sdk/coul/long/omp command
========================================

pair\_style lj/sdk/coul/msm command
===================================

pair\_style lj/sdk/coul/msm/omp command
=======================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *lj/sdk* or *lj/sdk/coul/long*
* args = list of arguments for a particular style


.. parsed-literal::

     *lj/sdk* args = cutoff
       cutoff = global cutoff for Lennard Jones interactions (distance units)
     *lj/sdk/coul/long* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style lj/sdk 2.5
   pair_coeff 1 1 lj12_6 1 1.1 2.8

   pair_style lj/sdk/coul/long 10.0
   pair_style lj/sdk/coul/long 10.0 12.0
   pair_coeff 1 1 lj9_6 100.0 3.5 12.0

   pair_style lj/sdk/coul/msm 10.0
   pair_style lj/sdk/coul/msm 10.0 12.0
   pair_coeff 1 1 lj9_6 100.0 3.5 12.0

Description
"""""""""""

The *lj/sdk* styles compute a 9/6, 12/4, or 12/6 Lennard-Jones potential,
given by

.. math::

 E = & \frac{27}{4} \epsilon \left[ \left(\frac{\sigma}{r}\right)^{9} - 
                       \left(\frac{\sigma}{r}\right)^6 \right] &
                       \qquad r < r_c \\
 E = & \frac{3\sqrt{3}}{2} \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - 
                       \left(\frac{\sigma}{r}\right)^4 \right] &
                       \qquad r < r_c \\
 E = &  4 \epsilon  \left[ \left(\frac{\sigma}{r}\right)^{12} - 
                       \left(\frac{\sigma}{r}\right)^6 \right] &
                       \qquad r < r_c


as required for the SDK Coarse-grained MD parameterization discussed in
:ref:`(Shinoda) <Shinoda3>` and :ref:`(DeVane) <DeVane>`.  Rc is the cutoff.

Style *lj/sdk/coul/long* computes the adds Coulombic interactions
with an additional damping factor applied so it can be used in
conjunction with the :doc:`kspace\_style <kspace_style>` command and
its *ewald* or *pppm* or *pppm/cg* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within
this distance are computed directly; interactions outside that
distance are computed in reciprocal space.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands, or by mixing as described below:

* cg\_type (lj9\_6, lj12\_4, or lj12\_6)
* epsilon (energy units)
* sigma (distance units)
* cutoff1 (distance units)

Note that sigma is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum. The prefactors
are chosen so that the potential minimum is at -epsilon.

The latter 2 coefficients are optional.  If not specified, the global
LJ and Coulombic cutoffs specified in the pair\_style command are used.
If only one cutoff is specified, it is used as the cutoff for both LJ
and Coulombic interactions for this type pair.  If both coefficients
are specified, they are used as the LJ and Coulombic cutoffs for this
type pair.

For *lj/sdk/coul/long* and *lj/sdk/coul/msm* only the LJ cutoff can be
specified since a Coulombic cutoff cannot be specified for an
individual I,J type pair.  All type pairs use the same global
Coulombic cutoff specified in the pair\_style command.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp* or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP, and OPT packages respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


**Mixing, shift, table, tail correction, restart, and rRESPA info**\ :

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/sdk pair styles *cannot* be mixed,
since different pairs may have different exponents. So all parameters
for all pairs have to be specified explicitly through the "pair\_coeff"
command. Defining then in a data file is also not supported, due to
limitations of that file format.

All of the lj/sdk pair styles support the
:doc:`pair\_modify <pair_modify>` shift option for the energy of the
Lennard-Jones portion of the pair interaction.

The *lj/sdk/coul/long* pair styles support the
:doc:`pair\_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the lj/sdk pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do
not need to be specified in an input script that reads a restart file.

The lj/sdk and lj/cut/coul/long pair styles do not support
the use of the *inner*\ , *middle*\ , and *outer* keywords of the :doc:`run\_style respa <run_style>` command.


----------


Restrictions
""""""""""""


All of the lj/sdk pair styles are part of the USER-CGSDK package.  The
*lj/sdk/coul/long* style also requires the KSPACE package to be built
(which is enabled by default).  They are only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`angle\_style sdk <angle_sdk>`

**Default:** none


----------


.. _Shinoda3:



**(Shinoda)** Shinoda, DeVane, Klein, Mol Sim, 33, 27 (2007).

.. _DeVane:



**(DeVane)**  Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
