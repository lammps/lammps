.. index:: pair_style lj/spica
.. index:: pair_style lj/spica/gpu
.. index:: pair_style lj/spica/kk
.. index:: pair_style lj/spica/omp
.. index:: pair_style lj/spica/coul/long
.. index:: pair_style lj/spica/coul/long/gpu
.. index:: pair_style lj/spica/coul/long/omp
.. index:: pair_style lj/spica/coul/msm
.. index:: pair_style lj/spica/coul/msm/omp

pair_style lj/spica command
===========================

Accelerator Variants: *lj/spica/gpu*, *lj/spica/kk*, *lj/spica/omp*

pair_style lj/spica/coul/long command
=====================================

Accelerator Variants: *lj/spica/coul/long/gpu*, *lj/spica/coul/long/omp*

pair_style lj/spica/coul/msm command
====================================

Accelerator Variants: *lj/spica/coul/msm/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/spica* or *lj/spica/coul/long*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/spica* args = cutoff
       cutoff = global cutoff for Lennard Jones interactions (distance units)
     *lj/spica/coul/long* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/spica 2.5
   pair_coeff 1 1 lj12_6 1 1.1 2.8

   pair_style lj/spica/coul/long 10.0
   pair_style lj/spica/coul/long 10.0 12.0
   pair_coeff 1 1 lj9_6 100.0 3.5 12.0

   pair_style lj/spica/coul/msm 10.0
   pair_style lj/spica/coul/msm 10.0 12.0
   pair_coeff 1 1 lj9_6 100.0 3.5 12.0

Description
"""""""""""

The *lj/spica* styles compute a 9/6, 12/4, 12/5, or 12/6 Lennard-Jones potential,
given by

.. math::

   E = & \frac{27}{4} \epsilon \left[ \left(\frac{\sigma}{r}\right)^{9} -
                         \left(\frac{\sigma}{r}\right)^6 \right]
                         \qquad r < r_c \\
   E = & \frac{3\sqrt{3}}{2} \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                         \left(\frac{\sigma}{r}\right)^4 \right]
                         \qquad r < r_c \\
   E = & \frac{12}{7}\left(\frac{12}{5}\right)^{\left(\frac{5}{7}\right)} \epsilon
                         \left[ \left(\frac{\sigma}{r}\right)^{12} -
                         \left(\frac{\sigma}{r}\right)^5 \right]
                         \qquad r < r_c \\
   E = &  4 \epsilon  \left[ \left(\frac{\sigma}{r}\right)^{12} -
                         \left(\frac{\sigma}{r}\right)^6 \right]
                         \qquad r < r_c

as required for the SPICA (formerly called SDK) and the pSPICA Coarse-grained MD parameterization discussed in
:ref:`(Shinoda) <Shinoda3>`, :ref:`(DeVane) <DeVane>`, :ref:`(Seo) <Seo>`, and :ref:`(Miyazaki) <Miyazaki>`.
Rc is the cutoff.
Summary information on these force fields can be found at https://www.spica-ff.org

Style *lj/spica/coul/long* computes the adds Coulombic interactions
with an additional damping factor applied so it can be used in
conjunction with the :doc:`kspace_style <kspace_style>` command and
its *ewald* or *pppm* or *pppm/cg* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within
this distance are computed directly; interactions outside that
distance are computed in reciprocal space.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cg_type (lj9_6, lj12_4, lj12_5, or lj12_6)
* epsilon (energy units)
* sigma (distance units)
* cutoff1 (distance units)

Note that sigma is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum. The prefactors
are chosen so that the potential minimum is at -epsilon.

The latter 2 coefficients are optional.  If not specified, the global
LJ and Coulombic cutoffs specified in the pair_style command are used.
If only one cutoff is specified, it is used as the cutoff for both LJ
and Coulombic interactions for this type pair.  If both coefficients
are specified, they are used as the LJ and Coulombic cutoffs for this
type pair.

For *lj/spica/coul/long* and *lj/spica/coul/msm* only the LJ cutoff can be
specified since a Coulombic cutoff cannot be specified for an
individual I,J type pair.  All type pairs use the same global
Coulombic cutoff specified in the pair_style command.

The original implementation of the above styles are
style *lj/sdk*, *lj/sdk/coul/long*, and *lj/sdk/coul/msm*,
and available for backward compatibility.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/spica pair styles *cannot* be mixed,
since different pairs may have different exponents. So all parameters
for all pairs have to be specified explicitly through the "pair_coeff"
command. Defining then in a data file is also not supported, due to
limitations of that file format.

All of the lj/spica pair styles support the
:doc:`pair_modify <pair_modify>` shift option for the energy of the
Lennard-Jones portion of the pair interaction.

The *lj/spica/coul/long* pair styles support the
:doc:`pair_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the lj/spica pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do
not need to be specified in an input script that reads a restart file.

The lj/spica and lj/cut/coul/long pair styles do not support
the use of the *inner*, *middle*, and *outer* keywords of the :doc:`run_style respa <run_style>` command.

----------

Restrictions
""""""""""""

All of the lj/spica pair styles are part of the CG-SPICA package.  The
*lj/spica/coul/long* style also requires the KSPACE package to be built
(which is enabled by default).  They are only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`angle_style spica <angle_spica>`

Default
"""""""

none

----------

.. _Shinoda3:

**(Shinoda)** Shinoda, DeVane, Klein, Mol Sim, 33, 27-36 (2007).

.. _DeVane:

**(DeVane)**  Shinoda, DeVane, Klein, Soft Matter, 4, 2453-2462 (2008).

.. _Seo:

**(Seo)**  Seo, Shinoda, J Chem Theory Comput, 15, 762-774 (2019).

.. _Miyazaki:

**(Miyazaki)**  Miyazaki, Okazaki, Shinoda, J Chem Theory Comput, 16, 782-793 (2020).
