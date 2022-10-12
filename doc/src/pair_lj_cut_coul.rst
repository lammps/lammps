.. index:: pair_style lj/cut/coul/cut
.. index:: pair_style lj/cut/coul/cut/gpu
.. index:: pair_style lj/cut/coul/cut/kk
.. index:: pair_style lj/cut/coul/cut/omp
.. index:: pair_style lj/cut/coul/debye
.. index:: pair_style lj/cut/coul/debye/gpu
.. index:: pair_style lj/cut/coul/debye/kk
.. index:: pair_style lj/cut/coul/debye/omp
.. index:: pair_style lj/cut/coul/dsf
.. index:: pair_style lj/cut/coul/dsf/gpu
.. index:: pair_style lj/cut/coul/dsf/kk
.. index:: pair_style lj/cut/coul/dsf/omp
.. index:: pair_style lj/cut/coul/long
.. index:: pair_style lj/cut/coul/long/gpu
.. index:: pair_style lj/cut/coul/long/kk
.. index:: pair_style lj/cut/coul/long/intel
.. index:: pair_style lj/cut/coul/long/opt
.. index:: pair_style lj/cut/coul/long/omp
.. index:: pair_style lj/cut/coul/msm
.. index:: pair_style lj/cut/coul/msm/gpu
.. index:: pair_style lj/cut/coul/msm/omp
.. index:: pair_style lj/cut/coul/wolf
.. index:: pair_style lj/cut/coul/wolf/omp

pair_style lj/cut/coul/cut command
==================================

Accelerator Variants: *lj/cut/coul/cut/gpu*, *lj/cut/coul/cut/kk*, *lj/cut/coul/cut/omp*

pair_style lj/cut/coul/debye command
====================================

Accelerator Variants: *lj/cut/coul/debye/gpu*, *lj/cut/coul/debye/kk*, *lj/cut/coul/debye/omp*

pair_style lj/cut/coul/dsf command
==================================

Accelerator Variants: *lj/cut/coul/dsf/gpu*, *lj/cut/coul/dsf/kk*, *lj/cut/coul/dsf/omp*

pair_style lj/cut/coul/long command
===================================

Accelerator Variants: *lj/cut/coul/long/gpu*, *lj/cut/coul/long/kk*, *lj/cut/coul/long/intel*, *lj/cut/coul/long/opt*, *lj/cut/coul/long/omp*

pair_style lj/cut/coul/msm command
==================================

Accelerator Variants: *lj/cut/coul/msm/gpu*, *lj/cut/coul/msm/omp*

pair_style lj/cut/coul/wolf command
===================================

Accelerator Variants: *lj/cut/coul/wolf/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/cut/coul/cut* or *lj/cut/coul/debye* or *lj/cut/coul/dsf* or *lj/cut/coul/long* *lj/cut/coul/msm* or *lj/cut/coul/wolf*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/cut/coul/cut* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/coul/debye* args = kappa cutoff (cutoff2)
       kappa = inverse of the Debye length (inverse distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/coul/dsf* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (distance units)
     *lj/cut/coul/long* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/coul/msm* args = cutoff (cutoff2)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/coul/wolf* args = alpha cutoff (cutoff2)
       alpha = damping parameter (inverse distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 2 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/cut/coul/cut 10.0
   pair_style lj/cut/coul/cut 10.0 8.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0
   pair_coeff 1 1 100.0 3.5 9.0 9.0

   pair_style lj/cut/coul/debye 1.5 3.0
   pair_style lj/cut/coul/debye 1.5 2.5 5.0
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.5 2.5
   pair_coeff 1 1 1.0 1.5 2.5 5.0

   pair_style lj/cut/coul/dsf 0.05 2.5 10.0
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.0 2.5

   pair_style lj/cut/coul/long 10.0
   pair_style lj/cut/coul/long 10.0 8.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/coul/msm 10.0
   pair_style lj/cut/coul/msm 10.0 8.0
   pair_coeff * * 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/coul/wolf 0.2 5. 10.0
   pair_coeff * * 1.0 1.0
   pair_coeff 1 1 1.0 1.0 2.5

Description
"""""""""""

The *lj/cut/coul* styles compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
       \left(\frac{\sigma}{r}\right)^6 \right]
                       \qquad r < r_c

:math:`r_c` is the cutoff.

Style *lj/cut/coul/cut* adds a Coulombic pairwise interaction given by

.. math::

   E = \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c

where :math:`C` is an energy-conversion constant, :math:`q_i` and :math:`q_j`
are the charges on the 2 atoms, and :math:`\epsilon` is the dielectric
constant which can be set by the :doc:`dielectric <dielectric>` command.
If one cutoff is specified in the pair_style command, it is used for
both the LJ and Coulombic terms.  If two cutoffs are specified, they are
used as cutoffs for the LJ and Coulombic terms respectively.

Style *lj/cut/coul/debye* adds an additional exp() damping factor
to the Coulombic term, given by

.. math::

   E = \frac{C q_i q_j}{\epsilon  r} \exp(- \kappa r) \qquad r < r_c

where :math:`\kappa` is the inverse of the Debye length.  This potential
is another way to mimic the screening effect of a polar solvent.

Style *lj/cut/coul/dsf* computes the Coulombic term via the damped
shifted force model described in :ref:`Fennell <Fennell2>`, given by:

.. math::

   E =
    q_iq_j \left[ \frac{\mbox{erfc} (\alpha r)}{r} -  \frac{\mbox{erfc} (\alpha r_c)}{r_c} +
   \left( \frac{\mbox{erfc} (\alpha r_c)}{r_c^2} +  \frac{2\alpha}{\sqrt{\pi}}\frac{\exp (-\alpha^2    r^2_c)}{r_c} \right)(r-r_c) \right] \qquad r < r_c

where :math:`\alpha` is the damping parameter and erfc() is the complementary
error-function. This potential is essentially a short-range,
spherically-truncated, charge-neutralized, shifted, pairwise *1/r*
summation.  The potential is based on Wolf summation, proposed as an
alternative to Ewald summation for condensed phase systems where
charge screening causes electrostatic interactions to become
effectively short-ranged. In order for the electrostatic sum to be
absolutely convergent, charge neutralization within the cutoff radius
is enforced by shifting the potential through placement of image
charges on the cutoff sphere. Convergence can often be improved by
setting :math:`\alpha` to a small non-zero value.

Styles *lj/cut/coul/long* and *lj/cut/coul/msm* compute the same
Coulombic interactions as style *lj/cut/coul/cut* except that an
additional damping factor is applied to the Coulombic term so it can
be used in conjunction with the :doc:`kspace_style <kspace_style>`
command and its *ewald* or *pppm* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within this
distance are computed directly; interactions outside that distance are
computed in reciprocal space.

Style *coul/wolf* adds a Coulombic pairwise interaction via the Wolf
summation method, described in :ref:`Wolf <Wolf3>`, given by:

.. math::

   E_i = \frac{1}{2} \sum_{j \neq i}
   \frac{q_i q_j {\rm erfc}(\alpha r_{ij})}{r_{ij}} +
   \frac{1}{2} \sum_{j \neq i}
   \frac{q_i q_j {\rm erf}(\alpha r_{ij})}{r_{ij}} \qquad r < r_c

where :math:`\alpha` is the damping parameter, and erfc() is the
complementary error-function terms.  This potential is essentially a
short-range, spherically-truncated, charge-neutralized, shifted,
pairwise *1/r* summation.  With a manipulation of adding and subtracting
a self term (for i = j) to the first and second term on the
right-hand-side, respectively, and a small enough :math:`\alpha` damping
parameter, the second term shrinks and the potential becomes a
rapidly-converging real-space summation.  With a long enough cutoff and
small enough :math:`\alpha` parameter, the energy and forces calculated by the
Wolf summation method approach those of the Ewald sum.  So it is a means
of getting effective long-range interactions with a short-range
potential.

Coefficients
""""""""""""

For all of the *lj/cut/coul* pair styles, the following coefficients must
be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands, or by mixing as
described below:

* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* cutoff1 (distance units)
* cutoff2 (distance units)

Note that :math:`\sigma` is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at :math:`2^{\frac{1}{6}} \sigma`.

The latter 2 coefficients are optional.  If not specified, the global
LJ and Coulombic cutoffs specified in the pair_style command are used.
If only one cutoff is specified, it is used as the cutoff for both LJ
and Coulombic interactions for this type pair.  If both coefficients
are specified, they are used as the LJ and Coulombic cutoffs for this
type pair.

For *lj/cut/coul/long* and *lj/cut/coul/msm* only the LJ cutoff can be
specified since a Coulombic cutoff cannot be specified for an individual I,J
type pair.  All type pairs use the same global Coulombic cutoff specified in
the pair_style command.

----------

A version of these styles with a soft core, *lj/cut/coul/soft*\ and
*lj/cut/coul/long/soft*, suitable for use in free energy calculations, is
part of the FEP package and is documented with the :doc:`pair_style */soft <pair_fep_soft>` styles.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed.
The default mix value is *geometric*\ .  See the "pair_modify" command
for details.

All of the *lj/cut* pair styles support the
:doc:`pair_modify <pair_modify>` shift option for the energy of the
Lennard-Jones portion of the pair interaction.

The *lj/cut/coul/long* pair styles support the
:doc:`pair_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the *lj/cut* pair styles support the
:doc:`pair_modify <pair_modify>` tail option for adding a long-range
tail correction to the energy and pressure for the Lennard-Jones
portion of the pair interaction.

All of the *lj/cut* pair styles write their information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do
not need to be specified in an input script that reads a restart file.

The *lj/cut/coul/long* pair styles support the use of the
*inner*, *middle*, and *outer* keywords of the :doc:`run_style respa <run_style>` command, meaning the pairwise forces can be
partitioned by distance at different levels of the rRESPA hierarchy.
The other styles only support the *pair* keyword of run_style respa.
See the :doc:`run_style <run_style>` command for details.

----------

Restrictions
""""""""""""

The *lj/cut/coul/long* and *lj/cut/coul/msm* styles are part of the KSPACE package.

The *lj/cut/coul/debye*, *lj/cut/coul/dsf*, and *lj/cut/coul/wolf* styles are part
of the EXTRA-PAIR package.

These styles are only enabled if LAMMPS was built with those respective
packages.  See the :doc:`Build package <Build_package>` page for
more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Wolf3:

**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J Chem
Phys, 110, 8254 (1999).

.. _Fennell2:

**(Fennell)** C. J. Fennell, J. D. Gezelter, J Chem Phys, 124,
234104 (2006).
