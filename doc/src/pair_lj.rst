.. index:: pair\_style lj/cut

pair\_style lj/cut command
==========================

pair\_style lj/cut/gpu command
==============================

pair\_style lj/cut/intel command
================================

pair\_style lj/cut/kk command
=============================

pair\_style lj/cut/opt command
==============================

pair\_style lj/cut/omp command
==============================

pair\_style lj/cut/coul/cut command
===================================

pair\_style lj/cut/coul/cut/gpu command
=======================================

pair\_style lj/cut/coul/cut/kk command
======================================

pair\_style lj/cut/coul/cut/omp command
=======================================

pair\_style lj/cut/coul/debye command
=====================================

pair\_style lj/cut/coul/debye/gpu command
=========================================

pair\_style lj/cut/coul/debye/kk command
========================================

pair\_style lj/cut/coul/debye/omp command
=========================================

pair\_style lj/cut/coul/dsf command
===================================

pair\_style lj/cut/coul/dsf/gpu command
=======================================

pair\_style lj/cut/coul/dsf/kk command
======================================

pair\_style lj/cut/coul/dsf/omp command
=======================================

pair\_style lj/cut/coul/long command
====================================

pair\_style lj/cut/coul/long/gpu command
========================================

pair\_style lj/cut/coul/long/kk command
=======================================

pair\_style lj/cut/coul/long/intel command
==========================================

pair\_style lj/cut/coul/long/opt command
========================================

pair\_style lj/cut/coul/long/omp command
========================================

pair\_style lj/cut/coul/msm command
===================================

pair\_style lj/cut/coul/msm/gpu command
=======================================

pair\_style lj/cut/coul/msm/omp command
=======================================

pair\_style lj/cut/coul/wolf command
====================================

pair\_style lj/cut/coul/wolf/omp command
========================================

pair\_style lj/cut/tip4p/cut command
====================================

pair\_style lj/cut/tip4p/cut/omp command
========================================

pair\_style lj/cut/tip4p/long command
=====================================

pair\_style lj/cut/tip4p/long/gpu command
========================================

pair\_style lj/cut/tip4p/long/omp command
=========================================

pair\_style lj/cut/tip4p/long/opt command
=========================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *lj/cut* or *lj/cut/coul/cut* or *lj/cut/coul/debye* or *lj/cut/coul/dsf* or *lj/cut/coul/long* *lj/cut/coul/msm* or *lj/cut/tip4p/long*
* args = list of arguments for a particular style


.. parsed-literal::

     *lj/cut* args = cutoff
       cutoff = global cutoff for Lennard Jones interactions (distance units)
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
     *lj/cut/tip4p/cut* args = otype htype btype atype qdist cutoff (cutoff2)
       otype,htype = atom types for TIP4P O and H
       btype,atype = bond and angle types for TIP4P waters
       qdist = distance from O atom to massless charge (distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)
     *lj/cut/tip4p/long* args = otype htype btype atype qdist cutoff (cutoff2)
       otype,htype = atom types for TIP4P O and H
       btype,atype = bond and angle types for TIP4P waters
       qdist = distance from O atom to massless charge (distance units)
       cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
       cutoff2 = global cutoff for Coulombic (optional) (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style lj/cut 2.5
   pair_coeff \* \* 1 1
   pair_coeff 1 1 1 1.1 2.8

   pair_style lj/cut/coul/cut 10.0
   pair_style lj/cut/coul/cut 10.0 8.0
   pair_coeff \* \* 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0
   pair_coeff 1 1 100.0 3.5 9.0 9.0

   pair_style lj/cut/coul/debye 1.5 3.0
   pair_style lj/cut/coul/debye 1.5 2.5 5.0
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 1.0 1.5 2.5
   pair_coeff 1 1 1.0 1.5 2.5 5.0

   pair_style lj/cut/coul/dsf 0.05 2.5 10.0
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 1.0 1.0 2.5

   pair_style lj/cut/coul/long 10.0
   pair_style lj/cut/coul/long 10.0 8.0
   pair_coeff \* \* 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/coul/msm 10.0
   pair_style lj/cut/coul/msm 10.0 8.0
   pair_coeff \* \* 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/tip4p/cut 1 2 7 8 0.15 12.0
   pair_style lj/cut/tip4p/cut 1 2 7 8 0.15 12.0 10.0
   pair_coeff \* \* 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

   pair_style lj/cut/coul/wolf 0.2 5. 10.0
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 1.0 1.0 2.5

   pair_style lj/cut/tip4p/long 1 2 7 8 0.15 12.0
   pair_style lj/cut/tip4p/long 1 2 7 8 0.15 12.0 10.0
   pair_coeff \* \* 100.0 3.0
   pair_coeff 1 1 100.0 3.5 9.0

Description
"""""""""""

The *lj/cut* styles compute the standard 12/6 Lennard-Jones potential,
given by

.. image:: Eqs/pair_lj.jpg
   :align: center

Rc is the cutoff.

Style *lj/cut/coul/cut* adds a Coulombic pairwise interaction given by

.. image:: Eqs/pair_coulomb.jpg
   :align: center

where C is an energy-conversion constant, Qi and Qj are the charges on
the 2 atoms, and epsilon is the dielectric constant which can be set
by the :doc:`dielectric <dielectric>` command.  If one cutoff is
specified in the pair\_style command, it is used for both the LJ and
Coulombic terms.  If two cutoffs are specified, they are used as
cutoffs for the LJ and Coulombic terms respectively.

Style *lj/cut/coul/debye* adds an additional exp() damping factor
to the Coulombic term, given by

.. image:: Eqs/pair_debye.jpg
   :align: center

where kappa is the inverse of the Debye length.  This potential is
another way to mimic the screening effect of a polar solvent.

Style *lj/cut/coul/dsf* computes the Coulombic term via the damped
shifted force model described in :ref:`Fennell <Fennell2>`, given by:

.. image:: Eqs/pair_coul_dsf.jpg
   :align: center

where *alpha* is the damping parameter and erfc() is the complementary
error-function. This potential is essentially a short-range,
spherically-truncated, charge-neutralized, shifted, pairwise *1/r*
summation.  The potential is based on Wolf summation, proposed as an
alternative to Ewald summation for condensed phase systems where
charge screening causes electrostatic interactions to become
effectively short-ranged. In order for the electrostatic sum to be
absolutely convergent, charge neutralization within the cutoff radius
is enforced by shifting the potential through placement of image
charges on the cutoff sphere. Convergence can often be improved by
setting *alpha* to a small non-zero value.

Styles *lj/cut/coul/long* and *lj/cut/coul/msm* compute the same
Coulombic interactions as style *lj/cut/coul/cut* except that an
additional damping factor is applied to the Coulombic term so it can
be used in conjunction with the :doc:`kspace\_style <kspace_style>`
command and its *ewald* or *pppm* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within this
distance are computed directly; interactions outside that distance are
computed in reciprocal space.

Style *coul/wolf* adds a Coulombic pairwise interaction via the Wolf
summation method, described in :ref:`Wolf <Wolf1>`, given by:

.. image:: Eqs/pair_coul_wolf.jpg
   :align: center

where *alpha* is the damping parameter, and erfc() is the
complementary error-function terms.  This potential
is essentially a short-range, spherically-truncated,
charge-neutralized, shifted, pairwise *1/r* summation.  With a
manipulation of adding and subtracting a self term (for i = j) to the
first and second term on the right-hand-side, respectively, and a
small enough *alpha* damping parameter, the second term shrinks and
the potential becomes a rapidly-converging real-space summation.  With
a long enough cutoff and small enough alpha parameter, the energy and
forces calculated by the Wolf summation method approach those of the
Ewald sum.  So it is a means of getting effective long-range
interactions with a short-range potential.

Styles *lj/cut/tip4p/cut* and *lj/cut/tip4p/long* implement the TIP4P
water model of :ref:`(Jorgensen) <Jorgensen2>`, which introduces a massless
site located a short distance away from the oxygen atom along the
bisector of the HOH angle.  The atomic types of the oxygen and
hydrogen atoms, the bond and angle types for OH and HOH interactions,
and the distance to the massless charge site are specified as
pair\_style arguments.  Style *lj/cut/tip4p/cut* uses a cutoff for
Coulomb interactions; style *lj/cut/tip4p/long* is for use with a
long-range Coulombic solver (Ewald or PPPM).

.. note::

   For each TIP4P water molecule in your system, the atom IDs for
   the O and 2 H atoms must be consecutive, with the O atom first.  This
   is to enable LAMMPS to "find" the 2 H atoms associated with each O
   atom.  For example, if the atom ID of an O atom in a TIP4P water
   molecule is 500, then its 2 H atoms must have IDs 501 and 502.

See the :doc:`Howto tip4p <Howto_tip4p>` doc page for more information
on how to use the TIP4P pair styles and lists of parameters to set.
Note that the neighbor list cutoff for Coulomb interactions is
effectively extended by a distance 2\*qdist when using the TIP4P pair
style, to account for the offset distance of the fictitious charges on
O atoms in water molecules.  Thus it is typically best in an
efficiency sense to use a LJ cutoff >= Coulombic cutoff + 2\*qdist, to
shrink the size of the neighbor list.  This leads to slightly larger
cost for the long-range calculation, so you can test the trade-off for
your model.

For all of the *lj/cut* pair styles, the following coefficients must
be defined for each pair of atoms types via the
:doc:`pair\_coeff <pair_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands, or by mixing as
described below:

* epsilon (energy units)
* sigma (distance units)
* cutoff1 (distance units)
* cutoff2 (distance units)

Note that sigma is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at 2\^(1/6)
sigma.

The latter 2 coefficients are optional.  If not specified, the global
LJ and Coulombic cutoffs specified in the pair\_style command are used.
If only one cutoff is specified, it is used as the cutoff for both LJ
and Coulombic interactions for this type pair.  If both coefficients
are specified, they are used as the LJ and Coulombic cutoffs for this
type pair.  You cannot specify 2 cutoffs for style *lj/cut*\ , since it
has no Coulombic terms.

For *lj/cut/coul/long* and *lj/cut/coul/msm* and *lj/cut/tip4p/cut*
and *lj/cut/tip4p/long* only the LJ cutoff can be specified since a
Coulombic cutoff cannot be specified for an individual I,J type pair.
All type pairs use the same global Coulombic cutoff specified in the
pair\_style command.


----------


A version of these styles with a soft core, *lj/cut/soft*\ , suitable for use in
free energy calculations, is part of the USER-FEP package and is documented with
the :doc:`pair\_fep\_soft <pair_fep_soft>` styles. The version with soft core is
only available if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.


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

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed.
The default mix value is *geometric*\ .  See the "pair\_modify" command
for details.

All of the *lj/cut* pair styles support the
:doc:`pair\_modify <pair_modify>` shift option for the energy of the
Lennard-Jones portion of the pair interaction.

The *lj/cut/coul/long* and *lj/cut/tip4p/long* pair styles support the
:doc:`pair\_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the *lj/cut* pair styles support the
:doc:`pair\_modify <pair_modify>` tail option for adding a long-range
tail correction to the energy and pressure for the Lennard-Jones
portion of the pair interaction.

All of the *lj/cut* pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do
not need to be specified in an input script that reads a restart file.

The *lj/cut* and *lj/cut/coul/long* pair styles support the use of the
*inner*\ , *middle*\ , and *outer* keywords of the :doc:`run\_style respa <run_style>` command, meaning the pairwise forces can be
partitioned by distance at different levels of the rRESPA hierarchy.
The other styles only support the *pair* keyword of run\_style respa.
See the :doc:`run\_style <run_style>` command for details.


----------


Restrictions
""""""""""""


The *lj/cut/coul/long* and *lj/cut/tip4p/long* styles are part of the
KSPACE package. The *lj/cut/tip4p/cut* style is part of the MOLECULE
package. These styles are only enabled if LAMMPS was built with those
packages.  See the :doc:`Build package <Build_package>` doc page for
more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`

**Default:** none


----------


.. _Jorgensen2:



**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).

.. _Fennell2:



**(Fennell)** C. J. Fennell, J. D. Gezelter, J Chem Phys, 124,
234104 (2006).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
