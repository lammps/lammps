.. index:: pair\_style coul/cut

pair\_style coul/cut command
============================

pair\_style coul/cut/gpu command
================================

pair\_style coul/cut/kk command
===============================

pair\_style coul/cut/omp command
================================

pair\_style coul/debye command
==============================

pair\_style coul/debye/gpu command
==================================

pair\_style coul/debye/kk command
=================================

pair\_style coul/debye/omp command
==================================

pair\_style coul/dsf command
============================

pair\_style coul/dsf/gpu command
================================

pair\_style coul/dsf/kk command
===============================

pair\_style coul/dsf/omp command
================================

pair\_style coul/long command
=============================

pair\_style coul/long/omp command
=================================

pair\_style coul/long/gpu command
=================================

pair\_style coul/long/kk command
================================

pair\_style coul/msm command
============================

pair\_style coul/msm/omp command
================================

pair\_style coul/streitz command
================================

pair\_style coul/wolf command
=============================

pair\_style coul/wolf/kk command
================================

pair\_style coul/wolf/omp command
=================================

pair\_style tip4p/cut command
=============================

pair\_style tip4p/long command
==============================

pair\_style tip4p/cut/omp command
=================================

pair\_style tip4p/long/omp command
==================================

Syntax
""""""


.. parsed-literal::

   pair_style coul/cut cutoff
   pair_style coul/debye kappa cutoff
   pair_style coul/dsf alpha cutoff
   pair_style coul/long cutoff
   pair_style coul/long/gpu cutoff
   pair_style coul/wolf alpha cutoff
   pair_style coul/streitz cutoff keyword alpha
   pair_style tip4p/cut otype htype btype atype qdist cutoff
   pair_style tip4p/long otype htype btype atype qdist cutoff

* cutoff = global cutoff for Coulombic interactions
* kappa = Debye length (inverse distance units)
* alpha = damping parameter (inverse distance units)

Examples
""""""""


.. parsed-literal::

   pair_style coul/cut 2.5
   pair_coeff \* \*
   pair_coeff 2 2 3.5

   pair_style coul/debye 1.4 3.0
   pair_coeff \* \*
   pair_coeff 2 2 3.5

   pair_style coul/dsf 0.05 10.0
   pair_coeff \* \*

   pair_style coul/long 10.0
   pair_coeff \* \*

   pair_style coul/msm 10.0
   pair_coeff \* \*

   pair_style coul/wolf 0.2 9.0
   pair_coeff \* \*

   pair_style coul/streitz 12.0 ewald
   pair_style coul/streitz 12.0 wolf 0.30
   pair_coeff \* \* AlO.streitz Al O

   pair_style tip4p/cut 1 2 7 8 0.15 12.0
   pair_coeff \* \*

   pair_style tip4p/long 1 2 7 8 0.15 10.0
   pair_coeff \* \*

Description
"""""""""""

The *coul/cut* style computes the standard Coulombic interaction
potential given by

.. math::

  E = \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c


where C is an energy-conversion constant, Qi and Qj are the charges on
the 2 atoms, and :math:`\epsilon` is the dielectric constant which can be set
by the :doc:`dielectric <dielectric>` command.  The cutoff :math:`r_c` truncates
the interaction distance.


----------


Style *coul/debye* adds an additional exp() damping factor to the
Coulombic term, given by

.. math::

  E = \frac{C q_i q_j}{\epsilon  r} \exp(- \kappa r) \qquad r < r_c


where :math:`\kappa` is the Debye length.  This potential is another way to
mimic the screening effect of a polar solvent.


----------


Style *coul/dsf* computes Coulombic interactions via the damped
shifted force model described in :ref:`Fennell <Fennell1>`, given by:

.. math::

  E = 
  q_iq_j \left[ \frac{\mbox{erfc} (\alpha r)}{r} -  \frac{\mbox{erfc} (\alpha r_c)}{r_c} + 
  \left( \frac{\mbox{erfc} (\alpha r_c)}{r_c^2} +  \frac{2\alpha}{\sqrt{\pi}}\frac{\exp (-\alpha^2    r^2_c)}{r_c} \right)(r-r_c) \right] \qquad r < r_c 


where :math:`\alpha` is the damping parameter and erfc() is the
complementary error-function. The potential corrects issues in the
Wolf model (described below) to provide consistent forces and energies
(the Wolf potential is not differentiable at the cutoff) and smooth
decay to zero.


----------


Style *coul/wolf* computes Coulombic interactions via the Wolf
summation method, described in :ref:`Wolf <Wolf1>`, given by:

.. math::

  E_i = \frac{1}{2} \sum_{j \neq i} 
  \frac{q_i q_j {\rm erfc}(\alpha r_{ij})}{r_{ij}} + 
  \frac{1}{2} \sum_{j \neq i} 
  \frac{q_i q_j {\rm erf}(\alpha r_{ij})}{r_{ij}} \qquad r < r_c


where :math:`\alpha` is the damping parameter, and erc() and erfc() are
error-function and complementary error-function terms.  This potential
is essentially a short-range, spherically-truncated,
charge-neutralized, shifted, pairwise *1/r* summation.  With a
manipulation of adding and subtracting a self term (for i = j) to the
first and second term on the right-hand-side, respectively, and a
small enough :math:`\alpha` damping parameter, the second term shrinks and
the potential becomes a rapidly-converging real-space summation.  With
a long enough cutoff and small enough alpha parameter, the energy and
forces calculated by the Wolf summation method approach those of the
Ewald sum.  So it is a means of getting effective long-range
interactions with a short-range potential.


----------


Style *coul/streitz* is the Coulomb pair interaction defined as part
of the Streitz-Mintmire potential, as described in :ref:`this paper <Streitz2>`, in which charge distribution about an atom is modeled
as a Slater 1\ *s* orbital.  More details can be found in the referenced
paper.  To fully reproduce the published Streitz-Mintmire potential,
which is a variable charge potential, style *coul/streitz* must be
used with :doc:`pair_style eam/alloy <pair_eam>` (or some other
short-range potential that has been parameterized appropriately) via
the :doc:`pair_style hybrid/overlay <pair_hybrid>` command.  Likewise,
charge equilibration must be performed via the :doc:`fix qeq/slater <fix_qeq>` command. For example:


.. parsed-literal::

   pair_style hybrid/overlay coul/streitz 12.0 wolf 0.31 eam/alloy
   pair_coeff \* \* coul/streitz AlO.streitz Al O
   pair_coeff \* \* eam/alloy AlO.eam.alloy Al O
   fix 1 all qeq/slater 1 12.0 1.0e-6 100 coul/streitz

The keyword *wolf* in the coul/streitz command denotes computing
Coulombic interactions via Wolf summation.  An additional damping
parameter is required for the Wolf summation, as described for the
coul/wolf potential above.  Alternatively, Coulombic interactions can
be computed via an Ewald summation.  For example:


.. parsed-literal::

   pair_style hybrid/overlay coul/streitz 12.0 ewald eam/alloy
   kspace_style ewald 1e-6

Keyword *ewald* does not need a damping parameter, but a
:doc:`kspace_style <kspace_style>` must be defined, which can be style
*ewald* or *pppm*\ .  The Ewald method was used in Streitz and
Mintmire's original paper, but a Wolf summation offers a speed-up in
some cases.

For the fix qeq/slater command, the *qfile* can be a filename that
contains QEq parameters as discussed on the :doc:`fix qeq <fix_qeq>`
command doc page.  Alternatively *qfile* can be replaced by
"coul/streitz", in which case the fix will extract QEq parameters from
the coul/streitz pair style itself.

See the examples/strietz directory for an example input script that
uses the Streitz-Mintmire potential.  The potentials directory has the
AlO.eam.alloy and AlO.streitz potential files used by the example.

Note that the Streiz-Mintmire potential is generally used for oxides,
but there is no conceptual problem with extending it to nitrides and
carbides (such as SiC, TiN).  Pair coul/strietz used by itself or with
any other pair style such as EAM, MEAM, Tersoff, or LJ in
hybrid/overlay mode.  To do this, you would need to provide a
Streitz-Mintmire parameterization for the material being modeled.


----------


Styles *coul/long* and *coul/msm* compute the same Coulombic
interactions as style *coul/cut* except that an additional damping
factor is applied so it can be used in conjunction with the
:doc:`kspace_style <kspace_style>` command and its *ewald* or *pppm*
option.  The Coulombic cutoff specified for this style means that
pairwise interactions within this distance are computed directly;
interactions outside that distance are computed in reciprocal space.

Styles *tip4p/cut* and *tip4p/long* implement the Coulomb part of
the TIP4P water model of :ref:`(Jorgensen) <Jorgensen3>`, which introduces
a massless site located a short distance away from the oxygen atom
along the bisector of the HOH angle.  The atomic types of the oxygen and
hydrogen atoms, the bond and angle types for OH and HOH interactions,
and the distance to the massless charge site are specified as
pair\_style arguments.  Style *tip4p/cut* uses a global cutoff for
Coulomb interactions; style *tip4p/long* is for use with a long-range
Coulombic solver (Ewald or PPPM).

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


----------


Note that these potentials are designed to be combined with other pair
potentials via the :doc:`pair_style hybrid/overlay <pair_hybrid>`
command.  This is because they have no repulsive core.  Hence if they
are used by themselves, there will be no repulsion to keep two
oppositely charged particles from moving arbitrarily close to each
other.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* cutoff (distance units)

For *coul/cut* and *coul/debye*\ , the cutoff coefficient is optional.
If it is not used (as in some of the examples above), the default
global value specified in the pair\_style command is used.

For *coul/long* and *coul/msm* no cutoff can be specified for an
individual I,J type pair via the pair\_coeff command.  All type pairs
use the same global Coulomb cutoff specified in the pair\_style
command.


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

For atom type pairs I,J and I != J, the cutoff distance for the
*coul/cut* style can be mixed.  The default mix value is *geometric*\ .
See the "pair\_modify" command for details.

The :doc:`pair_modify <pair_modify>` shift option is not relevant
for these pair styles.

The *coul/long* style supports the :doc:`pair_modify <pair_modify>`
table option for tabulation of the short-range portion of the
long-range Coulombic interaction.

These pair styles do not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

These pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


The *coul/long*\ , *coul/msm* and *tip4p/long* styles are part of the
KSPACE package.  They are only enabled if LAMMPS was built with that
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style, hybrid/overlay <pair_hybrid>`, :doc:`kspace_style <kspace_style>`

**Default:** none


----------


.. _Wolf1:



**(Wolf)** D. Wolf, P. Keblinski, S. R. Phillpot, J. Eggebrecht, J Chem
Phys, 110, 8254 (1999).

.. _Fennell1:



**(Fennell)** C. J. Fennell, J. D. Gezelter, J Chem Phys, 124,
234104 (2006).

.. _Streitz2:



**(Streitz)** F. H. Streitz, J. W. Mintmire, Phys Rev B, 50, 11996-12003
(1994).

.. _Jorgensen3:



**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).
