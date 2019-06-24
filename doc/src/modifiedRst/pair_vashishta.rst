.. index:: pair\_style vashishta

pair\_style vashishta command
=============================

pair\_style vashishta/gpu command
=================================

pair\_style vashishta/omp command
=================================

pair\_style vashishta/kk command
================================

pair\_style vashishta/table command
===================================

pair\_style vashishta/table/omp command
=======================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *vashishta* or *vashishta/table* or *vashishta/omp* or *vashishta/table/omp*
* args = list of arguments for a particular style


.. parsed-literal::

     *vashishta* or *vashishta/omp* args = none
     *vashishta/table* or *vashishta/table/omp* args = Ntable cutinner
       Ntable = # of tabulation points
       cutinner = tablulate from cutinner to cutoff

Examples
""""""""


.. parsed-literal::

   pair_style vashishta
   pair_coeff \* \* SiC.vashishta Si C

   pair_style vashishta/table 100000 0.2
   pair_coeff \* \* SiC.vashishta Si C

Description
"""""""""""

The *vashishta* and *vashishta/table* styles compute the combined
2-body and 3-body family of potentials developed in the group of Priya
Vashishta and collaborators.  By combining repulsive, screened
Coulombic, screened charge-dipole, and dispersion interactions with a
bond-angle energy based on the Stillinger-Weber potential, this
potential has been used to describe a variety of inorganic compounds,
including SiO2 :ref:`Vashishta1990 <Vashishta1990>`, SiC
:ref:`Vashishta2007 <Vashishta2007>`, and InP :ref:`Branicio2009 <Branicio2009>`.

The potential for the energy U of a system of atoms is

.. math source doc: src/Eqs/pair_vashishta.tex
.. math::

   U & = & \sum_i^N \sum_{j > i}^N U_{ij}^{(2)} (r_{ij}) + 
   \sum_i^N \sum_{j \neq i}^N \sum_{k > j, k \neq i}^N 
   U_{ijk}^{(3)} (r_{ij}, r_{ik}, \theta_{ijk}) 
   \\
   U_{ij}^{(2)} (r) & = &  \frac{H_{ij}}{r^{\eta_{ij}}}
   	                        + \frac{Z_i Z_j}{r}\exp(-r/\lambda_{1,ij})
   	                        - \frac{D_{ij}}{r^4}\exp(-r/\lambda_{4,ij})
   	                        - \frac{W_{ij}}{r^6}, r < r_{c,{ij}}
   \\
   U_{ijk}^{(3)}(r_{ij},r_{ik},\theta_{ijk}) & = & B_{ijk} 
   \frac{\left[ \cos \theta_{ijk} - \cos \theta_{0ijk} \right]^2}
   {1+C_{ijk}\left[ \cos \theta_{ijk} - \cos \theta_{0ijk} \right]^2} \times \\
   &  &                     \exp \left( \frac{\gamma_{ij}}{r_{ij} - r_{0,ij}} \right)
   \exp \left( \frac{\gamma_{ik}}{r_{ik} - r_{0,ik}} \right), r_{ij} < r_{0,ij}, r_{ik} < r_{0,ik}


where we follow the notation used in :ref:`Branicio2009 <Branicio2009>`.
U2 is a two-body term and U3 is a three-body term.  The
summation over two-body terms is over all neighbors J within
a cutoff distance = *rc*\ .  The twobody terms are shifted and
tilted by a linear function so that the energy and force are
both zero at *rc*\ . The summation over three-body terms
is over all neighbors J and K within a cut-off distance = *r0*\ ,
where the exponential screening function becomes zero.

The *vashishta* style computes these formulas analytically.  The
*vashishta/table* style tabulates the analytic values for *Ntable*
points from cutinner to the cutoff of the potential.  The points are
equally spaced in R\^2 space from cutinner\^2 to cutoff\^2.  For the
two-body term in the above equation, a linear interpolation for each
pairwise distance between adjacent points in the table.  In practice
the tabulated version can run 3-5x faster than the analytic version
with moderate to little loss of accuracy for Ntable values
between 10000 and 1000000. It is not recommended to use less than
5000 tabulation points.

Only a single pair\_coeff command is used with either style which
specifies a Vashishta potential file with parameters for all needed
elements.  These are mapped to LAMMPS atom types by specifying N
additional arguments after the filename in the pair\_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of Vashishta elements to atom types

See the :doc:`pair\_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine a file SiC.vashishta has parameters for
Si and C.  If your LAMMPS simulation has 4 atoms types and you want
the 1st 3 to be Si, and the 4th to be C, you would use the following
pair\_coeff command:


.. parsed-literal::

   pair_coeff \* \* SiC.vashishta Si Si Si C

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the file.  The final C argument maps LAMMPS atom type 4
to the C element in the file.  If a mapping value is specified as
NULL, the mapping is not performed.  This can be used when a *vashishta*
potential is used as part of the *hybrid* pair style.  The NULL values
are placeholders for atom types that will be used with other
potentials.

Vashishta files in the *potentials* directory of the LAMMPS
distribution have a ".vashishta" suffix.  Lines that are not blank or
comments (starting with #) define parameters for a triplet of
elements.  The parameters in a single entry correspond to the two-body
and three-body coefficients in the formulae above:

* element 1 (the center atom in a 3-body interaction)
* element 2
* element 3
* H (energy units)
* eta
* Zi (electron charge units)
* Zj (electron charge units)
* lambda1 (distance units)
* D (energy units)
* lambda4 (distance units)
* W (energy units)
* rc (distance units)
* B (energy units)
* gamma
* r0 (distance units)
* C
* costheta0

The non-annotated parameters are unitless.  The Vashishta potential
file must contain entries for all the elements listed in the
pair\_coeff command.  It can also contain entries for additional
elements not being used in a particular simulation; LAMMPS ignores
those entries.  For a single-element simulation, only a single entry
is required (e.g. SiSiSi).  For a two-element simulation, the file
must contain 8 entries (for SiSiSi, SiSiC, SiCSi, SiCC, CSiSi, CSiC,
CCSi, CCC), that specify parameters for all permutations of the two
elements interacting in three-body configurations.  Thus for 3
elements, 27 entries would be required, etc.

Depending on the particular version of the Vashishta potential, the
values of these parameters may be keyed to the identities of zero,
one, two, or three elements.  In order to make the input file format
unambiguous, general, and simple to code, LAMMPS uses a slightly
confusing method for specifying parameters.  All parameters are
divided into two classes: two-body and three-body.  Two-body and
three-body parameters are handled differently, as described below.
The two-body parameters are H, eta, lambda1, D, lambda4, W, rc, gamma,
and r0.  They appear in the above formulae with two subscripts.  The
parameters Zi and Zj are also classified as two-body parameters, even
though they only have 1 subscript.  The three-body parameters are B,
C, costheta0.  They appear in the above formulae with three
subscripts.  Two-body and three-body parameters are handled
differently, as described below.

The first element in each entry is the center atom in a three-body
interaction, while the second and third elements are two neighbor
atoms. Three-body parameters for a central atom I and two neighbors J
and K are taken from the IJK entry.  Note that even though three-body
parameters do not depend on the order of J and K, LAMMPS stores
three-body parameters for both IJK and IKJ.  The user must ensure that
these values are equal.  Two-body parameters for an atom I interacting
with atom J are taken from the IJJ entry, where the 2nd and 3rd
elements are the same. Thus the two-body parameters for Si interacting
with C come from the SiCC entry. Note that even though two-body
parameters (except possibly gamma and r0 in U3) do not depend on the
order of the two elements, LAMMPS will get the Si-C value from the
SiCC entry and the C-Si value from the CSiSi entry. The user must
ensure that these values are equal. Two-body parameters appearing in
entries where the 2nd and 3rd elements are different are stored but
never used. It is good practice to enter zero for these values. Note
that the three-body function U3 above contains the two-body parameters
gamma and r0. So U3 for a central C atom bonded to an Si atom and a
second C atom will take three-body parameters from the CSiC entry, but
two-body parameters from the CCC and CSiSi entries.


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

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as
described above from values in the potential file.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


These pair style are part of the MANYBODY package.  They is only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

These pair styles requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The Vashishta potential files provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`.  You can
use the Vashishta potential with any LAMMPS units, but you would need
to create your own potential file with coefficients listed in the
appropriate units if your simulation doesn't use "metal" units.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`

**Default:** none


----------


.. _Vashishta1990:



**(Vashishta1990)** P. Vashishta, R. K. Kalia, J. P. Rino, Phys. Rev. B
41, 12197 (1990).

.. _Vashishta2007:



**(Vashishta2007)** P. Vashishta, R. K. Kalia, A. Nakano,
J. P. Rino. J. Appl. Phys. 101, 103515 (2007).

.. _Branicio2009:



**(Branicio2009)** Branicio, Rino, Gan and Tsuzuki, J. Phys Condensed
Matter 21 (2009) 095002


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
