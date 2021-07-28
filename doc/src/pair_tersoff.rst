.. index:: pair_style tersoff
.. index:: pair_style tersoff/table
.. index:: pair_style tersoff/gpu
.. index:: pair_style tersoff/intel
.. index:: pair_style tersoff/kk
.. index:: pair_style tersoff/omp
.. index:: pair_style tersoff/table/omp

pair_style tersoff command
==========================

Accelerator Variants: *tersoff/gpu*, *tersoff/intel*, *tersoff/kk*, *tersoff/omp*

pair_style tersoff/table command
================================

Accelerator Variants: *tersoff/table/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style keywords values

* style = *tersoff* or *tersoff/table* or *tersoff/gpu* or *tersoff/omp* or *tersoff/table/omp*
* keyword = *shift*

  .. parsed-literal::

       *shift* value = delta
         delta = negative shift in equilibrium bond length

Examples
""""""""

.. code-block:: LAMMPS

   pair_style tersoff
   pair_coeff * * Si.tersoff Si
   pair_coeff * * SiC.tersoff Si C Si

   pair_style tersoff/table
   pair_coeff * * SiCGe.tersoff Si(D)

   pair_style tersoff shift 0.05
   pair_coeff * * Si.tersoff Si

Description
"""""""""""

The *tersoff* style computes a 3-body Tersoff potential
:ref:`(Tersoff_1) <Tersoff_11>` for the energy E of a system of atoms as

.. math::

  E & = \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij} \\
  V_{ij} & = f_C(r_{ij} + \delta) \left[ f_R(r_{ij} + \delta) + b_{ij} f_A(r_{ij} + \delta) \right] \\
  f_C(r) & = \left\{ \begin{array} {r@{\quad:\quad}l}
    1 & r < R - D \\
    \frac{1}{2} - \frac{1}{2} \sin \left( \frac{\pi}{2} \frac{r-R}{D} \right) &
      R-D < r < R + D \\
    0 & r > R + D
    \end{array} \right. \\
  f_R(r) & =  A \exp (-\lambda_1 r) \\
  f_A(r) & =  -B \exp (-\lambda_2 r) \\
  b_{ij} & =  \left( 1 + \beta^n {\zeta_{ij}}^n \right)^{-\frac{1}{2n}} \\
  \zeta_{ij} & =  \sum_{k \neq i,j} f_C(r_{ik} + \delta) g \left[ \theta_{ijk}(r_{ij}, r_{ik}) \right]
                   \exp \left[ {\lambda_3}^m (r_{ij} - r_{ik})^m \right] \\
  g(\theta) & =  \gamma_{ijk} \left( 1 + \frac{c^2}{d^2} -
                  \frac{c^2}{\left[ d^2 + (\cos \theta - \cos \theta_0)^2\right]} \right)

where :math:`f_R` is a two-body term and :math:`f_A` includes three-body
interactions.  The summations in the formula are over all neighbors
J and K of atom I within a cutoff distance = R + D. :math:`\delta` is
an optional negative shift of the equilibrium bond length, as described below.

The *tersoff/table* style uses tabulated forms for the two-body,
environment and angular functions. Linear interpolation is performed
between adjacent table entries. The table length is chosen to be
accurate within 10\^-6 with respect to the *tersoff* style energy.
The *tersoff/table* should give better performance in terms of speed.

Only a single pair_coeff command is used with the *tersoff* style
which specifies a Tersoff potential file with parameters for all
needed elements.  These are mapped to LAMMPS atom types by specifying
N additional arguments after the filename in the pair_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of Tersoff elements to atom types

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the potential file.

As an example, imagine the SiC.tersoff file has Tersoff values for Si
and C.  If your LAMMPS simulation has 4 atoms types and you want the
first 3 to be Si, and the fourth to be C, you would use the following
pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * SiC.tersoff Si Si Si C

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the Tersoff file.  The final C argument maps LAMMPS atom
type 4 to the C element in the Tersoff file.  If a mapping value is
specified as NULL, the mapping is not performed.  This can be used
when a *tersoff* potential is used as part of the *hybrid* pair style.
The NULL values are placeholders for atom types that will be used with
other potentials.

Tersoff files in the *potentials* directory of the LAMMPS distribution
have a ".tersoff" suffix.  Lines that are not blank or comments
(starting with #) define parameters for a triplet of elements.  The
parameters in a single entry correspond to coefficients in the formula
above:

* element 1 (the center atom in a 3-body interaction)
* element 2 (the atom bonded to the center atom)
* element 3 (the atom influencing the 1-2 bond in a bond-order sense)
* m
* :math:`\gamma`
* :math:`\lambda_3` (1/distance units)
* c
* d
* :math:`\cos\theta_0` (can be a value < -1 or > 1)
* n
* :math:`\beta`
* :math:`\lambda_2` (1/distance units)
* B (energy units)
* R (distance units)
* D (distance units)
* :math:`\lambda_1` (1/distance units)
* A (energy units)

The n, :math:`\beta`, :math:`\lambda_2`, B, :math:`\lambda_1`, and A
parameters are only used for
two-body interactions.  The m, :math:`\gamma`, :math:`\lambda_3`, c, d,
and :math:`\cos\theta_0`
parameters are only used for three-body interactions. The R and D
parameters are used for both two-body and three-body interactions. The
non-annotated parameters are unitless.  The value of m must be 3 or 1.

The Tersoff potential file must contain entries for all the elements
listed in the pair_coeff command.  It can also contain entries for
additional elements not being used in a particular simulation; LAMMPS
ignores those entries.

For a single-element simulation, only a single entry is required
(e.g. SiSiSi).  For a two-element simulation, the file must contain 8
entries (for SiSiSi, SiSiC, SiCSi, SiCC, CSiSi, CSiC, CCSi, CCC), that
specify Tersoff parameters for all permutations of the two elements
interacting in three-body configurations.  Thus for 3 elements, 27
entries would be required, etc.

As annotated above, the first element in the entry is the center atom
in a three-body interaction and it is bonded to the second atom and the
bond is influenced by the third atom.  Thus an entry for SiCC means Si
bonded to a C with another C atom influencing the bond.  Thus
three-body parameters for SiCSi and SiSiC entries will not, in
general, be the same.  The parameters used for the two-body
interaction come from the entry where the second element is repeated.
Thus the two-body parameters for Si interacting with C, comes from the
SiCC entry.

The parameters used for a particular
three-body interaction come from the entry with the corresponding
three elements.  The parameters used only for two-body interactions
(n, :math:`\beta`, :math:`\lambda_2`, B, :math:`\lambda_1`, and A)
in entries whose second and third
element are different (e.g. SiCSi) are not used for anything and can
be set to 0.0 if desired.

Note that the twobody parameters in entries such as SiCC and CSiSi
are often the same, due to the common use of symmetric mixing rules,
but this is not always the case. For example, the beta and n parameters in
Tersoff_2 :ref:`(Tersoff_2) <Tersoff_21>` are not symmetric.
Similarly, the threebody parameters in entries such as SiCSi and SiSiC
are often the same, but this is not always the case, particularly
the value of R, which is sometimes typed on the
first and second elements, sometimes on the first and third elements.
Hence the need to specify R and D explicitly for all element triples.
For example, while Tersoff's notation
in Tersoff_2 :ref:`(Tersoff_2) <Tersoff_21>` is ambiguous on this point,
and properties of the zincblende lattice are the same for either choice,
Tersoff's results for rocksalt are consistent with typing on the first
and third elements. :ref:`Albe et al. <Albe>` adopts the same convention.
Conversely, the potential for B/N/C from the Cagin group
uses the opposite convention, typing on the first and second elements.

We chose the above form so as to enable users to define all commonly
used variants of the Tersoff potential.  In particular, our form
reduces to the original Tersoff form when m = 3 and gamma = 1, while
it reduces to the form of :ref:`Albe et al. <Albe>` when beta = 1 and m = 1.
Note that in the current Tersoff implementation in LAMMPS, m must be
specified as either 3 or 1.  Tersoff used a slightly different but
equivalent form for alloys, which we will refer to as Tersoff_2
potential :ref:`(Tersoff_2) <Tersoff_21>`.
The *tersoff/table* style implements
Tersoff_2 parameterization only.

LAMMPS parameter values for Tersoff_2 can be obtained as follows:
:math:`\gamma_{ijk} = \omega_{ik}`, :math:`\lambda_3 = 0` and the value of
m has no effect.  The parameters for species i and j can be calculated
using the Tersoff_2 mixing rules:

.. math::

   \lambda_1^{i,j} & = \frac{1}{2}(\lambda_1^i + \lambda_1^j)\\
   \lambda_2^{i,j} & = \frac{1}{2}(\lambda_2^i + \lambda_2^j)\\
   A_{i,j} & = (A_{i}A_{j})^{1/2}\\
   B_{i,j} & = \chi_{ij}(B_{i}B_{j})^{1/2}\\
   R_{i,j} & = (R_{i}R_{j})^{1/2}\\
   S_{i,j} & = (S_{i}S_{j})^{1/2}

Tersoff_2 parameters R and S must be converted to the LAMMPS
parameters R and D (R is different in both forms), using the following
relations: R=(R'+S')/2 and D=(S'-R')/2, where the primes indicate the
Tersoff_2 parameters.

In the potentials directory, the file SiCGe.tersoff provides the
LAMMPS parameters for Tersoff's various versions of Si, as well as his
alloy parameters for Si, C, and Ge. This file can be used for pure Si,
(three different versions), pure C, pure Ge, binary SiC, and binary
SiGe.  LAMMPS will generate an error if this file is used with any
combination involving C and Ge, since there are no entries for the GeC
interactions (Tersoff did not publish parameters for this
cross-interaction.)  Tersoff files are also provided for the SiC alloy
(SiC.tersoff) and the GaN (GaN.tersoff) alloys.

Many thanks to Rutuparna Narulkar, David Farrell, and Xiaowang Zhou
for helping clarify how Tersoff parameters for alloys have been
defined in various papers.

The *shift* keyword computes the energy E of a system of atoms, whose
formula is the same as the Tersoff potential. The only modification is
that the original equilibrium bond length ( :math:`r_0`) of the system
is shifted to :math:`r_0-\delta`.  The minus sign arises because each
radial distance :math:`r` is replaced by :math:`r+\delta`.

The *shift* keyword is designed for simulations of closely matched van
der Waals heterostructures.  For instance, consider the case of a system
with few-layers graphene atop a thick hexagonal boron nitride (h-BN)
substrate simulated using periodic boundary conditions. The experimental
lattice mismatch of ~1.8% between graphene and h-BN is not well captured
by the equilibrium lattice constants of available potentials, thus a
small in-plane strain will be introduced in the system when building a
periodic supercell.  To minimize the effect of strain on simulation
results, the *shift* keyword allows adjusting the equilibrium bond
length of one of the two materials (e.g., h-BN). Validation, benchmark
tests, and applications of the *shift* keyword can be found in
:ref:`(Mandelli_1) <Mandelli1>` and :ref:`(Ouyang_1) <Ouyang5>`.

For the specific case discussed above, the force field can be defined as

.. code-block:: LAMMPS

   pair_style  hybrid/overlay rebo tersoff shift -0.00407 ilp/graphene/hbn 16.0 coul/shield 16.0
   pair_coeff  * * rebo              CH.rebo     NULL NULL C
   pair_coeff  * * tersoff           BNC.tersoff B    N    NULL
   pair_coeff  * * ilp/graphene/hbn  BNCH.ILP    B    N    C
   pair_coeff  1 1 coul/shield 0.70
   pair_coeff  1 2 coul/shield 0.695
   pair_coeff  2 2 coul/shield 0.69

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as
described above from values in the potential file.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`, since it is stored in potential files.  Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the MANYBODY package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The *shift* keyword is not supported by the *tersoff/gpu*,
*tersoff/intel*, *tersoff/kk*, *tersoff/table* or *tersoff/table/omp*
variants.

The *tersoff/intel* pair style is only available when compiling LAMMPS
with the Intel compilers.

The Tersoff potential files provided with LAMMPS (see the potentials
directory) are parameterized for :doc:`"metal" units <units>`.  In addition
the pair style supports converting potential parameters on-the-fly between
"metal" and "real" units. You can use the *tersoff* pair style variants
with any LAMMPS units setting, but you would need to
create your own Tersoff potential file with coefficients listed in the
appropriate units if your simulation does not use "metal" or "real" units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

shift delta = 0.0

----------

.. _Tersoff_11:

**(Tersoff_1)** J. Tersoff, Phys Rev B, 37, 6991 (1988).

.. _Albe:

**(Albe)** J. Nord, K. Albe, P. Erhart, and K. Nordlund, J. Phys.:
Condens. Matter, 15, 5649(2003).

.. _Tersoff_21:

**(Tersoff_2)** J. Tersoff, Phys Rev B, 39, 5566 (1989); errata (PRB 41, 3248)

.. _Mandelli1:

**(Mandelli_1)** D. Mandelli, W. Ouyang, M. Urbakh, and O. Hod, ACS Nano 13(7), 7603-7609 (2019).

.. _Ouyang5:

**(Ouyang_1)** W. Ouyang et al., J. Chem. Theory Comput. 16(1), 666-676 (2020).
