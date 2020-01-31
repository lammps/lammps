.. index:: pair\_style tersoff/zbl

pair\_style tersoff/zbl command
===============================

pair\_style tersoff/zbl/gpu command
===================================

pair\_style tersoff/zbl/kk command
==================================

pair\_style tersoff/zbl/omp command
===================================

Syntax
""""""


.. parsed-literal::

   pair_style tersoff/zbl

Examples
""""""""


.. parsed-literal::

   pair_style tersoff/zbl
   pair_coeff \* \* SiC.tersoff.zbl Si C Si

Description
"""""""""""

The *tersoff/zbl* style computes a 3-body Tersoff potential
:ref:`(Tersoff\_1) <zbl-Tersoff_1>` with a close-separation pairwise modification
based on a Coulomb potential and the Ziegler-Biersack-Littmark
universal screening function :ref:`(ZBL) <zbl-ZBL>`, giving the energy E of a
system of atoms as

.. image:: Eqs/pair_tersoff_zbl.jpg
   :align: center

The f\_F term is a fermi-like function used to smoothly connect the ZBL
repulsive potential with the Tersoff potential.  There are 2
parameters used to adjust it: A\_F and r\_C.  A\_F controls how "sharp"
the transition is between the two, and r\_C is essentially the cutoff
for the ZBL potential.

For the ZBL portion, there are two terms. The first is the Coulomb
repulsive term, with Z1, Z2 as the number of protons in each nucleus,
e as the electron charge (1 for metal and real units) and epsilon0 as
the permittivity of vacuum.  The second part is the ZBL universal
screening function, with a0 being the Bohr radius (typically 0.529
Angstroms), and the remainder of the coefficients provided by the
original paper.  This screening function should be applicable to most
systems.  However, it is only accurate for small separations
(i.e. less than 1 Angstrom).

For the Tersoff portion, f\_R is a two-body term and f\_A includes
three-body interactions. The summations in the formula are over all
neighbors J and K of atom I within a cutoff distance = R + D.

Only a single pair\_coeff command is used with the *tersoff/zbl* style
which specifies a Tersoff/ZBL potential file with parameters for all
needed elements.  These are mapped to LAMMPS atom types by specifying
N additional arguments after the filename in the pair\_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of Tersoff/ZBL elements to atom types

See the :doc:`pair_coeff <pair_coeff>` doc page for alternate ways
to specify the path for the potential file.

As an example, imagine the SiC.tersoff.zbl file has Tersoff/ZBL values
for Si and C.  If your LAMMPS simulation has 4 atoms types and you
want the 1st 3 to be Si, and the 4th to be C, you would use the
following pair\_coeff command:


.. parsed-literal::

   pair_coeff \* \* SiC.tersoff Si Si Si C

The 1st 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the Tersoff/ZBL file.  The final C argument maps LAMMPS
atom type 4 to the C element in the Tersoff/ZBL file.  If a mapping
value is specified as NULL, the mapping is not performed.  This can be
used when a *tersoff/zbl* potential is used as part of the *hybrid*
pair style.  The NULL values are placeholders for atom types that will
be used with other potentials.

Tersoff/ZBL files in the *potentials* directory of the LAMMPS
distribution have a ".tersoff.zbl" suffix.  Lines that are not blank
or comments (starting with #) define parameters for a triplet of
elements.  The parameters in a single entry correspond to coefficients
in the formula above:

* element 1 (the center atom in a 3-body interaction)
* element 2 (the atom bonded to the center atom)
* element 3 (the atom influencing the 1-2 bond in a bond-order sense)
* m
* gamma
* lambda3 (1/distance units)
* c
* d
* costheta0 (can be a value < -1 or > 1)
* n
* beta
* lambda2 (1/distance units)
* B (energy units)
* R (distance units)
* D (distance units)
* lambda1 (1/distance units)
* A (energy units)
* Z\_i
* Z\_j
* ZBLcut (distance units)
* ZBLexpscale (1/distance units)

The n, beta, lambda2, B, lambda1, and A parameters are only used for
two-body interactions.  The m, gamma, lambda3, c, d, and costheta0
parameters are only used for three-body interactions. The R and D
parameters are used for both two-body and three-body interactions. The
Z\_i,Z\_j, ZBLcut, ZBLexpscale parameters are used in the ZBL repulsive
portion of the potential and in the Fermi-like function.  The
non-annotated parameters are unitless.  The value of m must be 3 or 1.

The Tersoff/ZBL potential file must contain entries for all the
elements listed in the pair\_coeff command.  It can also contain
entries for additional elements not being used in a particular
simulation; LAMMPS ignores those entries.

For a single-element simulation, only a single entry is required
(e.g. SiSiSi).  For a two-element simulation, the file must contain 8
entries (for SiSiSi, SiSiC, SiCSi, SiCC, CSiSi, CSiC, CCSi, CCC), that
specify Tersoff parameters for all permutations of the two elements
interacting in three-body configurations.  Thus for 3 elements, 27
entries would be required, etc.

As annotated above, the first element in the entry is the center atom
in a three-body interaction and it is bonded to the 2nd atom and the
bond is influenced by the 3rd atom.  Thus an entry for SiCC means Si
bonded to a C with another C atom influencing the bond.  Thus
three-body parameters for SiCSi and SiSiC entries will not, in
general, be the same.  The parameters used for the two-body
interaction come from the entry where the 2nd element is repeated.
Thus the two-body parameters for Si interacting with C, comes from the
SiCC entry.

The parameters used for a particular
three-body interaction come from the entry with the corresponding
three elements.  The parameters used only for two-body interactions
(n, beta, lambda2, B, lambda1, and A) in entries whose 2nd and 3rd
element are different (e.g. SiCSi) are not used for anything and can
be set to 0.0 if desired.

Note that the twobody parameters in entries such as SiCC and CSiSi
are often the same, due to the common use of symmetric mixing rules,
but this is not always the case. For example, the beta and n parameters in
Tersoff\_2 :ref:`(Tersoff\_2) <zbl-Tersoff_2>` are not symmetric.

We chose the above form so as to enable users to define all commonly
used variants of the Tersoff portion of the potential.  In particular,
our form reduces to the original Tersoff form when m = 3 and gamma =
1, while it reduces to the form of :ref:`Albe et al. <zbl-Albe>` when beta = 1
and m = 1.  Note that in the current Tersoff implementation in LAMMPS,
m must be specified as either 3 or 1.  Tersoff used a slightly
different but equivalent form for alloys, which we will refer to as
Tersoff\_2 potential :ref:`(Tersoff\_2) <zbl-Tersoff_2>`.

LAMMPS parameter values for Tersoff\_2 can be obtained as follows:
gamma = omega\_ijk, lambda3 = 0 and the value of
m has no effect.  The parameters for species i and j can be calculated
using the Tersoff\_2 mixing rules:

.. image:: Eqs/pair_tersoff_2.jpg
   :align: center

Tersoff\_2 parameters R and S must be converted to the LAMMPS
parameters R and D (R is different in both forms), using the following
relations: R=(R'+S')/2 and D=(S'-R')/2, where the primes indicate the
Tersoff\_2 parameters.

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
defined in various papers.  Also thanks to Ram Devanathan for
providing the base ZBL implementation.


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

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair\_style and pair\_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


This pair style is part of the MANYBODY package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The Tersoff/ZBL potential files provided with LAMMPS (see the
potentials directory) are parameterized for metal :doc:`units <units>`.
You can use the Tersoff potential with any LAMMPS units, but you would
need to create your own Tersoff potential file with coefficients
listed in the appropriate units if your simulation doesn't use "metal"
units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

**Default:** none


----------


.. _zbl-Tersoff\_1:



**(Tersoff\_1)** J. Tersoff, Phys Rev B, 37, 6991 (1988).

.. _zbl-ZBL:



**(ZBL)** J.F. Ziegler, J.P. Biersack, U. Littmark, 'Stopping and Ranges
of Ions in Matter' Vol 1, 1985, Pergamon Press.

.. _zbl-Albe:



**(Albe)** J. Nord, K. Albe, P. Erhart and K. Nordlund, J. Phys.:
Condens. Matter, 15, 5649(2003).

.. _zbl-Tersoff\_2:



**(Tersoff\_2)** J. Tersoff, Phys Rev B, 39, 5566 (1989); errata (PRB 41, 3248)


