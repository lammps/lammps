.. index:: pair_style threebody/table

pair_style threebody/table command
==================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *threebody/table*


Examples
""""""""

.. code-block:: LAMMPS

   pair_style threebody/table
   pair_coeff * * spce2.3b type1 type2

   pair_style hybrid/overlay table linear 1200 threebody/table
   pair_coeff      1 1 table table_CG_CG.txt VOTCA
   pair_coeff      * * threebody/table spce.3b type

Used in example input scripts:

.. parsed-literal::

   examples/PACKAGES/manybody_table/in.spce
   examples/PACKAGES/manybody_table/in.spce2

Description
"""""""""""

.. versionadded:: 2Jun2022

The *threebody/table* style is a pair style for generic tabulated
three-body interactions.  It has been developed for (coarse-grained)
simulations (of water) with Kernel-based machine learning (ML)
potentials (:ref:`Scherer2 <Scherer2>`).  As for many other MANYBODY
package pair styles the energy of a system is computed as a sum over
three-body terms:

.. math::

   E =  \sum_i \sum_{j \neq i} \sum_{k > j} \phi_3 (r_{ij}, r_{ik}, \theta_{ijk})

The summations in the formula are over all neighbors J and K of atom I
within a cutoff distance :math:`cut`.  In contrast to the
Stillinger-Weber potential, all forces are not calculated analytically,
but read in from a three-body force/energy table which can be generated
with the csg_ml app of VOTCA as available at:
https://gitlab.mpcdf.mpg.de/votca/votca.

Only a single pair_coeff command is used with the *threebody/table*
style which specifies a threebody potential (".3b") file with parameters
for all needed elements. These are then mapped to LAMMPS atom types by
specifying N_el additional arguments after the ".3b" filename in the
pair_coeff command, where N_el is the number of LAMMPS atom types:

* ".3b" filename
* N_el element names = mapping of threebody elements to atom types

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways to
specify the path for the potential file.

As an example, imagine a file SiC.3b has three-body values for Si and C.
If your LAMMPS simulation has 4 atoms types and you want the first 3 to
be Si, and the fourth to be C, you would use the following pair_coeff
command:

.. code-block:: LAMMPS

   pair_coeff * * SiC.3b Si Si Si C

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the ".3b" file.  The final C argument maps LAMMPS atom type 4
to the C element in the threebody file. If a mapping value is specified
as NULL, the mapping is not performed. This can be used when a
*threebody/table* potential is used as part of the *hybrid* pair style.
The NULL values are placeholders for atom types that will be used with
other potentials.

The three-body files have a ".3b" suffix. Lines that are not blank or
comments (starting with #) define parameters for a triplet of
elements. The parameters in a single entry specify to the (three-body)
cutoff distance and the tabulated three-body interaction. A single entry
then contains:

* element 1 (the center atom in a 3-body interaction)
* element 2
* element 3
* cut (distance units)
* filename
* keyword
* style
* N

The parameter :math:`cut` is the (three-body) cutoff distance.  When set
to 0, no interaction is calculated for this element triplet.  The
parameters *filename*, *keyword*, *style*, and *N* refer to the
tabulated three-body potential.

The tabulation is done on a three-dimensional grid of the two distances
:math:`r_{ij}` and :math:`r_{ik}` as well as the angle
:math:`\theta_{ijk}` which is constructed in the following way. There
are two different cases.  If element 2 and element 3 are of the same
type (e.g. SiCC), the distance :math:`r_{ij}` is varied in "N" steps
from rmin to rmax and the distance :math:`r_{ik}` is varied from
:math:`r_{ij}` to rmax. This can be done, due to the symmetry of the
triplet. If element 2 and element 3 are not of the same type
(e.g. SiCSi), there is no additional symmetry and the distance
:math:`r_{ik}` is also varied from rmin to rmax in "N" steps.  The angle
:math:`\theta_{ijk}` is always varied in "2N" steps from (0.0 +
180.0/(4N)) to (180.0 - 180.0/(4N)). Therefore, the total number of
table entries is "M = N * N * (N+1)" for the symmetric (element 2 and
element 3 are of the same type) and "M = 2 * N * N * N" for the general
case (element 2 and element 3 are not of the same type).

The forces on all three particles I, J, and K of a triplet of this type
of three-body interaction potential (:math:`\phi_3 (r_{ij}, r_{ik},
\theta_{ijk})`) lie within the plane defined by the three inter-particle
distance vectors :math:`{\mathbf r}_{ij}`, :math:`{\mathbf r}_{ik}`, and
:math:`{\mathbf r}_{jk}`.  This property is used to project the forces
onto the inter-particle distance vectors as follows

.. math::

   \begin{pmatrix}
      {\mathbf f}_{i} \\
      {\mathbf f}_{j} \\
      {\mathbf f}_{k} \\
   \end{pmatrix} =
   \begin{pmatrix}
      f_{i1} & f_{i2} & 0 \\
      f_{j1} & 0 & f_{j2} \\
      0 & f_{k1} & f_{k2} \\
   \end{pmatrix}
   \begin{pmatrix}
      {\mathbf r}_{ij} \\
      {\mathbf r}_{ik} \\
      {\mathbf r}_{jk} \\
   \end{pmatrix}

and then tabulate the 6 force constants :math:`f_{i1}`, :math:`f_{i2}`,
:math:`f_{j1}`, :math:`f_{j2}`, :math:`f_{k1}`, and :math:`f_{k2}`, as
well as the energy of a triplet e.  Due to symmetry reasons, the
following relations hold: :math:`f_{i1}=-f_{j1}`,
:math:`f_{i2}=-f_{k1}`, and :math:`f_{j2}=-f_{k2}`. As in this pair
style the forces are read in directly, a correct MD simulation is also
performed in the case that the triplet energies are set to e=0.

The *filename* specifies the file containing the tabulated energy and
derivative values of :math:`\phi_3 (r_{ij}, r_{ik}, \theta_{ijk})`.  The
*keyword* then specifies a section of the file. The format of this file
is as follows (without the parenthesized comments):

.. parsed-literal::

   # Tabulated three-body potential for spce water (one or more comment or blank lines)

   ENTRY1                                                                      (keyword is the first text on line)
   N 12 rmin 2.55 rmax 3.65                                                    (N, rmin, rmax parameters)
                                                                               (blank line)
   1 2.55 2.55 3.75 -867.212 -611.273 867.212 21386.8 611.273 -21386.8  0.0    (index, r_ij, r_ik, theta, f_i1, f_i2, f_j1, f_j2, f_k1, f_k2, e)
   2 2.55 2.55 11.25 -621.539 -411.189 621.539 5035.95 411.189 -5035.95  0.0
   ...
   1872 3.65 3.65 176.25 -0.00215132 -0.00412886 0.00215137 0.00111754 0.00412895 -0.00111757  0.0

A section begins with a non-blank line whose first character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which identifies
the section. The next line lists (in any order) one or more parameters
for the table.  Each parameter is a keyword followed by one or more
numeric values.

The parameter "N" is required. It should be the same than the parameter
"N" of the ".3b" file, otherwise its value is overwritten. "N"
determines the number of table entries "M" that follow: "M = N * N *
(N+1)" (symmetric triplet, e.g. SiCC) or "M = 2 * N * N * N" (asymmetric
triplet, e.g. SiCSi). Therefore "M = 12 * 12 * 13 = 1872" in the above
symmetric example. The parameters "rmin" and "rmax" are also required
and determine the minimum and maximum of the inter-particle distances
:math:`r_{ij}` and :math:`r_{ik}`.

Following a blank line, the next M lines of the angular table file list
the tabulated values.  On each line, the first value is the index from 1
to M, the second value is the distance :math:`r_{ij}`, the third value
is the distance :math:`r_{ik}`, the fourth value is the angle
:math:`\theta_{ijk})`, the next six values are the force constants
:math:`f_{i1}`, :math:`f_{i2}`, :math:`f_{j1}`, :math:`f_{j2}`,
:math:`f_{k1}`, and :math:`f_{k2}`, and the last value is the energy e.

Note that one three-body potential file can contain many sections, each
with a tabulated potential. LAMMPS reads the file section by section
until it finds one that matches the specified *keyword* of appropriate
section of the ".3b" file.

At the moment, only the *style* *linear* is allowed and
implemented. After reading in the force table, it is internally stored
in LAMMPS as a lookup table. For each triplet configuration occurring in
the simulation within the cutoff distance, the next nearest tabulated
triplet configuration is looked up. No interpolation is done.  This
allows for a very efficient force calculation with the stored force
constants and energies. Due to the know table structure, the lookup can
be done efficiently. It has been tested (:ref:`Scherer2 <Scherer2>`)
that with a reasonably small bin size, the accuracy and speed is
comparable to that of a Stillinger-Weber potential with tabulated
three-body interactions (:doc:`pair_style sw/angle/table
<pair_sw_angle_table>`) while the table format of this pair style allows
for more flexible three-body interactions.

As for the Stillinger-Weber potential, the three-body potential file
must contain entries for all the elements listed in the pair_coeff
command.  It can also contain entries for additional elements not being
used in a particular simulation; LAMMPS ignores those entries.

For a single-element simulation, only a single entry is required
(e.g. SiSiSi).  For a two-element simulation, the file must contain 8
entries (for SiSiSi, SiSiC, SiCSi, SiCC, CSiSi, CSiC, CCSi, CCC), that
specify threebody parameters for all permutations of the two elements
interacting in three-body configurations.  Thus for 3 elements, 27
entries would be required, etc.

As annotated above, the first element in the entry is the center atom in
a three-body interaction. Thus an entry for SiCC means a Si atom with 2
C atoms as neighbors. The tabulated three-body forces can in principle
be specific to the three elements of the configuration.  However, the
user must ensure that it makes physically sense.  E.g., the tabulated
three-body forces for the entries CSiC and CCSi should be the same
exchanging :math:`r_{ij}` with r_{ik}, :math:`f_{j1}` with
:math:`f_{k1}`, and :math:`f_{j2}` with :math:`f_{k2}`.

Additional input files and reference data can be found at:
https://gitlab.mpcdf.mpg.de/votca/votca/-/tree/master/csg-tutorials/ml

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As all interactions are tabulated, no mixing is performed.

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
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair sw/angle/table <pair_sw_angle_table>`


----------

.. _Scherer2:

**(Scherer2)** C. Scherer, R. Scheid, D. Andrienko, and T. Bereau, J. Chem. Theor. Comp. 16, 3194-3204 (2020).

