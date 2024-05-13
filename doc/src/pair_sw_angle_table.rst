.. index:: pair_style sw/angle/table

pair_style sw/angle/table command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *sw/angle/table*


Examples
""""""""

.. code-block:: LAMMPS

   pair_style sw/angle/table
   pair_coeff * * spce.sw type

Used in example input script:

.. parsed-literal::

   examples/PACKAGES/manybody_table/in.spce_sw


Description
"""""""""""

.. versionadded:: 2Jun2022

The *sw/angle/table* style is a modification of the original
:doc:`pair_style sw <pair_sw>`. It has been developed for coarse-grained
simulations (of water) (:ref:`Scherer1 <Scherer1>`), but can be employed
for all kinds of systems. It computes a modified 3-body
:ref:`Stillinger-Weber <Stillinger3>` potential for the energy E of a
system of atoms as

.. math::

   E & =  \sum_i \sum_{j > i} \phi_2 (r_{ij}) +
          \sum_i \sum_{j \neq i} \sum_{k > j}
          \phi_3 (r_{ij}, r_{ik}, \theta_{ijk}) \\
  \phi_2(r_{ij}) & =  A_{ij} \epsilon_{ij} \left[ B_{ij} (\frac{\sigma_{ij}}{r_{ij}})^{p_{ij}} -
                    (\frac{\sigma_{ij}}{r_{ij}})^{q_{ij}} \right]
                    \exp \left( \frac{\sigma_{ij}}{r_{ij} - a_{ij} \sigma_{ij}} \right) \\
  \phi_3(r_{ij},r_{ik},\theta_{ijk}) & = f^{\textrm{3b}}\left(\theta_{ijk}\right)
                    \exp \left( \frac{\gamma_{ij} \sigma_{ij}}{r_{ij} - a_{ij} \sigma_{ij}} \right)
                    \exp \left( \frac{\gamma_{ik} \sigma_{ik}}{r_{ik} - a_{ik} \sigma_{ik}} \right)

where :math:`\phi_2` is a two-body term and :math:`\phi_3` is a
three-body term.  The summations in the formula are over all neighbors J
and K of atom I within a cutoff distance :math:`a \sigma`.  In contrast
to the original *sw* style, *sw/angle/table* allows for a flexible
three-body term :math:`f^{\textrm{3b}}\left(\theta_{ijk}\right)` which
is read in as a tabulated interaction. It can be parameterized with the
csg_fmatch app of VOTCA as available at:
https://gitlab.mpcdf.mpg.de/votca/votca.

Only a single pair_coeff command is used with the *sw/angle/table* style
which specifies a modified Stillinger-Weber potential file with
parameters for all needed elements.  These are mapped to LAMMPS atom
types by specifying N_el additional arguments after the ".sw" filename
in the pair_coeff command, where N_el is the number of LAMMPS atom
types:

* ".sw" filename
* N_el element names = mapping of SW elements to atom types

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways to
specify the path for the potential file.

As an example, imagine a file SiC.sw has Stillinger-Weber values for Si
and C.  If your LAMMPS simulation has 4 atoms types and you want the
first 3 to be Si, and the fourth to be C, you would use the following
pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * SiC.sw Si Si Si C

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The first three Si arguments map LAMMPS atom types 1,2,3 to the Si
element in the SW file.  The final C argument maps LAMMPS atom type 4 to
the C element in the SW file.  If a mapping value is specified as NULL,
the mapping is not performed.  This can be used when a *sw/angle/table*
potential is used as part of the *hybrid* pair style.  The NULL values
are placeholders for atom types that will be used with other potentials.

The (modified) Stillinger-Weber files have a ".sw" suffix. Lines that
are not blank or comments (starting with #) define parameters for a
triplet of elements. The parameters in a single entry correspond to the
two-body and three-body coefficients in the formula above. Here, also
the suffix ".sw" is used though the original Stillinger-Weber file
format is supplemented with four additional lines per parameter block to
specify the tabulated three-body interaction. A single entry then
contains:

* element 1 (the center atom in a 3-body interaction)
* element 2
* element 3
* :math:`\epsilon` (energy units)
* :math:`\sigma` (distance units)
* a
* :math:`\lambda`
* :math:`\gamma`
* :math:`\cos\theta_0`
* A
* B
* p
* q
* tol
* filename
* keyword
* style
* N

The A, B, p, and q parameters are used only for two-body interactions.
The :math:`\lambda` and :math:`\cos\theta_0` parameters, only used for
three-body interactions in the original Stillinger-Weber style, are read
in but ignored in this modified pair style. The :math:`\epsilon`
parameter is only used for two-body interactions in this modified pair
style and not for the three-body terms. The :math:`\sigma` and *a*
parameters are used for both two-body and three-body
interactions. :math:`\gamma` is used only in the three-body
interactions, but is defined for pairs of atoms. The non-annotated
parameters are unitless.

LAMMPS introduces an additional performance-optimization parameter tol
that is used for both two-body and three-body interactions.  In the
Stillinger-Weber potential, the interaction energies become negligibly
small at atomic separations substantially less than the theoretical
cutoff distances.  LAMMPS therefore defines a virtual cutoff distance
based on a user defined tolerance tol.  The use of the virtual cutoff
distance in constructing atom neighbor lists can significantly reduce
the neighbor list sizes and therefore the computational cost.  LAMMPS
provides a *tol* value for each of the three-body entries so that they
can be separately controlled. If tol = 0.0, then the standard
Stillinger-Weber cutoff is used.

The additional parameters *filename*, *keyword*, *style*, and *N* refer
to the tabulated angular potential
:math:`f^{\textrm{3b}}\left(\theta_{ijk}\right)`.  The tabulated angular
potential has to be of the format as used in the :doc:`angle_style table
<angle_table>` command:

An interpolation tables of length *N* is created. The interpolation is
done in one of 2 *styles*: *linear* or *spline*.  For the *linear*
style, the angle is used to find 2 surrounding table values from which
an energy or its derivative is computed by linear interpolation. For the
*spline* style, a cubic spline coefficients are computed and stored at
each of the *N* values in the table.  The angle is used to find the
appropriate set of coefficients which are used to evaluate a cubic
polynomial which computes the energy or derivative.

The *filename* specifies the file containing the tabulated energy and
derivative values of :math:`f^{\textrm{3b}}\left(\theta_{ijk}\right)`.
The *keyword* then specifies a section of the file.  The format of this
file is as follows (without the parenthesized comments):

.. parsed-literal::

   # Angle potential for harmonic (one or more comment or blank lines)

   HAM                           (keyword is the first text on line)
   N 181 FP 0 0 EQ 90.0          (N, FP, EQ parameters)
                                 (blank line)
   1 0.0 200.5 2.5               (index, angle, energy, derivative)
   2 1.0 198.0 2.5
   ...
   181 180.0 0.0 0.0

A section begins with a non-blank line whose first character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which identifies
the section. The next line lists (in any order) one or more parameters
for the table.  Each parameter is a keyword followed by one or more
numeric values.

The parameter "N" is required and its value is the number of table
entries that follow.  Note that this may be different than the *N*
specified in the Stillinger-Weber potential file. Let Nsw = *N* in the
".sw" file, and Nfile = "N" in the tabulated angular file.  What LAMMPS
does is a preliminary interpolation by creating splines using the Nfile
tabulated values as nodal points.  It uses these to interpolate as
needed to generate energy and derivative values at Ntable different
points.  The resulting tables of length Nsw are then used as described
above, when computing energy and force for individual angles and their
atoms.  This means that if you want the interpolation tables of length
Nsw to match exactly what is in the tabulated file (with effectively no
preliminary interpolation), you should set Nsw = Nfile.

The "FP" parameter is optional.  If used, it is followed by two values
fplo and fphi, which are the second derivatives at the innermost and
outermost angle settings.  These values are needed by the spline
construction routines.  If not specified by the "FP" parameter, they are
estimated (less accurately) by the first two and last two derivative
values in the table.

The "EQ" parameter is also optional.  If used, it is followed by a the
equilibrium angle value, which is used, for example, by the :doc:`fix
shake <fix_shake>` command. If not used, the equilibrium angle is set to
180.0.

Following a blank line, the next N lines of the angular table file list
the tabulated values.  On each line, the first value is the index from 1
to N, the second value is the angle value (in degrees), the third value
is the energy (in energy units), and the fourth is -dE/d(theta) (also in
energy units).  The third term is the energy of the 3-atom configuration
for the specified angle.  The last term is the derivative of the energy
with respect to the angle (in degrees, not radians).  Thus the units of
the last term are still energy, not force.  The angle values must
increase from one line to the next.  The angle values must also begin
with 0.0 and end with 180.0, i.e. span the full range of possible
angles.

Note that one angular potential file can contain many sections, each
with a tabulated potential.  LAMMPS reads the file section by section
until it finds one that matches the specified *keyword* of appropriate
section of the ".sw" file.

The Stillinger-Weber potential file must contain entries for all the
elements listed in the pair_coeff command.  It can also contain entries
for additional elements not being used in a particular simulation;
LAMMPS ignores those entries.

For a single-element simulation, only a single entry is required
(e.g. SiSiSi).  For a two-element simulation, the file must contain 8
entries (for SiSiSi, SiSiC, SiCSi, SiCC, CSiSi, CSiC, CCSi, CCC), that
specify SW parameters for all permutations of the two elements
interacting in three-body configurations.  Thus for 3 elements, 27
entries would be required, etc.

As annotated above, the first element in the entry is the center atom in
a three-body interaction.  Thus an entry for SiCC means a Si atom with 2
C atoms as neighbors.  The parameter values used for the two-body
interaction come from the entry where the second and third elements are
the same.  Thus the two-body parameters for Si interacting with C, comes
from the SiCC entry.  The three-body angular potential
:math:`f^{\textrm{3b}}\left(\theta_{ijk}\right)` can in principle be
specific to the three elements of the configuration. However, the user
must ensure that it makes physically sense.  Note also that the function
:math:`\phi_3` contains two exponential screening factors with parameter
values from the ij pair and ik pairs. So :math:`\phi_3` for a C atom
bonded to a Si atom and a second C atom will depend on the three-body
parameters for the CSiC entry, and also on the two-body parameters for
the CCC and CSiSi entries. Since the order of the two neighbors is
arbitrary, the three-body parameters and the tabulated angular potential
for entries CSiC and CCSi should be the same.  Similarly, the two-body
parameters for entries SiCC and CSiSi should also be the same.  The
parameters used only for two-body interactions (A, B, p, and q) in
entries whose second and third element are different (e.g. SiCSi) are
not used for anything and can be set to 0.0 if desired.  This is also
true for the parameters in :math:`\phi_3` that are taken from the ij and
ik pairs (:math:`\sigma`, *a*, :math:`\gamma`)

Additional input files and reference data can be found at:
https://gitlab.mpcdf.mpg.de/votca/votca/-/tree/master/csg-tutorials/spce/3body_sw

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as described
above from values in the potential file, but not for the tabulated
angular potential file.

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

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style sw <pair_sw>`,
:doc:`pair_style threebody/table <pair_threebody_table>`


----------

.. _Stillinger3:

**(Stillinger)** Stillinger and Weber, Phys Rev B, 31, 5262 (1985).

.. _Scherer1:

**(Scherer1)** C. Scherer and D. Andrienko, Phys. Chem. Chem. Phys. 20, 22387-22394 (2018).

