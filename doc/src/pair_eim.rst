.. index:: pair_style eim
.. index:: pair_style eim/omp

pair_style eim command
======================

Accelerator Variants: *eim/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *eim*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style eim
   pair_coeff * * Na Cl ../potentials/ffield.eim Na Cl
   pair_coeff * * Na Cl ffield.eim  Na Na Na Cl
   pair_coeff * * Na Cl ../potentials/ffield.eim Cl NULL Na

Description
"""""""""""

Style *eim* computes pairwise interactions for ionic compounds
using embedded-ion method (EIM) potentials :ref:`(Zhou) <Zhou2>`.  The
energy of the system E is given by

.. math::

   E = \frac{1}{2} \sum_{i=1}^{N} \sum_{j=i_1}^{i_N} \phi_{ij} \left(r_{ij}\right) + \sum_{i=1}^{N}E_i\left(q_i,\sigma_i\right)

The first term is a double pairwise sum over the J neighbors of all I
atoms, where :math:`\phi_{ij}` is a pair potential.  The second term sums over
the embedding energy E_i of atom I, which is a function of its charge
q_i and the electrical potential :math:`\sigma_i` at its location.  E_i, q_i,
and :math:`sigma_i` are calculated as

.. math::

   q_i  = & \sum_{j=i_1}^{i_N} \eta_{ji}\left(r_{ij}\right) \\
   \sigma_i  = & \sum_{j=i_1}^{i_N} q_j \cdot \psi_{ij} \left(r_{ij}\right) \\
   E_i\left(q_i,\sigma_i\right)  = & \frac{1}{2} \cdot q_i \cdot \sigma_i

where :math:`\eta_{ji}` is a pairwise function describing electron flow from atom
I to atom J, and :math:`\psi_{ij}` is another pairwise function.  The multi-body
nature of the EIM potential is a result of the embedding energy term.
A complete list of all the pair functions used in EIM is summarized
below

.. math::

   \phi_{ij}\left(r\right) = & \left\{ \begin{array}{lr}
   \left[\frac{E_{b,ij}\beta_{ij}}{\beta_{ij}-\alpha_{ij}}\exp\left(-\alpha_{ij} \frac{r-r_{e,ij}}{r_{e,ij}}\right)-\frac{E_{b,ij}\alpha_{ij}}{\beta_{ij}-\alpha_{ij}}\exp\left(-\beta_{ij} \frac{r-r_{e,ij}}{r_{e,ij}}\right)\right]f_c\left(r,r_{e,ij},r_{c,\phi,ij}\right),& p_{ij}=1 \\
   \left[\frac{E_{b,ij}\beta_{ij}}{\beta_{ij}-\alpha_{ij}} \left(\frac{r_{e,ij}}{r}\right)^{\alpha_{ij}}  -\frac{E_{b,ij}\alpha_{ij}}{\beta_{ij}-\alpha_{ij}} \left(\frac{r_{e,ij}}{r}\right)^{\beta_{ij}}\right]f_c\left(r,r_{e,ij},r_{c,\phi,ij}\right),& p_{ij}=2
   \end{array}
   \right.\\
   \eta_{ji} = & A_{\eta,ij}\left(\chi_j-\chi_i\right)f_c\left(r,r_{s,\eta,ij},r_{c,\eta,ij}\right) \\
   \psi_{ij}\left(r\right) = & A_{\psi,ij}\exp\left(-\zeta_{ij}r\right)f_c\left(r,r_{s,\psi,ij},r_{c,\psi,ij}\right) \\
   f_{c}\left(r,r_p,r_c\right) = & 0.510204 \cdot \mathrm{erfc}\left[\frac{1.64498\left(2r-r_p-r_c\right)}{r_c-r_p}\right] - 0.010204

Here :math:`E_b, r_e, r_(c,\phi), \alpha, \beta, A_(\psi), \zeta, r_(s,\psi),
r_(c,\psi), A_(\eta), r_(s,\eta), r_(c,\eta), \chi,` and pair function type
*p* are parameters, with subscripts *ij* indicating the two species of
atoms in the atomic pair.

.. note::

   Even though the EIM potential is treating atoms as charged ions,
   you should not use a LAMMPS :doc:`atom_style <atom_style>` that stores a
   charge on each atom and thus requires you to assign a charge to each
   atom, e.g. the *charge* or *full* atom styles.  This is because the
   EIM potential infers the charge on an atom from the equation above for
   q_i; you do not assign charges explicitly.

----------

All the EIM parameters are listed in a potential file which is
specified by the :doc:`pair_coeff <pair_coeff>` command.  This is an
ASCII text file in a format described below.  The "ffield.eim" file
included in the "potentials" directory of the LAMMPS distribution
currently includes nine elements Li, Na, K, Rb, Cs, F, Cl, Br, and I.
A system with any combination of these elements can be modeled.  This
file is parameterized in terms of LAMMPS :doc:`metal units <units>`.

Note that unlike other potentials, cutoffs for EIM potentials are not
set in the pair_style or pair_coeff command; they are specified in the
EIM potential file itself.  Likewise, the EIM potential file lists
atomic masses; thus you do not need to use the :doc:`mass <mass>`
command to specify them.

Only a single pair_coeff command is used with the *eim* style which
specifies an EIM potential file and the element(s) to extract
information for.  The EIM elements are mapped to LAMMPS atom types by
specifying N additional arguments after the filename in the pair_coeff
command, where N is the number of LAMMPS atom types:

* Elem1, Elem2, ...
* EIM potential file
* N element names = mapping of EIM elements to atom types

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the potential file.

As an example like one of those above, suppose you want to model a
system with Na and Cl atoms.  If your LAMMPS simulation has 4 atoms
types and you want the first 3 to be Na, and the fourth to be Cl, you would
use the following pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * Na Cl ffield.eim Na Na Na Cl

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The filename is the EIM potential file.  The Na and Cl arguments
(before the file name) are the two elements for which info will be
extracted from the potential file.  The first three trailing Na
arguments map LAMMPS atom types 1,2,3 to the EIM Na element.  The
final Cl argument maps LAMMPS atom type 4 to the EIM Cl element.

If a mapping value is specified as NULL, the mapping is not performed.
This can be used when an *eim* potential is used as part of the
*hybrid* pair style.  The NULL values are placeholders for atom types
that will be used with other potentials.

The ffield.eim file in the *potentials* directory of the LAMMPS
distribution is formatted as follows:

Lines starting with # are comments and are ignored by LAMMPS.  Lines
starting with "global:" include three global values. The first value
divides the cations from anions, i.e., any elements with
electronegativity above this value are viewed as anions, and any
elements with electronegativity below this value are viewed as
cations. The second and third values are related to the cutoff
function - i.e. the 0.510204, 1.64498, and 0.010204 shown in the above
equation can be derived from these values.

Lines starting with "element:" are formatted as follows: name of
element, atomic number, atomic mass, electronic negativity, atomic
radius (LAMMPS ignores it), ionic radius (LAMMPS ignores it), cohesive
energy (LAMMPS ignores it), and q0 (must be 0).

Lines starting with "pair:" are entered as: element 1, element 2,
r_(c,phi), r_(c,phi) (redundant for historical reasons), E_b, r_e,
alpha, beta, r_(c,eta), A_(eta), r_(s,eta), r_(c,psi), A_(psi), zeta,
r_(s,psi), and p.

The lines in the file can be in any order; LAMMPS extracts the info it
needs.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This style is part of the MANYBODY package.  It is only enabled if
LAMMPS was built with that package.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Zhou2:

**(Zhou)** Zhou, submitted for publication (2010).  Please contact
Xiaowang Zhou (Sandia) for details via email at xzhou at sandia.gov.
