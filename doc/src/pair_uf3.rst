.. index:: pair_style uf3
.. index:: pair_style uf3/kk

pair_style uf3 command
======================

Accelerator Variants: *uf3/kk*

Syntax
""""""

.. code-block:: LAMMPS

    pair_style style BodyFlag

* style = *uf3* or *uf3/kk*

  .. parsed-literal::

       BodyFlag = Indicates whether to calculate only 2-body or 2 and 3-body interactions. Possible values: 2 or 3

Examples
""""""""

.. code-block:: LAMMPS

    pair_style uf3 3
    pair_coeff * * Nb.uf3 Nb

    pair_style uf3 2
    pair_coeff * * NbSn.uf3 Nb Sn

    pair_style uf3 3
    pair_coeff * * NbSn.uf3 Nb Sn

Description
"""""""""""

.. versionadded:: 27June2024

The *uf3* style computes the :ref:`Ultra-Fast Force Fields (UF3)
<Xie23>` potential, a machine-learning interatomic potential. In UF3,
the total energy of the system is defined via two- and three-body
interactions:

.. math::

   E & = \sum_{i,j} V_2(r_{ij}) + \sum_{i,j,k} V_3 (r_{ij},r_{ik},r_{jk}) \\
   V_2(r_{ij}) & = \sum_{n=0}^N c_n B_n(r_{ij}) \\
   V_3 (r_{ij},r_{ik},r_{jk}) & = \sum_{l=0}^{N_l} \sum_{m=0}^{N_m} \sum_{n=0}^{N_n} c_{l,m,n} B_l(r_{ij}) B_m(r_{ik}) B_n(r_{jk})

where :math:`V_2(r_{ij})` and :math:`V_3 (r_{ij},r_{ik},r_{jk})` are the
two- and three-body interactions, respectively. For the two-body the
summation is over all neighbors J and for the three-body the summation
is over all neighbors J and K of atom I within a cutoff distance
determined from the potential files. :math:`B_n(r_{ij})` are the cubic
b-spline basis, :math:`c_n` and :math:`c_{l,m,n}` are the machine-learned
interaction parameters and :math:`N`, :math:`N_l`, :math:`N_m`, and
:math:`N_n` denote the number of basis functions per spline or tensor
spline dimension.

With *uf3* style only a single pair_coeff command is used to indicate the
UF3 LAMMPS potential file containing all the two- and three-body interactions
followed by N additional arguments specifying the mapping of UF3 elements to
LAMMPS atom types, where N is the number of LAMMPS atom types:

* UF3 LAMMPS potential file
* N elements names = mapping of UF3 elements to atom types

As an example, if a LAMMPS simulation contains 2 atom types (elements
'A' and 'B'), the pair_coeff command will be:

.. code-block:: LAMMPS

   pair_style uf3 3
   pair_coeff * * AB.uf3 A B

The AB.uf3 file should contain all two-body (A-A, A-B, B-B) and three-body
(A-A-A, A-A-B, A-B-B, B-A-A, B-A-B, B-B-B).

If a value of "2" is specified in the :code:`pair_style uf3` command,
only the two-body potentials are needed. For 3-body interaction the
first atom type is the central atom. We recommend using the
:code:`generate_uf3_lammps_pots.py` script (found `here
<https://github.com/uf3/uf3/tree/develop/lammps_plugin/scripts>`_) for
generating the UF3 LAMMPS potential file from the UF3 JSON potentials.

----------

UF3 LAMMPS potential file in the *potentials* directory of the LAMMPS
distribution have a ".uf3" suffix. The interaction block in UF3 LAMMPS potential
file should start with :code:`#UF3 POT` and end with :code:`#` characters.
Following shows the format of a generic 2-body and 3-body potential block in
UF3 LAMMPS potential file-

.. code-block:: LAMMPS

   #UF3 POT UNITS: units DATE: POT_GEN_DATE AUTHOR: AUTHOR_NAME CITATION: CITE
   2B ELEMENT1 ELEMENT2 LEADING_TRIM TRAILING_TRIM
   Rij_CUTOFF NUM_OF_KNOTS
   BSPLINE_KNOTS
   NUM_OF_COEFF
   COEFF
   #
   #UF3 POT UNITS: units DATE: POT_GEN_DATE AUTHOR: AUTHOR_NAME CITATION: CITE
   3B ELEMENT1 ELEMENT2 ELEMENT3 LEADING_TRIM TRAILING_TRIM
   Rjk_CUTOFF Rik_CUTOFF Rij_CUTOFF NUM_OF_KNOTS_JK NUM_OF_KNOTS_IK NUM_OF_KNOTS_IJ
   BSPLINE_KNOTS_FOR_JK
   BSPLINE_KNOTS_FOR_IK
   BSPLINE_KNOTS_FOR_IJ
   SHAPE_OF_COEFF_MATRIX[I][J][K]
   COEFF_MATRIX[0][0][K]
   COEFF_MATRIX[0][1][K]
   COEFF_MATRIX[0][2][K]
   .
   .
   .
   COEFF_MATRIX[1][0][K]
   COEFF_MATRIX[1][1][K]
   COEFF_MATRIX[1][2][K]
   .
   .
   .
   #

The second line indicates whether the block contains data for 2-body
(:code:`2B`) or 3-body (:code:`3B`) interaction. This is followed by element
combination interaction, :code:`LEADING_TRIM` and :code:`TRAILING_TRIM`
number on the same line. The current implementation is only tested for
:code:`LEADING_TRIM=0` and :code:`TRAILING_TRIM=3`.
If other values are used LAMMPS is terminated after issuing an error message.
The :code:`Rij_CUTOFF` sets the 2-body cutoff for the interaction described
by the potential block. :code:`NUM_OF_KNOTS` is the number of knots
(or the length of the knot vector) present on the very next line. The
:code:`BSPLINE_KNOTS` line should contain all the knots in ascending order.
:code:`NUM_OF_COEFF` is the number of coefficients in the :code:`COEFF` line.
All the numbers in the BSPLINE_KNOTS and COEFF line should be space-separated.
Similar to the 2-body potential block, the third line sets the cutoffs and
length of the knots. The cutoff distance between atom-type I and J is
:code:`Rij_CUTOFF`, atom-type I and K is :code:`Rik_CUTOFF` and between
J and K is :code:`Rjk_CUTOFF`.

.. note::

   The current implementation only works for UF3 potentials with cutoff
   distances for 3-body interactions that follows
   :code:`2Rij_CUTOFF=2Rik_CUTOFF=Rjk_CUTOFF` relation.

The :code:`BSPLINE_KNOTS_FOR_JK`, :code:`BSPLINE_KNOTS_FOR_IK`, and
:code:`BSPLINE_KNOTS_FOR_IJ` lines (note the order) contain the knots in
increasing order for atoms J and K, I and K, and atoms I and J
respectively. The number of knots is defined by the
:code:`NUM_OF_KNOTS_*` characters in the previous line.  The shape of
the coefficient matrix is defined on the
:code:`SHAPE_OF_COEFF_MATRIX[I][J][K]` line followed by the columns of
the coefficient matrix, one per line, as shown above. For example, if
the coefficient matrix has the shape of 8x8x13, then
:code:`SHAPE_OF_COEFF_MATRIX[I][J][K]` will be :code:`8 8 13` followed
by 64 (8x8) lines each containing 13 coefficients separated by space.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as described
above from values in the potential file.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`, since it is stored in potential file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

The 'uf3' pair style is part of the ML-UF3 package. It is only enabled
if LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on".

The UF3 LAMMPS potential file provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`.

The single() function of 'uf3' pair style only return the 2-body
interaction energy.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Xie23:

**(Xie23)** Xie, S.R., Rupp, M. & Hennig, R.G. Ultra-fast interpretable machine-learning potentials. npj Comput Mater 9, 162 (2023). https://doi.org/10.1038/s41524-023-01092-7
