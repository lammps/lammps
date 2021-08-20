.. index:: pair_style tersoff/mod
.. index:: pair_style tersoff/mod/c
.. index:: pair_style tersoff/mod/gpu
.. index:: pair_style tersoff/mod/kk
.. index:: pair_style tersoff/mod/omp
.. index:: pair_style tersoff/mod/c/omp

pair_style tersoff/mod command
==============================

Accelerator Variants: *tersoff/mod/gpu*, *tersoff/mod/kk*, *tersoff/mod/omp*

pair_style tersoff/mod/c command
================================

Accelerator Variants: *tersoff/mod/c/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style keywords values

* style = *tersoff/mod* or *tersoff/mod/c*
* keyword = *shift*

  .. parsed-literal::

       *shift* value = delta
         delta = negative shift in equilibrium bond length

Examples
""""""""

.. code-block:: LAMMPS

   pair_style tersoff/mod
   pair_coeff * * Si.tersoff.mod Si Si

   pair_style tersoff/mod/c
   pair_coeff * * Si.tersoff.modc Si Si

Description
"""""""""""

The *tersoff/mod* and *tersoff/mod/c* styles computes a bond-order type
interatomic potential :ref:`(Kumagai) <Kumagai>` based on a 3-body Tersoff
potential :ref:`(Tersoff_1) <Tersoff_12>`, :ref:`(Tersoff_2) <Tersoff_22>` with
modified cutoff function and angular-dependent term, giving the energy
E of a system of atoms as

.. math::

   E & = \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij} \\
   V_{ij} & = f_C(r_{ij} + \delta) \left[ f_R(r_{ij} + \delta) + b_{ij} f_A(r_{ij} + \delta) \right] \\
   f_C(r) & = \left\{ \begin{array} {r@{\quad:\quad}l}
     1 & r < R - D \\
     \frac{1}{2} - \frac{9}{16} \sin \left( \frac{\pi}{2} \frac{r-R}{D} \right) - \frac{1}{16} \sin \left( \frac{3\pi}{2} \frac{r-R}{D} \right) &
       R-D < r < R + D \\
     0 & r > R + D
     \end{array} \right. \\
   f_R(r) & = A \exp (-\lambda_1 r) \\
   f_A(r) & = -B \exp (-\lambda_2 r) \\
   b_{ij} & = \left( 1 + {\zeta_{ij}}^\eta \right)^{-\frac{1}{2n}} \\
   \zeta_{ij} & = \sum_{k \neq i,j} f_C(r_{ik} + \delta) g(\theta_{ijk})
                    \exp \left[ \alpha (r_{ij} - r_{ik})^\beta \right] \\
   g(\theta) & = c_1 + g_o(\theta) g_a(\theta) \\
   g_o(\theta) & = \frac{c_2 (h - \cos \theta)^2}{c_3 + (h - \cos \theta)^2} \\
   g_a(\theta) & = 1 + c_4 \exp \left[ -c_5 (h - \cos \theta)^2 \right] \\

where :math:`f_R` is a two-body term and :math:`f_A` includes three-body interactions.
:math:`\delta` is an optional negative shift of the
equilibrium bond length, as described below.

The summations in the formula are over all neighbors J and K of atom I
within a cutoff distance = R + D.
The *tersoff/mod/c* style differs from *tersoff/mod* only in the
formulation of the V_ij term, where it contains an additional c0 term.

.. math::

   V_{ij} = f_C(r_{ij} + \delta) \left[ f_R(r_{ij} + \delta) + b_{ij} f_A(r_{ij} + \delta) + c_0 \right] \\

The modified cutoff function :math:`f_C` proposed by :ref:`(Murty) <Murty>` and
having a continuous second-order differential is employed. The
angular-dependent term :math:`g(\theta)` was modified to increase the
flexibility of the potential.

The *tersoff/mod* potential is fitted to both the elastic constants
and melting point by employing the modified Tersoff potential function
form in which the angular-dependent term is improved. The model
performs extremely well in describing the crystalline, liquid, and
amorphous phases :ref:`(Schelling) <Schelling>`.

Only a single pair_coeff command is used with the *tersoff/mod* style
which specifies a Tersoff/MOD potential file with parameters for all
needed elements.  These are mapped to LAMMPS atom types by specifying
N additional arguments after the filename in the pair_coeff command,
where N is the number of LAMMPS atom types:

* filename
* N element names = mapping of Tersoff/MOD elements to atom types

As an example, imagine the Si.tersoff_mod file has Tersoff values for Si.
If your LAMMPS simulation has 3 Si atoms types, you would use the following
pair_coeff command:

.. code-block:: LAMMPS

   pair_coeff * * Si.tersoff_mod Si Si Si

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The three Si arguments map LAMMPS atom types 1,2,3 to the Si element
in the Tersoff/MOD file. If a mapping value is specified as NULL, the
mapping is not performed.  This can be used when a *tersoff/mod*
potential is used as part of the *hybrid* pair style. The NULL values
are placeholders for atom types that will be used with other
potentials.

Tersoff/MOD file in the *potentials* directory of the LAMMPS
distribution have a ".tersoff.mod" suffix. Potential files for the
*tersoff/mod/c* style have the suffix ".tersoff.modc". Lines that are
not blank or comments (starting with #) define parameters for a triplet
of elements.  The parameters in a single entry correspond to
coefficients in the formulae above:

* element 1 (the center atom in a 3-body interaction)
* element 2 (the atom bonded to the center atom)
* element 3 (the atom influencing the 1-2 bond in a bond-order sense)
* :math:`\beta`
* :math:`\alpha`
* h
* :math:`\eta`
* :math:`\beta_{ters}` = 1 (dummy parameter)
* :math:`\lambda_2` (1/distance units)
* B (energy units)
* R (distance units)
* D (distance units)
* :math:`\lambda_1` (1/distance units)
* A (energy units)
* n
* c1
* c2
* c3
* c4
* c5
* c0 (energy units, tersoff/mod/c only):ul

The n, :math:`\eta`, :math:`\lambda_2`, B, :math:`\lambda_1`, and A parameters are only used for
two-body interactions.  The :math:`\beta`, :math:`\alpha`, c1, c2, c3, c4, c5, h
parameters are only used for three-body interactions. The R and D
parameters are used for both two-body and three-body interactions.
The c0 term applies to *tersoff/mod/c* only. The non-annotated
parameters are unitless.

The Tersoff/MOD potential file must contain entries for all the elements
listed in the pair_coeff command.  It can also contain entries for
additional elements not being used in a particular simulation; LAMMPS
ignores those entries.

For a single-element simulation, only a single entry is required
(e.g. SiSiSi). As annotated above, the first element in the entry is
the center atom in a three-body interaction and it is bonded to the
second atom and the bond is influenced by the third atom.  Thus an entry
for SiSiSi means Si bonded to a Si with another Si atom influencing the bond.

The *shift* keyword computes the energy E of a system of atoms, whose formula
is the same as the Tersoff potential. The only modification is that the original
equilibrium bond length ( :math:`r_0`) of the system is shifted to :math:`r_0-\delta`.
The minus sign arises because each radial distance :math:`r` is replaced by :math:`r+\delta`.
More information on this option is given on the main :doc:`pair_tersoff <pair_tersoff>` page.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the MANYBODY package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This pair style requires the :doc:`newton <newton>` setting to be "on"
for pair interactions.

The *shift* keyword is not supported by the *tersoff/gpu*,
*tersoff/intel*, *tersoff/kk*, *tersoff/table* or *tersoff/table/omp*
variants.

The *tersoff/mod* potential files provided with LAMMPS (see the potentials
directory) are parameterized for metal :doc:`units <units>`.  You can
use the *tersoff/mod* pair style with any LAMMPS units, but you would need to
create your own Tersoff/MOD potential file with coefficients listed in the
appropriate units if your simulation does not use "metal" units.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Kumagai:

**(Kumagai)** T. Kumagai, S. Izumi, S. Hara, S. Sakai,
Comp. Mat. Science, 39, 457 (2007).

.. _Tersoff_12:

**(Tersoff_1)** J. Tersoff, Phys Rev B, 37, 6991 (1988).

.. _Tersoff_22:

**(Tersoff_2)** J. Tersoff, Phys Rev B, 38, 9902 (1988).

.. _Murty:

**(Murty)** M.V.R. Murty, H.A. Atwater, Phys Rev B, 51, 4889 (1995).

.. _Schelling:

**(Schelling)** Patrick K. Schelling, Comp. Mat. Science, 44, 274 (2008).
