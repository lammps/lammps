.. index:: pair_style oxrna2/excv
.. index:: pair_style oxrna2/stk
.. index:: pair_style oxrna2/hbond
.. index:: pair_style oxrna2/xstk
.. index:: pair_style oxrna2/coaxstk
.. index:: pair_style oxrna2/dh

pair_style oxrna2/excv command
==============================

pair_style oxrna2/stk command
=============================

pair_style oxrna2/hbond command
===============================

pair_style oxrna2/xstk command
==============================

pair_style oxrna2/coaxstk command
=================================

pair_style oxrna2/dh command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style1

   pair_coeff * * style2 args

* style1 = *hybrid/overlay oxrna2/excv oxrna2/stk oxrna2/hbond oxrna2/xstk oxrna2/coaxstk oxrna2/dh*

* style2 = *oxrna2/excv* or *oxrna2/stk* or *oxrna2/hbond* or *oxrna2/xstk* or *oxrna2/coaxstk* or *oxrna2/dh*
* args = list of arguments for these particular styles

.. parsed-literal::

     *oxrna2/stk* args = seq T xi kappa 6.0 0.43 0.93 0.35 0.78 0.9 0 0.95 0.9 0 0.95 1.3 0 0.8 1.3 0 0.8 2.0 0.65 2.0 0.65
       seq = seqav (for average sequence stacking strength) or seqdep (for sequence-dependent stacking strength)
       T = temperature (LJ units: 0.1 = 300 K, real units: 300 = 300 K)
       xi = 1.40206 (LJ units) or 8.35864576375849 (real units), temperature-independent coefficient in stacking strength
       kappa = 2.77 (LJ units) or 0.005504556 (real units), coefficient of linear temperature dependence in stacking strength
     *oxrna2/hbond* args = seq eps 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
       seq = seqav (for average sequence base-pairing strength) or seqdep (for sequence-dependent base-pairing strength)
       eps = 0.870439 (LJ units) or 5.18928666388042 (real units), average hydrogen bonding strength between A-U and C-G Watson-Crick and G-U wobble base pairs, 0 between all other pairs
     *oxrna2/dh* args = T rhos qeff
       T = temperature (LJ units: 0.1 = 300 K, real units: 300 = 300 K)
       rhos = salt concentration (mole per litre)
       qeff = 1.02455 (effective charge in elementary charges)

Examples
""""""""

.. code-block:: LAMMPS

   # LJ units
   pair_style hybrid/overlay oxrna2/excv oxrna2/stk oxrna2/hbond oxrna2/xstk oxrna2/coaxstk oxrna2/dh
   pair_coeff * * oxrna2/excv    2.0 0.7 0.675 2.0 0.515 0.5 2.0 0.33 0.32
   pair_coeff * * oxrna2/stk     seqdep 0.1 1.40206 2.77 6.0 0.43 0.93 0.35 0.78 0.9 0 0.95 0.9 0 0.95 1.3 0 0.8 1.3 0 0.8 2.0 0.65 2.0 0.65
   pair_coeff * * oxrna2/hbond   seqdep 0.0 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff 1 4 oxrna2/hbond   seqdep 0.870439 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff 2 3 oxrna2/hbond   seqdep 0.870439 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff 3 4 oxrna2/hbond   seqdep 0.870439 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff * * oxrna2/xstk    59.9626 0.5 0.6 0.42 0.58 2.25 0.505 0.58 1.7 1.266 0.68 1.7 1.266 0.68 1.7 0.309 0.68 1.7 0.309 0.68
   pair_coeff * * oxrna2/coaxstk 80 0.5 0.6 0.42 0.58 2.0 2.592 0.65 1.3 0.151 0.8 0.9 0.685 0.95 0.9 0.685 0.95 2.0 -0.65 2.0 -0.65
   pair_coeff * * oxrna2/dh      0.1 0.5 1.02455

   pair_style hybrid/overlay oxrna2/excv oxrna2/stk oxrna2/hbond oxrna2/xstk oxrna2/coaxstk oxrna2/dh
   pair_coeff * * oxrna2/excv    oxrna2_lj.cgdna
   pair_coeff * * oxrna2/stk     seqdep 0.1 oxrna2_lj.cgdna
   pair_coeff * * oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff 1 4 oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff 2 3 oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff 3 4 oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff * * oxrna2/xstk    oxrna2_lj.cgdna
   pair_coeff * * oxrna2/coaxstk oxrna2_lj.cgdna
   pair_coeff * * oxrna2/dh      0.1 0.5 oxrna2_lj.cgdna

   # Real units
   pair_style hybrid/overlay oxrna2/excv oxrna2/stk oxrna2/hbond oxrna2/xstk oxrna2/coaxstk oxrna2/dh
   pair_coeff * * oxrna2/excv    11.92337812042065 5.9626 5.74965 11.92337812042065 4.38677 4.259 11.92337812042065 2.81094 2.72576
   pair_coeff * * oxrna2/stk     seqdep 300.0 8.35864576375849 0.005504556 0.70439070204273 3.66274 7.92174 2.9813 6.64404 0.9 0.0 0.95 0.9 0.0 0.95 1.3 0.0 0.8 1.3 0.0 0.8 2.0 0.65 2.0 0.65
   pair_coeff * * oxrna2/hbond   seqdep 0.0 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff 1 4 oxrna2/hbond   seqdep 5.18928666388042 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff 2 3 oxrna2/hbond   seqdep 5.18928666388042 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff 3 4 oxrna2/hbond   seqdep 5.18928666388042 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45
   pair_coeff * * oxrna2/xstk    4.92690859644113 4.259 5.1108 3.57756 4.94044 2.25 0.505 0.58 1.7 1.266 0.68 1.7 1.266 0.68 1.7 0.309 0.68 1.7 0.309 0.68
   pair_coeff * * oxrna2/coaxstk 6.57330882442206 4.259 5.1108 3.57756 4.94044 2.0 2.592 0.65 1.3 0.151 0.8 0.9 0.685 0.95 0.9 0.685 0.95 2.0 -0.65 2.0 -0.65
   pair_coeff * * oxrna2/dh      300.0 0.5 1.02455

   pair_style hybrid/overlay oxrna2/excv oxrna2/stk oxrna2/hbond oxrna2/xstk oxrna2/coaxstk oxrna2/dh
   pair_coeff * * oxrna2/excv    oxrna2_real.cgdna
   pair_coeff * * oxrna2/stk     seqdep 300.0 oxrna2_real.cgdna
   pair_coeff * * oxrna2/hbond   seqdep oxrna2_real.cgdna
   pair_coeff 1 4 oxrna2/hbond   seqdep oxrna2_real.cgdna
   pair_coeff 2 3 oxrna2/hbond   seqdep oxrna2_real.cgdna
   pair_coeff 3 4 oxrna2/hbond   seqdep oxrna2_real.cgdna
   pair_coeff * * oxrna2/xstk    oxrna2_real.cgdna
   pair_coeff * * oxrna2/coaxstk oxrna2_real.cgdna
   pair_coeff * * oxrna2/dh      300.0 0.5 oxrna2_real.cgdna

.. note::

   The coefficients in the above examples are provided in forms
   compatible with both *units lj* and *units real* (see documentation
   of :doc:`units <units>`).  These can also be read from a potential
   file with correct unit style by specifying the name of the
   file. Several potential files for each unit style are included in the
   ``potentials`` directory of the LAMMPS distribution.

Description
"""""""""""

The *oxrna2* pair styles compute the pairwise-additive parts of the
oxDNA force field for coarse-grained modelling of RNA. The effective
interaction between the nucleotides consists of potentials for the
excluded volume interaction *oxrna2/excv*, the stacking *oxrna2/stk*,
cross-stacking *oxrna2/xstk* and coaxial stacking interaction
*oxrna2/coaxstk*, electrostatic Debye-Hueckel interaction *oxrna2/dh* as
well as the hydrogen-bonding interaction *oxrna2/hbond* between
complementary pairs of nucleotides on opposite strands. Average sequence
or sequence-dependent stacking and base-pairing strengths are supported
:ref:`(Sulc2) <Sulc32>`.

The exact functional form of the pair styles is rather complex.  The
individual potentials consist of products of modulation factors, which
themselves are constructed from a number of more basic potentials
(Morse, Lennard-Jones, harmonic angle and distance) as well as quadratic
smoothing and modulation terms.  We refer to :ref:`(Sulc1) <Sulc31>` and
the original oxDNA publications :ref:`(Ouldridge-DPhil)
<Ouldridge-DPhil3>` and :ref:`(Ouldridge) <Ouldridge3>` for a detailed
description of the oxRNA2 force field.

.. note::

   These pair styles have to be used together with the related oxDNA2
   bond style *oxrna2/fene* for the connectivity of the phosphate
   backbone (see also documentation of :doc:`bond_style oxrna2/fene
   <bond_oxdna>`). Most of the coefficients in the above example have to
   be kept fixed and cannot be changed without reparameterizing the
   entire model.  Exceptions are the first two coefficients after
   *oxrna2/stk* (seq=seqdep and T=0.1 and
   corresponding *real unit* equivalents in the above examples), the
   first coefficient after *oxrna2/hbond* (seq=seqdep in the above
   example) and the two coefficients after *oxrna2/dh* (T=0.1 and
   rhos=0.5 in the above example). When using a Langevin
   thermostat e.g. through :doc:`fix langevin <fix_langevin>` or
   :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>` the temperature
   coefficients have to be matched to the one used in the fix.

.. note::

   These pair styles have to be used with the *atom_style hybrid bond
   ellipsoid oxdna* (see documentation of :doc:`atom_style
   <atom_style>`). The *atom_style oxdna* stores the 3'-to-5' polarity
   of the nucleotide strand, which is set through the bond topology in
   the data file. The first (second) atom in a bond definition is
   understood to point towards the 3'-end (5'-end) of the strand.

.. warning::

   If data files are produced with :doc:`write_data <write_data>`, then
   the :doc:`newton <newton>` command should be set to *newton on*.
   Otherwise the data files will not have the same 3'-to-5' polarity
   as the initial data file. This limitation does not apply to
   binary restart files produced with :doc:`write_restart <write_restart>`.

Example input and data files for DNA duplexes can be found in
``examples/PACKAGES/cgdna/examples/lj_units/oxRNA2/`` or in the
corresponding folder for real units.
A simple python setup tool which creates single straight or helical
DNA strands, DNA duplexes or arrays of DNA duplexes can be found in
``examples/PACKAGES/cgdna/util/``.

Please cite :ref:`(Henrich) <Henrich3>` in any publication that uses
this implementation.  The article contains general information on the
model, its implementation and performance as well as the structure of
the data and input file. The preprint version of the article can be
found `here <PDF/CG-DNA.pdf>`_.  Please cite also the relevant oxRNA2
publications :ref:`(Sulc1) <Sulc31>` and :ref:`(Sulc2) <Sulc32>`.

----------

Potential file reading
""""""""""""""""""""""

For each pair style above the first non-modifiable argument can be a
filename (with exception of Debye-Hueckel, for which the effective
charge argument can be a filename), and if it is, no further arguments
should be supplied.  Therefore the following command:

.. code-block:: LAMMPS

   pair_coeff 3 4 oxrna2/hbond   seqdep oxrna2_lj.cgdna

will be interpreted as a request to read the corresponding hydrogen
bonding potential parameters from the file with the given name.  The
file can define multiple potential parameters for both bonded and pair
interactions, but for the example pair interaction above there must
exist in the file a line of the form:

.. code-block:: LAMMPS

  3 4 hbond     <coefficients>

If potential customization is required, the potential file reading can
be mixed with the manual specification of the potential parameters. For
example, the following command:

.. code-block:: LAMMPS

   pair_style hybrid/overlay oxrna2/excv oxrna2/stk oxrna2/hbond oxrna2/xstk oxrna2/coaxstk oxrna2/dh
   pair_coeff * * oxrna2/excv    2.0 0.7 0.675 2.0 0.515 0.5 2.0 0.33 0.32
   pair_coeff * * oxrna2/stk     seqdep 0.1 oxrna2_lj.cgdna
   pair_coeff * * oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff 1 4 oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff 2 3 oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff 3 4 oxrna2/hbond   seqdep oxrna2_lj.cgdna
   pair_coeff * * oxrna2/xstk    oxrna2_lj.cgdna
   pair_coeff * * oxrna2/coaxstk oxrna2_lj.cgdna
   pair_coeff * * oxrna2/dh      0.1 0.5 1.02455

will read the excluded volume and Debye-Hueckel effective charge *qeff*
parameters from the manual specification and all others from the
potential file *oxrna2_lj.cgdna*.

There are sample potential files for each unit style in the
``potentials`` directory of the LAMMPS distribution. The potential file
unit system must align with the units defined via the :doc:`units
<units>` command. For conversion between different *LJ* and *real* unit
systems for oxDNA, the python tool *lj2real.py* located in the
``examples/PACKAGES/cgdna/util/`` directory can be used. This tool
assumes similar file structure to the examples found in
``examples/PACKAGES/cgdna/examples/``.

----------

Unique base pairing
""""""""""""""""""""""

Unique base pairing describes the restriction on the specific complementary
nucleotide with which a particular base can pair. This can be used to prevent
asymmetric base pairs or to simplify the free energy landscape. With unique
base pairing enabled base pairs can only form between complementary nucleotides
with specific atom IDs. This functionality draws on :doc:`fix property/atom <fix_property_atom>`
and a modified :doc:`read_data <read_data>` command.

To use unique base pairing, the data file of a system with N nucleotides must contain a section like

.. code-block:: LAMMPS

   Basepairs # i_idc

   1 idc1
   2 idc2
   3 idc3
   4 idc4
   ...
   N idcN

where idc is the non-negative atom ID of a complementary nucleotide that binds uniquely
to the preceding atom ID.

Unique base pairing can be combined with normal base pairing by setting a zero or negative value for idc.
For instance, in a 4-mer with 8 nucleotides consisting of a ssDNA strand 3'-A-A-A-A-5' with atom IDs 3'-1-2-3-4-5'
and a complementary strand 5'-T-T-T-T-3' with atom IDs 5'-8-7-6-5-3' set up as

.. code-block:: LAMMPS

   Basepairs # i_idc

   1 8
   2 -1
   3 -1
   4 5
   5 4
   6 -1
   7 -1
   8 1

the A nucleotide with ID 1 can only hybridize with the T nucleotide with ID 8 and
the A nucleotide with ID 4 can only hybridize with the T nucleotide with ID 5,
whereas the A nucleotides with ID 2 and 3 can hybridize with either T nucleotide with ID 6 and 7.

The input file requires an instance of the :doc:`fix property/atom <fix_property_atom>` and a
:doc:`read_data <read_data>` command as follows:

.. code-block:: LAMMPS

   fix Basepairs all property/atom i_idc ghost yes
   read_data file fix Basepairs NULL Basepairs

where *file* is the name of the data file and the only modifiable argument.
An example input and data file for a dsDNA ring can be found in
``examples/PACKAGES/cgdna/examples/lj_units/oxDNA3/unique_bp``
or in the corresponding folder for real units.

----------

Restrictions
""""""""""""

These pair styles can only be used if LAMMPS was built with the
CG-DNA package and the MOLECULE and ASPHERE package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`bond_style oxrna2/fene <bond_oxdna>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`bond_style oxdna/fene <bond_oxdna>`, :doc:`pair_style oxdna/excv <pair_oxdna>`,
:doc:`bond_style oxdna2/fene <bond_oxdna>`, :doc:`pair_style oxdna2/excv <pair_oxdna2>`,
:doc:`atom_style oxdna <atom_style>`, :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`

Default
"""""""


none

----------

.. _Henrich3:

**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E. Ouldridge, Eur. Phys. J. E 41, 57 (2018).

.. _Sulc31:

**(Sulc1)** P. Sulc, F. Romano, T. E. Ouldridge, et al., J. Chem. Phys. 140, 235102 (2014).

.. _Sulc32:

**(Sulc2)** P. Sulc, F. Romano, T.E. Ouldridge, L. Rovigatti, J.P.K. Doye, A.A. Louis, J. Chem. Phys. 137, 135101 (2012).

.. _Ouldridge-DPhil3:

**(Ouldridge-DPhil)** T.E. Ouldridge, Coarse-grained modelling of DNA and DNA self-assembly, DPhil. University of Oxford (2011).

.. _Ouldridge3:

**(Ouldridge)** T.E. Ouldridge, A.A. Louis, J.P.K. Doye, J. Chem. Phys. 134, 085101 (2011).
