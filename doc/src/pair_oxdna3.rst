.. index:: pair_style oxdna3/excv
.. index:: pair_style oxdna3/stk
.. index:: pair_style oxdna3/hbond
.. index:: pair_style oxdna3/xstk
.. index:: pair_style oxdna3/coaxstk
.. index:: pair_style oxdna3/dh

pair_style oxdna3/excv command
==============================

pair_style oxdna3/stk command
=============================

pair_style oxdna3/hbond command
===============================

pair_style oxdna3/xstk command
==============================

pair_style oxdna3/coaxstk command
=================================

pair_style oxdna3/dh command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style1

   pair_coeff * * style2 args (keyword value)

* style1 = *hybrid/overlay oxdna3/excv oxdna3/stk oxdna3/hbond oxdna3/xstk oxdna3/coaxstk oxdna3/dh*
* style2 = *oxdna3/excv* or *oxdna3/stk* or *oxdna3/hbond* or *oxdna3/xstk* or *oxdna3/coaxstk* or *oxdna3/dh*
* args = list of arguments for these particular styles
* zero or one keyword/value pair may be appended to *oxdna3/dh*
* keyword = *half_charged_ends*

.. parsed-literal::

     *oxdna3/excv* args = oxdna3_lj.cgdna or oxdna3_real.cgdna
     *oxdna3/stk* args = T oxdna3_lj.cgdna or oxdna3_real.cgdna
       T = temperature (LJ units: 0.1 = 300 K, real units: 300 = 300 K)
     *oxdna3/hbond* args = oxdna3_lj.cgdna or oxdna3_real.cgdna
     *oxdna3/xstk* args = oxdna3_lj.cgdna or oxdna3_real.cgdna
     *oxdna3/coaxstk* args = oxdna3_lj.cgdna or oxdna3_real.cgdna
     *oxdna3/dh* args [keyword value] = T rhos oxdna3_lj.cgdna or oxdna3_real.cgdna [half_charged_ends no|yes]
       T = temperature (LJ units: 0.1 = 300 K, real units: 300 = 300 K)
       rhos = salt concentration (mole per litre)
       half_charged_ends yes = set half charge at terminal nucleotides
       half_charged_ends no  = set full charge at terminal nucleotides

Examples
""""""""

.. code-block:: LAMMPS

   # LJ units
   pair_style hybrid/overlay oxdna3/excv oxdna3/stk oxdna3/hbond oxdna3/xstk oxdna3/coaxstk oxdna3/dh
   pair_coeff * * oxdna3/excv     oxdna3_lj.cgdna
   pair_coeff * * oxdna3/stk      0.1 oxdna3_lj.cgdna
   pair_coeff * * oxdna3/hbond    oxdna3_lj.cgdna
   pair_coeff 1 4 oxdna3/hbond    oxdna3_lj.cgdna
   pair_coeff 2 3 oxdna3/hbond    oxdna3_lj.cgdna
   pair_coeff * * oxdna3/xstk     oxdna3_lj.cgdna
   pair_coeff * * oxdna3/coaxstk  oxdna3_lj.cgdna
   pair_coeff * * oxdna3/dh       0.1 0.2 oxdna3_lj.cgdna

   # Real units
   pair_style hybrid/overlay oxdna3/excv oxdna3/stk oxdna3/hbond oxdna3/xstk oxdna3/coaxstk oxdna3/dh
   pair_coeff * * oxdna3/excv     oxdna3_real.cgdna
   pair_coeff * * oxdna3/stk      300.0 oxdna3_real.cgdna
   pair_coeff * * oxdna3/hbond    oxdna3_real.cgdna
   pair_coeff 1 4 oxdna3/hbond    oxdna3_real.cgdna
   pair_coeff 2 3 oxdna3/hbond    oxdna3_real.cgdna
   pair_coeff * * oxdna3/xstk     oxdna3_real.cgdna
   pair_coeff * * oxdna3/coaxstk  oxdna3_real.cgdna
   pair_coeff * * oxdna3/dh       300.0 0.2 oxdna3_real.cgdna

.. note::

   The coefficients are provided in forms compatible with both
   *units lj* and *units real*. The potential file unit system
   must align with the units defined via the :doc:`units <units>` command.
   In case of oxDNA3 almost all coefficients have to be read from a potential
   file with correct unit style by specifying the name of the file. The
   potential files for each unit style are included in the ``potentials``
   directory of the LAMMPS distribution.


Description
"""""""""""

.. versionadded:: 30Mar2026

The *oxdna3* pair styles compute the pairwise-additive parts of the
oxDNA force field for coarse-grained modelling of DNA. The effective
interaction between the nucleotides consists of potentials for the
excluded volume interaction *oxdna3/excv*, the stacking *oxdna3/stk*,
cross-stacking *oxdna3/xstk* and coaxial stacking interaction
*oxdna3/coaxstk*, electrostatic Debye-Hueckel interaction *oxdna3/dh* as
well as the hydrogen-bonding interaction *oxdna3/hbond* between
complementary pairs of nucleotides on opposite strands.

The exact functional form of the pair styles is rather complex.  The
individual potentials consist of products of modulation factors, which
themselves are constructed from a number of more basic potentials
(Morse, Lennard-Jones, harmonic angle and distance) as well as quadratic
smoothing and modulation terms.  We refer to :ref:`(Bonato) <Bonato2>`
and the original oxDNA publications :ref:`(Ouldridge-DPhil)
<Ouldridge-DPhil4>` and :ref:`(Ouldridge) <Ouldridge4>`
for a detailed description of the oxDNA3 force field.

.. note::

   These pair styles have to be used together with the related oxDNA3
   bond style *oxdna3/fene* for the connectivity of the phosphate
   backbone (see also documentation of :doc:`bond_style oxdna3/fene
   <bond_oxdna>`). All coefficients in the above mentioned potential files
   have to be kept fixed and cannot be changed without reparameterizing the
   entire model.  The first coefficient after *oxdna3/stk*
   (T=0.1 and corresponding *real unit* equivalents in the above examples)
   and the two coefficients after *oxdna3/dh* (T=0.1 and rhos=0.2 in the
   above example) have to be set to the temperature and salt concentration
   of the system.
   *oxdna3/dh* has the option to set half a charge at terminal nucleotides
   (half_charged_ends yes) to aid coaxial stacking. When using a
   Langevin thermostat e.g. through :doc:`fix langevin <fix_langevin>` or
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
``examples/PACKAGES/cgdna/examples/lj_units/oxDNA3/`` or in the
corresponding folder for real units.
A simple python setup tool which creates single straight or helical DNA
strands, DNA duplexes or arrays of DNA duplexes can be found in
``examples/PACKAGES/cgdna/util/``.

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

Please cite :ref:`(Henrich) <Henrich6>` in any publication that uses
this implementation. An updated documentation that contains general
information on the model, its implementation and performance as well as
the structure of the data and input file can be found `here
<PDF/CG-DNA.pdf>`_.

Please cite also the relevant oxDNA3 publication :ref:`(Bonato) <Bonato2>`.

----------

Restrictions
""""""""""""

These pair styles can only be used if LAMMPS was built with the
CG-DNA package and the MOLECULE and ASPHERE package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`bond_style oxdna3/fene <bond_oxdna>`,
:doc:`bond_style oxdna/fene <bond_oxdna>`, :doc:`pair_style oxdna/excv <pair_oxdna>`,
:doc:`bond_style oxdna2/fene <bond_oxdna>`, :doc:`pair_style oxdna2/excv <pair_oxdna2>`,
:doc:`bond_style oxrna2/fene <bond_oxdna>`, :doc:`pair_style oxrna2/excv <pair_oxrna2>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`atom_style oxdna <atom_style>`,
:doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`

Default
"""""""

The option default is half_charged_ends = no.

----------

.. _Bonato2:

**(Bonato)** A. Bonato, T.E. Ouldridge, A.A. Louis, J.P.K. Doye, L. Rovigatti, M. Matthies, O.Henrich, in preparation.

.. _Ouldridge-DPhil4:

**(Ouldridge-DPhil)** T.E. Ouldridge, Coarse-grained modelling of DNA and DNA self-assembly, DPhil. University of Oxford (2011).

.. _Ouldridge4:

**(Ouldridge)** T.E. Ouldridge, A.A. Louis, J.P.K. Doye, J. Chem. Phys. 134, 085101 (2011).

.. _Henrich6:

**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E. Ouldridge, Eur. Phys. J. E 41, 57 (2018).

