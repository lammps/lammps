.. index:: pair_style oxrna2/excv

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
       T = temperature (oxDNA units, 0.1 = 300 K)
       xi = 1.40206 (temperature-independent coefficient in stacking strength)
       kappa = 2.77 (coefficient of linear temperature dependence in stacking strength)
     *oxrna2/hbond* args = seq eps 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45 
       seq = seqav (for average sequence base-pairing strength) or seqdep (for sequence-dependent base-pairing strength)
       eps = 0.870439 (between base pairs A-T, C-G and G-T) or 0 (all other pairs)
     *oxrna2/dh* args = T rhos qeff
       T = temperature (oxDNA units, 0.1 = 300 K)
       rhos = salt concentration (mole per litre)
       qeff = 1.02455 (effective charge in elementary charges)

Examples
""""""""


.. code-block:: LAMMPS

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

Description
"""""""""""

The *oxrna2* pair styles compute the pairwise-additive parts of the oxDNA force field
for coarse-grained modelling of DNA. The effective interaction between the nucleotides consists of potentials for the
excluded volume interaction *oxrna2/excv*\ , the stacking *oxrna2/stk*\ , cross-stacking *oxrna2/xstk*
and coaxial stacking interaction *oxrna2/coaxstk*\ , electrostatic Debye-Hueckel interaction *oxrna2/dh*
as well as the hydrogen-bonding interaction *oxrna2/hbond* between complementary pairs of nucleotides on
opposite strands. Average sequence or sequence-dependent stacking and base-pairing strengths
are supported :ref:`(Sulc2) <Sulc32>`. Quasi-unique base-pairing between nucleotides can be achieved by using 
more complementary pairs of atom types like 5-8 and 6-7, 9-12 and 10-11, 13-16 and 14-15, etc. 
This prevents the hybridization of in principle complementary bases within Ntypes/4 bases 
up and down along the backbone.

The exact functional form of the pair styles is rather complex.
The individual potentials consist of products of modulation factors,
which themselves are constructed from a number of more basic potentials
(Morse, Lennard-Jones, harmonic angle and distance) as well as quadratic smoothing and modulation terms.
We refer to :ref:`(Sulc1) <Sulc31>` and the original oxDNA publications :ref:`(Ouldridge-DPhil) <Ouldridge-DPhil3>`
and  :ref:`(Ouldridge) <Ouldridge3>` for a detailed description of the oxRNA2 force field.

.. note::

   These pair styles have to be used together with the related oxDNA2 bond style
   *oxrna2/fene* for the connectivity of the phosphate backbone (see also documentation of
   :doc:`bond_style oxrna2/fene <bond_oxdna>`). Most of the coefficients
   in the above example have to be kept fixed and cannot be changed without reparameterizing the entire model.
   Exceptions are the first four coefficients after *oxrna2/stk* (seq=seqdep, T=0.1, xi=1.40206 and kappa=2.77 in the above example),
   the first coefficient after *oxrna2/hbond* (seq=seqdep in the above example) and the three coefficients
   after *oxrna2/dh* (T=0.1, rhos=0.5, qeff=1.02455 in the above example). When using a Langevin thermostat
   e.g. through :doc:`fix langevin <fix_langevin>` or :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`
   the temperature coefficients have to be matched to the one used in the fix.

Example input and data files for DNA duplexes can be found in examples/USER/cgdna/examples/oxDNA/ and /oxDNA2/.
A simple python setup tool which creates single straight or helical DNA strands,
DNA duplexes or arrays of DNA duplexes can be found in examples/USER/cgdna/util/.

Please cite :ref:`(Henrich) <Henrich3>` in any publication that uses
this implementation.  The article contains general information
on the model, its implementation and performance as well as the structure of
the data and input file. The preprint version of the article can be found
`here <PDF/USER-CGDNA.pdf>`_.
Please cite also the relevant oxRNA2 publications
:ref:`(Sulc1) <Sulc31>` and :ref:`(Sulc2) <Sulc32>`.

----------


Restrictions
""""""""""""


These pair styles can only be used if LAMMPS was built with the
USER-CGDNA package and the MOLECULE and ASPHERE package.  See the
:doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`bond_style oxrna2/fene <bond_oxdna>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`bond_style oxdna/fene <bond_oxdna>`, :doc:`pair_style oxdna/excv <pair_oxdna>`,
:doc:`bond_style oxdna2/fene <bond_oxdna>`, :doc:`pair_style oxdna2/excv <pair_oxdna2>`,
:doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`

**Default:**

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
