.. index:: bond_style oxdna/fene

bond_style oxdna/fene command
=============================

bond_style oxdna2/fene command
==============================

bond_style oxrna2/fene command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style oxdna/fene

   bond_style oxdna2/fene

   bond_style oxrna2/fene

Examples
""""""""

.. code-block:: LAMMPS

   bond_style oxdna/fene
   bond_coeff * 2.0 0.25 0.7525

   bond_style oxdna2/fene
   bond_coeff * 2.0 0.25 0.7564

   bond_style oxrna2/fene
   bond_coeff * 2.0 0.25 0.76107

Description
"""""""""""

The *oxdna/fene* , *oxdna2/fene* and *oxrna2/fene* bond styles use the potential

.. math::

   E = - \frac{\epsilon}{2} \ln \left[ 1 - \left(\frac{r-r_0}{\Delta}\right)^2\right]

to define a modified finite extensible nonlinear elastic (FENE)
potential :ref:`(Ouldridge) <Ouldridge0>` to model the connectivity of the
phosphate backbone in the oxDNA/oxRNA force field for coarse-grained
modelling of DNA/RNA.

The following coefficients must be defined for the bond type via the
:doc:`bond_coeff <bond_coeff>` command as given in the above example, or
in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`\epsilon` (energy)
* :math:`\Delta` (distance)
* :math:`r_0` (distance)

.. note::

   The oxDNA bond style has to be used together with the
   corresponding oxDNA pair styles for excluded volume interaction
   *oxdna/excv* , stacking *oxdna/stk* , cross-stacking *oxdna/xstk* and
   coaxial stacking interaction *oxdna/coaxstk* as well as
   hydrogen-bonding interaction *oxdna/hbond* (see also documentation of
   :doc:`pair_style oxdna/excv <pair_oxdna>`). For the oxDNA2
   :ref:`(Snodin) <Snodin0>` bond style the analogous pair styles
   *oxdna2/excv* , *oxdna2/stk* , *oxdna2/xstk* , *oxdna2/coaxstk* ,
   *oxdna2/hbond* and an additional Debye-Hueckel pair style
   *oxdna2/dh* have to be defined. The same applies to the oxRNA2
   :ref:`(Sulc1) <Sulc01>` styles.
   The coefficients in the above example have to be kept fixed and cannot
   be changed without reparameterizing the entire model.

Example input and data files for DNA and RNA duplexes can be found in
examples/USER/cgdna/examples/oxDNA/ , /oxDNA2/ and /oxRNA2/.  A simple python
setup tool which creates single straight or helical DNA strands, DNA/RNA
duplexes or arrays of DNA/RNA duplexes can be found in
examples/USER/cgdna/util/.

Please cite :ref:`(Henrich) <Henrich0>` in any publication that uses
this implementation.  The article contains general information
on the model, its implementation and performance as well as the structure of
the data and input file. The preprint version of the article can be found
`here <PDF/USER-CGDNA.pdf>`_.
Please cite also the relevant oxDNA/oxRNA publications. These are
:ref:`(Ouldridge) <Ouldridge0>` and
:ref:`(Ouldridge-DPhil) <Ouldridge-DPhil0>` for oxDNA,
:ref:`(Snodin) <Snodin0>` for oxDNA2,
:ref:`(Sulc1) <Sulc01>` for oxRNA2
and for sequence-specific hydrogen-bonding and stacking interactions
:ref:`(Sulc2) <Sulc02>`.

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
USER-CGDNA package and the MOLECULE and ASPHERE package.  See the
:doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_style oxdna/excv <pair_oxdna>`, :doc:`pair_style oxdna2/excv <pair_oxdna2>`, :doc:`pair_style oxrna2/excv <pair_oxrna2>`,
:doc:`bond_coeff <bond_coeff>`, :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`

**Default:**

none

----------

.. _Henrich0:

**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E. Ouldridge, Eur. Phys. J. E 41, 57 (2018).

.. _Ouldridge-DPhil0:

**(Ouldridge-DPhil)** T.E. Ouldridge, Coarse-grained modelling of DNA and DNA self-assembly, DPhil. University of Oxford (2011).

.. _Ouldridge0:

**(Ouldridge)** T.E. Ouldridge, A.A. Louis, J.P.K. Doye, J. Chem. Phys. 134, 085101 (2011).

.. _Snodin0:

**(Snodin)** B.E. Snodin, F. Randisi, M. Mosayebi, et al., J. Chem. Phys. 142, 234901 (2015).

.. _Sulc01:

**(Sulc1)** P. Sulc, F. Romano, T. E. Ouldridge, et al., J. Chem. Phys. 140, 235102 (2014).

.. _Sulc02:

**(Sulc2)** P. Sulc, F. Romano, T.E. Ouldridge, L. Rovigatti, J.P.K. Doye, A.A. Louis, J. Chem. Phys. 137, 135101 (2012).
