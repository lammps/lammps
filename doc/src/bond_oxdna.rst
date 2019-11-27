.. index:: bond_style oxdna/fene

bond_style oxdna/fene command
=============================

bond_style oxdna2/fene command
==============================

Syntax
""""""


.. code-block:: LAMMPS

   bond_style oxdna/fene

   bond_style oxdna2/fene

Examples
""""""""


.. code-block:: LAMMPS

   bond_style oxdna/fene
   bond_coeff * 2.0 0.25 0.7525

   bond_style oxdna2/fene
   bond_coeff * 2.0 0.25 0.7564

Description
"""""""""""

The *oxdna/fene* and *oxdna2/fene* bond styles use the potential

.. math::

   E = - \frac{\epsilon}{2} \ln \left[ 1 - \left(\frac{r-r_0}{\Delta}\right)^2\right]


to define a modified finite extensible nonlinear elastic (FENE)
potential :ref:`(Ouldridge) <oxdna_fene>` to model the connectivity of the
phosphate backbone in the oxDNA force field for coarse-grained
modelling of DNA.

The following coefficients must be defined for the bond type via the
:doc:`bond\_coeff <bond_coeff>` command as given in the above example, or
in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* :math:`\epsilon` (energy)
* :math:`\Delta` (distance)
* :math:`r_0` (distance)

.. note::

   The oxDNA bond style has to be used together with the
   corresponding oxDNA pair styles for excluded volume interaction
   *oxdna/excv*\ , stacking *oxdna/stk*\ , cross-stacking *oxdna/xstk* and
   coaxial stacking interaction *oxdna/coaxstk* as well as
   hydrogen-bonding interaction *oxdna/hbond* (see also documentation of
   :doc:`pair\_style oxdna/excv <pair_oxdna>`). For the oxDNA2
   :ref:`(Snodin) <oxdna2>` bond style the analogous pair styles and an
   additional Debye-Hueckel pair style *oxdna2/dh* have to be defined.
   The coefficients in the above example have to be kept fixed and cannot
   be changed without reparameterizing the entire model.

Example input and data files for DNA duplexes can be found in
examples/USER/cgdna/examples/oxDNA/ and /oxDNA2/.  A simple python
setup tool which creates single straight or helical DNA strands, DNA
duplexes or arrays of DNA duplexes can be found in
examples/USER/cgdna/util/.

Please cite :ref:`(Henrich) <Henrich2>` and the relevant oxDNA articles in
any publication that uses this implementation.  The article contains
more information on the model, the structure of the input file, the
setup tool and the performance of the LAMMPS-implementation of oxDNA.
The preprint version of the article can be found
`here <PDF/USER-CGDNA.pdf>`_.


----------


Restrictions
""""""""""""


This bond style can only be used if LAMMPS was built with the
USER-CGDNA package and the MOLECULE and ASPHERE package.  See the
:doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_style oxdna/excv <pair_oxdna>`, :doc:`pair\_style oxdna2/excv <pair_oxdna2>`, :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`,
:doc:`bond\_coeff <bond_coeff>`

**Default:** none


----------


.. _Henrich2:



**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk,
T. E. Ouldridge, Eur. Phys. J. E 41, 57 (2018).

.. _oxdna\_fene:



**(Ouldridge)** T.E. Ouldridge, A.A. Louis, J.P.K. Doye,
J. Chem. Phys. 134, 085101 (2011).

.. _oxdna2:



**(Snodin)** B.E. Snodin, F. Randisi, M. Mosayebi, et al.,
J. Chem. Phys. 142, 234901 (2015).
