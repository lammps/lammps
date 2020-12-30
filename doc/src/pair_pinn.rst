.. index:: pair_style pinn

pair_style pinn command
======================

Currently no Accelerator Variants

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style

* style = *pinn*

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pinn
   pair_coeff * * ../potentials/Al_2020.pinn Al Al

Description
"""""""""""

Style *eam* computes pairwise interactions of atoms using Physically 
Informed Neural Network (PINN) potentials. Details of the model can be 
found in :ref:`(Purja) <Purja>`. It supports multicomponent systems. The 
potential files in the *potentials* directory of the LAMMPS distribution 
have an ".pinn" suffix.

.. note::

   Note that unlike for other potentials, species and cutoffs for PINN
   potentials are not set in the pair_style or pair_coeff command; they
   are specified in the PINN potential files themselves.  Likewise, the
   PINN potential files list atomic masses; thus you do not need to use
   the :doc:`mass <mass>` command to specify them.

----------

As an example, the potentials/Al_2020.pinn file is a single element PINN 
file which has values for aluminum. If your LAMMPS simulation has 3 
atoms types, you would use the following pair_coeff
command:

.. code-block:: LAMMPS

   pair_coeff * * Al_2020.pinn Al Al Al
   
The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
Three Al arguments map LAMMPS atom types 1,2,3 to the Al
element in the *setfl* file.  Note that there is no requirement that 
your simulation use all the elements specified by the *setfl* file for 
a multicomponent system.

----------

Restrictions
""""""""""""

This style is part of the USER-PINN package.  It is only enabled if 
LAMMPS was built with that package.  See the :doc:`Build package 
<Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Purja:

**(Purja Pun)** G. P. Purja Pun, V. Yamakov , J. Hickman , E. H. Glaessgen 
and Y. Mishin, Phys Rev Mat 4, 113807 (2020)
 
