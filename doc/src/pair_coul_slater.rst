.. index:: pair_style coul/slater

pair_style coul/slater/cut command
============================

pair_style coul/slater/long command
================================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style coul/slater/cut lamda cutoff
   pair_style coul/slater/long lamda cutoff

lamda = decay length of the charge (distance units)
cutoff = cutoff (distance units)

Examples
""""""""


.. code-block:: LAMMPS

   pair_style coul/slater/cut 1.0 3.5
   pair_coeff * *
   pair_coeff 2 2 2.5

   pair_style coul/slater/long 1.0 12.0
   pair_coeff * *
   pair_coeff 1 1 5.0

Description
"""""""""""

Styles *coul/slater* compute electrostatic interactions in mesoscopic models
which employ potentials without explicit excluded-volume interactions. 
The goal is to prevent artificial ionic pair formation by including a charge 
distribution in the Coulomb potential, following the formulation of 
:ref:`(Melchor) <Melchor>`: 

.. math::

   E  = & \frac{Cq_iq_j}{\epsilon r} \left( 1- \left( 1 + \frac{r_{ij}}{\lamda} exp\left( -2r_{ij}/\lamda \right) \right) \right)                       \qquad r < r_c 


where :math:`r_c` is the cutoff distance and :math:`\lamda` is the decay length of the charge.
C is the same Coulomb conversion factor as in the pair\_styles coul/cut and coul/long. In this way the Coulomb
interaction between ions is corrected at small distances r. 
For the *coul/slater/cut* style, the potential energy for distances larger than the cutoff is zero, 
while for the *coul/slater/long*, the long-range interactions are computed either by the Ewald or the PPPM technique.

Phenomena that can be captured at a mesoscopic level using this type of electrostatic 
interactions include the formation of polyelectrolyte-surfactant aggregates, 
charge stabilization of colloidal suspensions, and the formation of
complexes driven by charged species in biological systems. :ref:`(Vaiwala) <Vaiwala>`.

The cutoff distance is optional. If it is not used,
the default global value specified in the pair_style command is used.
For each pair of atom types, a specific cutoff distance can be defined via the :doc:`pair_coeff <pair_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`\r_c` (distance units)

The global decay length of the charge (:math:`\lambda`) specified in the pair\_style command is used for all pairs.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the cutoff distance for the
*coul/slater* styles can be mixed.  The default mix value is *geometric*\ .
See the "pair\_modify" command for details.

The :doc:`pair_modify <pair_modify>` shift and table options are not relevant
for these pair styles.

These pair styles do not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

These pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.

Restrictions
""""""""""""

The  *coul/slater/long* style requires the long-range solvers included in the KSPACE package. 

These styles are part of the "USER-MISC" package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style, hybrid/overlay <pair_hybrid>`, :doc:`kspace_style <kspace_style>`

**Default:** none

----------


.. _Melchor:

**(Melchor)** González-Melchor, Mayoral, Velázquez, and Alejandre, J Chem Phys, 125, 224107 (2006).

.. _Vaiwala:

**(Vaiwala)** Vaiwala, Jadhav, and Thaokar, J Chem Phys, 146, 124904 (2017).


