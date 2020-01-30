.. index:: bond_style special

bond_style special command
=================================

Syntax
""""""


.. code-block:: LAMMPS

   bond_style special

Examples
""""""""


.. code-block:: LAMMPS

   bond_style special
   bond_coeff 0.5 0.5

Description
"""""""""""

The *special* bond style can be used to impose weighted Lennard Jones and/or
Coloumbic interactions on any two particles in the system. It can be used for
cases that cannot be handled in :doc:`special_bonds <special_bonds>`, such as
1-5 interactions or hyrbrid simulations that require the use of different
weighting factors for different molecules.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* `w_lj` weight (0.0 to 1.0) on pairwise Lennard-Jones interactions

* `w_coul` weight (0.0 to 1.0) on pairwise Coulombic interactions

This command will typically be used in conjunction with 

----------



Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`special_bonds <special_bonds>`,

**Default:** none
