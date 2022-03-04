.. index:: bond_style class2
.. index:: bond_style class2/omp
.. index:: bond_style class2/kk

bond_style class2 command
=========================

Accelerator Variants: *class2/omp*, *class2/kk*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style class2

Examples
""""""""

.. code-block:: LAMMPS

   bond_style class2
   bond_coeff 1 1.0 100.0 80.0 80.0

Description
"""""""""""

The *class2* bond style uses the potential

.. math::

   E = K_2 (r - r_0)^2 + K_3 (r - r_0)^3 + K_4 (r - r_0)^4

where :math:`r_0` is the equilibrium bond distance.

See :ref:`(Sun) <bond-Sun>` for a description of the COMPASS class2 force field.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`r_0` (distance)
* :math:`K_2` (energy/distance\^2)
* :math:`K_3` (energy/distance\^3)
* :math:`K_4` (energy/distance\^4)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the CLASS2
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none

----------

.. _bond-Sun:

**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).
