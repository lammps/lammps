.. index:: bond_style morse
.. index:: bond_style morse/omp

bond_style morse command
========================

Accelerator Variants: *morse/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style morse

Examples
""""""""

.. code-block:: LAMMPS

   bond_style morse
   bond_coeff 5 1.0 2.0 1.2

Description
"""""""""""

The *morse* bond style uses the potential

.. math::

   E = D \left[ 1 - e^{-\alpha (r - r_0)} \right]^2

where :math:`r_0` is the equilibrium bond distance, :math:`\alpha` is a stiffness
parameter, and :math:`D` determines the depth of the potential well.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`D` (energy)
* :math:`\alpha` (inverse distance)
* :math:`r_0` (distance)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none
