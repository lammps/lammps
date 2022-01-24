.. index:: bond_style gromos
.. index:: bond_style gromos/omp

bond_style gromos command
=========================

Accelerator Variants: *gromos/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style gromos

Examples
""""""""

.. code-block:: LAMMPS

   bond_style gromos
   bond_coeff 5 80.0 1.2

Description
"""""""""""

The *gromos* bond style uses the potential

.. math::

   E = K (r^2 - r_0^2)^2

where :math:`r_0` is the equilibrium bond distance.  Note that the usual 1/4
factor is included in :math:`K`.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/distance\^4)
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
