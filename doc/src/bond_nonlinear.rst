.. index:: bond_style nonlinear
.. index:: bond_style nonlinear/omp

bond_style nonlinear command
============================

Accelerator Variants: *nonlinear/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style nonlinear

Examples
""""""""

.. code-block:: LAMMPS

   bond_style nonlinear
   bond_coeff 2 100.0 1.1 1.4

Description
"""""""""""

The *nonlinear* bond style uses the potential

.. math::

   E = \frac{\epsilon (r - r_0)^2}{ [ \lambda^2 - (r - r_0)^2 ]}

to define an anharmonic spring :ref:`(Rector) <Rector>` of equilibrium
length :math:`r_0` and maximum extension lamda.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`\epsilon` (energy)
* :math:`r_0` (distance)
* :math:`\lambda` (distance)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the EXTRA-MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none

----------

.. _Rector:

**(Rector)** Rector, Van Swol, Henderson, Molecular Physics, 82, 1009 (1994).
