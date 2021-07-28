.. index:: bond_style gaussian

bond_style gaussian command
================================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style gaussian

Examples
""""""""

.. code-block:: LAMMPS

   bond_style gaussian
   bond_coeff 1 300.0 2 0.0128 0.375 3.37 0.0730 0.148 3.63

Description
"""""""""""

The *gaussian* bond style uses the potential:

.. math::

   E = -k_B T ln\left(\sum_{i=1}^{n} \frac{A_i}{w_i \sqrt{\pi/2}} exp\left( \frac{-(r-r_{i})^2}{w_i^2})\right) \right)

This analytical form is a suitable potential for obtaining
mesoscale effective force fields which can reproduce target atomistic distributions :ref:`(Milano) <Milano0>`

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* T temperature at which the potential was derived
* :math:`n` (integer >=1)
* :math:`A_1` (-)
* :math:`w_1` (-)
* :math:`r_1` (length)
* ...
* :math:`A_n` (-)
* :math:`w_n` (-)
* :math:`r_n` (length)


Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`

Default
"""""""

none

----------

.. _Milano0:

**(Milano)** G. Milano, S. Goudeau, F. Mueller-Plathe, J. Polym. Sci. B Polym. Phys. 43, 871 (2005).
