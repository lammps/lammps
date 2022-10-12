.. index:: bond_style harmonic
.. index:: bond_style harmonic/intel
.. index:: bond_style harmonic/kk
.. index:: bond_style harmonic/omp

bond_style harmonic command
===========================

Accelerator Variants: *harmonic/intel*, *harmonic/kk*, *harmonic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style harmonic

Examples
""""""""

.. code-block:: LAMMPS

   bond_style harmonic
   bond_coeff 5 80.0 1.2

Description
"""""""""""

The *harmonic* bond style uses the potential

.. math::

   E = K (r - r_0)^2

where :math:`r_0` is the equilibrium bond distance.  Note that the usual 1/2
factor is included in :math:`K`.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/distance\^2)
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
