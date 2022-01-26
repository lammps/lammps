.. index:: bond_style morse
.. index:: bond_style morse/omp

bond_style morse command
========================

Accelerator Variants: *morse/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style morse/lattice
   bond_coeff N D alpha b0 kappa nx ny nz


   * N = bond type (see asterisk form in :doc:`bond_coeff <bond_coeff>`)
   * D = well depth
   * alpha = stiffness
   * b0 = natural length (length units)
   * kappa = transverse strength (dimensionless)
   * nx, ny, nz = bond vector direction (will be normalised)

Examples
""""""""

.. code-block:: LAMMPS

   bond_style morse/lattice
   # for e.g. fcc - 6 bonding directions
   bond_coeff 1 1.0 2.0 1.2 0.0 1 1 0
   bond_coeff 2 1.0 2.0 1.2 0.0 1 -1 0
   bond_coeff 3 1.0 2.0 1.2 0.0 1 0 1
   bond_coeff 4 1.0 2.0 1.2 0.0 1 0 -1
   bond_coeff 5 1.0 2.0 1.2 0.0 0 1 1
   bond_coeff 6 1.0 2.0 1.2 0.0 0 1 -1

Description
"""""""""""

The *morse/lattice* bond style uses the potential

.. math::

   E = D\left[ 1 - e^{-\alpha (b_l - b_0)} \right]^2 + D\left[ 1 - e^{\alpha (b_l - b_0)} \right]^2 + \kappa D \alpha^2\left[b^2 - b_l^2\right]

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

This bond style can only be used if LAMMPS was built with the EXTRA-MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none
