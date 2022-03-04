.. index:: bond_style harmonic/shift/cut
.. index:: bond_style harmonic/shift/cut/omp

bond_style harmonic/shift/cut command
=====================================

Accelerator Variants: *harmonic/shift/cut/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style harmonic/shift/cut

Examples
""""""""

.. code-block:: LAMMPS

   bond_style harmonic/shift/cut
   bond_coeff 5 10.0 0.5 1.0

Description
"""""""""""

The *harmonic/shift/cut* bond style is a shifted harmonic bond that
uses the potential

.. math::

   E = \frac{U_{\text{min}}}{(r_0-r_c)^2} \left[ (r-r_0)^2-(r_c-r_0)^2 \right]

where :math:`r_0` is the equilibrium bond distance, and rc the critical distance.
The bond potential is zero for distances :math:`r > r_c`. The potential is :math:`-U_{\text{min}}`
at :math:`r_0` and zero at :math:`r_c`. The spring constant is :math:`k = U_{\text{min}} / [ 2 (r_0-r_c)^2]`.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`U_{\text{min}}` (energy)
* :math:`r_0` (distance)
* :math:`r_c` (distance)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`,
:doc:`bond_harmonic <bond_harmonic>`,
:doc:`bond_style harmonic/shift <bond_harmonic_shift>`

Default
"""""""

none
