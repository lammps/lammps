.. index:: bond_style quartic/exp
.. index:: bond_style quartic/exp/kk

bond_style quartic/exp command
===============================

Accelerator Variants: *quartic/exp/kk*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style quartic/exp

Examples
""""""""

.. code-block:: LAMMPS

   bond_style quartic/exp
   bond_coeff 1 1.54 200.0 -100.0 50.0 0.0 1.0
   bond_coeff 2 1.0 500.0 0.0 -200.0 1000.0 0.3

Description
"""""""""""

.. versionadded:: 4Jul2026

The *quartic/exp* bond style uses the potential

.. math::

   E = k_2 (r - r_0)^2 + k_3 (r - r_0)^3 + k_4 (r - r_0)^4
       + A \exp\left(-\frac{r}{B}\right)

where :math:`r_0` is the equilibrium bond distance.  The first three
terms form a quartic polynomial in the bond displacement
:math:`(r - r_0)` that can model anharmonic stretching.  The optional
exponential term (:math:`A \exp(-r/B)`) adds a repulsive or attractive
wall at short range.  Setting :math:`A = 0` disables the exponential
term entirely.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`r_0` (distance)
* :math:`k_2` (energy/distance\^2)
* :math:`k_3` (energy/distance\^3)
* :math:`k_4` (energy/distance\^4)
* :math:`A` (energy)
* :math:`B` (distance, must be non-zero if :math:`A \neq 0`)

Note that the :math:`k_2` coefficient plays the role of the spring
constant.  The :math:`k_3` and :math:`k_4` coefficients introduce
anharmonicity.  Setting :math:`k_3 = k_4 = 0` and :math:`A = 0`
reduces this style to :doc:`bond_style harmonic <bond_harmonic>` with
:math:`K = k_2` (the factor of 1/2 is *not* included in :math:`k_2`
here).

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>`
page for more info.

The parameter :math:`B` must be non-zero whenever :math:`A \neq 0`.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`,
:doc:`bond_style harmonic <bond_harmonic>`,
:doc:`bond_style nonlinear <bond_nonlinear>`

Default
"""""""

none
