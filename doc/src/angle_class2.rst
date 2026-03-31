.. index:: angle_style class2
.. index:: angle_style class2/kk
.. index:: angle_style class2/omp
.. index:: angle_style class2xe
.. index:: angle_style class2/p6

angle_style class2 command
==========================

Accelerator Variants: *class2/kk*, *class2/omp*

angle_style class2xe command
============================

angle_style class2/p6 command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style class2

Examples
""""""""

.. code-block:: LAMMPS

   angle_style class2
   angle_coeff * 75.0 25.0 0.3 0.002
   angle_coeff 1 bb 10.5872 1.0119 1.5228
   angle_coeff * ba 3.6551 24.895 1.0119 1.5228

Description
"""""""""""

The *class2* angle style uses the potential

.. math::

   E & = E_a + E_{bb} + E_{ba} \\
   E_a & = K_2 (\theta - \theta_0)^2 + K_3 (\theta - \theta_0)^3 + K_4(\theta - \theta_0)^4 \\
   E_{bb} & = M (r_{ij} - r_1) (r_{jk} - r_2) \\
   E_{ba} & = N_1 (r_{ij} - r_1) (\theta - \theta_0) + N_2(r_{jk} - r_2)(\theta - \theta_0)

where :math:`E_a` is the angle term, :math:`E_{bb}` is a bond-bond term,
and :math:`E_{ba}` is a bond-angle term.  :math:`\theta_0` is the
equilibrium angle and :math:`r_1` and :math:`r_2` are the equilibrium
bond lengths.

See :ref:`(Sun) <angle-Sun>` for a description of the COMPASS class2
force field.

Coefficients for the :math:`E_a`, :math:`E_{bb}`, and :math:`E_{ba}`
formulas must be defined for each angle type via the :doc:`angle_coeff
<angle_coeff>` command as in the example above, or in the data file or
restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands.

These are the 4 coefficients for the :math:`E_a` formula:

* :math:`\theta_0` (degrees)
* :math:`K_2` (energy)
* :math:`K_3` (energy)
* :math:`K_4` (energy)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence the various :math:`K` are effectively energy
per radian\^2 or radian\^3 or radian\^4.

For the :math:`E_{bb}` formula, each line in a :doc:`angle_coeff
<angle_coeff>` command in the input script lists 4 coefficients, the
first of which is "bb" to indicate they are BondBond coefficients.  In a
data file, these coefficients should be listed under a "BondBond Coeffs"
heading and you must leave out the "bb", i.e. only list 3 coefficients
after the angle type.

* bb
* :math:`M` (energy/distance\^2)
* :math:`r_1` (distance)
* :math:`r_2` (distance)

For the :math:`E_{ba}` formula, each line in a :doc:`angle_coeff
<angle_coeff>` command in the input script lists 5 coefficients, the
first of which is "ba" to indicate they are BondAngle coefficients.  In
a data file, these coefficients should be listed under a "BondAngle
Coeffs" heading and you must leave out the "ba", i.e. only list 4
coefficients after the angle type.

* ba
* :math:`N_1` (energy/distance)
* :math:`N_2` (energy/distance)
* :math:`r_1` (distance)
* :math:`r_2` (distance)

The :math:`\theta_0` value in the :math:`E_{ba}` formula is not specified,
since it is the same value from the :math:`E_a` formula.

.. note::

   It is important that the order of the I,J,K atoms in each angle
   listed in the Angles section of the data file read by the
   :doc:`read_data <read_data>` command be consistent with the order of
   the :math:`r_1` and :math:`r_2` BondBond and BondAngle coefficients.
   This is because the terms in the formulas for :math:`E_{bb}` and
   :math:`E_{ba}` will use the I,J atoms to compute :math:`r_{ij}` and
   the J,K atoms to compute :math:`r_{jk}`.

----------

.. versionadded:: 30Mar2026

The *class2xe* angle style uses the potential

.. math::

   \begin{aligned}
   E_a = & K_2 (\theta - \theta_0)^2 + K_3 (\theta - \theta_0)^3 + K_4(\theta - \theta_0)^4 \\
   E_{bb} = &  D \left[ 1 - e^{-\alpha (r_{ij} - r_1)} \right] \left[ 1 - e^{-\alpha (r_{jk} - r_2)} \right] \\
   E_{ba} = & D_1 \left[ 1 - e^{-\alpha_1 (r_{ij} - r_1)} \right] \left[\theta - \theta_0\right] + D_2 \left[ 1 - e^{-\alpha_2 (r_{jk} - r_2)} \right] \left[\theta - \theta_0\right]
   \end{aligned}

where :math:`E_a` is the angle term, :math:`E_{bb}` is a bond-bond term
(:math:`D` is the dissociation energy), and :math:`E_{ba}` is a
bond-angle term (:math:`D_1` and :math:`D_2` are the dissociation
energies). :math:`\theta_0` is the equilibrium angle and :math:`r_1` and
:math:`r_2` are the equilibrium bond lengths.

See :ref:`(Kemppainen) <angle-Kemppainen>` for a description of the
ClassII-xe force field and see :doc:`Howto bioFF <Howto_bioFF>` page for
a motivation for the ClassII-xe force field.

.. note::

   The *class2xe* angle style only describes the dissociation of a bond
   stretch. However once a bond is dissociated and stretched beyond the
   processor communication cutoff distance in parallel, the simulation
   will crash with atoms missing errors. This is often after the
   material fractures and thus for post-fracture phenomena the bonded
   interactions need to be removed for proper parallel communication.

   To disconnect the dissociated bond and remove higher order
   interactions (angles, dihedrals, and impropers) the following LAMMPS
   commands can be used with the *class2xe* angle style :doc:`fix
   bond/react <fix_bond_react>` or :doc:`fix bond/break
   <fix_bond_break>`. See the :doc:`Howto bioFF <Howto_bioFF>` page for
   more details.

Coefficients for the :math:`E_a`, :math:`E_{bb}`, and :math:`E_{ba}`
formulas must be defined for each angle type via the :doc:`angle_coeff
<angle_coeff>` command as in the example above, or in the data file or
restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands.

These are the 4 coefficients for the :math:`E_a` formula:

* :math:`\theta_0` (degrees)
* :math:`K_2` (energy)
* :math:`K_3` (energy)
* :math:`K_4` (energy)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence the various :math:`K` are effectively energy
per radian\^2 or radian\^3 or radian\^4.

For the :math:`E_{bb}` formula, each line in a :doc:`angle_coeff
<angle_coeff>` command in the input script lists 5 coefficients, the
first of which is "bb" to indicate they are BondBond coefficients.  In a
data file, these coefficients should be listed under a "BondBond Coeffs"
heading and you must leave out the "bb", i.e. only list 4 coefficients
after the angle type.

* bb
* :math:`D` (energy)
* :math:`\alpha` (inverse distance)
* :math:`r_1` (distance)
* :math:`r_2` (distance)

For the :math:`E_{ba}` formula, each line in a :doc:`angle_coeff
<angle_coeff>` command in the input script lists 7 coefficients, the
first of which is "ba" to indicate they are BondAngle coefficients.  In
a data file, these coefficients should be listed under a "BondAngle
Coeffs" heading and you must leave out the "ba", i.e. only list 6
coefficients after the angle type.

* ba
* :math:`D_1` (energy)
* :math:`D_2` (energy)
* :math:`\alpha_1` (inverse distance)
* :math:`\alpha_2` (inverse distance)
* :math:`r_1` (distance)
* :math:`r_2` (distance)

The :math:`\theta_0` value in the :math:`E_{ba}` formula is not specified,
since it is the same value from the :math:`E_a` formula.

----------

The *class2/p6* angle style uses the *class2* potential expanded to sixth order:

.. math::

   E_{a} = K_2\left(\theta - \theta_0\right)^2 + K_3\left(\theta - \theta_0\right)^3 + K_4\left(\theta - \theta_0\right)^4 + K_5\left(\theta - \theta_0\right)^5 + K_6\left(\theta - \theta_0\right)^6

In this expanded term 6 coefficients for the :math:`E_a` formula need to
be set:

* :math:`\theta_0` (degrees)
* :math:`K_2` (energy)
* :math:`K_3` (energy)
* :math:`K_4` (energy)
* :math:`K_5` (energy)
* :math:`K_6` (energy)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence the various :math:`K` are effectively energy
per radian\^2 or radian\^3 or radian\^4 or radian\^5 or radian\^6.

The bond-bond and bond-angle terms remain unchanged.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

The *class2* and *class2xe* angle styles can only be used if LAMMPS was
built with the CLASS2 package.  For the *class2/p6* angle style LAMMPS
needs to be built with the MOFFF package.  See the :doc:`Build package
<Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none

----------

.. _angle-Sun:

**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).

.. _angle-Kemppainen:

**(Kemppainen)** Kemppainen, npj Computational Materials 11, 341 (2025).
