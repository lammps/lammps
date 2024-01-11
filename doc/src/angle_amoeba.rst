.. index:: angle_style amoeba

angle_style amoeba command
==========================


Syntax
""""""

.. code-block:: LAMMPS

   angle_style amoeba

Examples
""""""""

.. code-block:: LAMMPS

   angle_style amoeba
   angle_coeff * 75.0 -25.0 1.0 0.3 0.02 0.003
   angle_coeff * ba 3.6551 24.895 1.0119 1.5228
   angle_coeff * ub -7.6 1.5537

Description
"""""""""""

The *amoeba* angle style uses the potential

.. math::

   E & = E_a + E_{ba} + E_{ub} \\
   E_a & = K_2\left(\theta - \theta_0\right)^2 + K_3\left(\theta - \theta_0\right)^3 + K_4\left(\theta - \theta_0\right)^4 + K_5\left(\theta - \theta_0\right)^5 + K_6\left(\theta - \theta_0\right)^6 \\
   E_{ba} & = N_1 (r_{ij} - r_1) (\theta - \theta_0) + N_2(r_{jk} - r_2)(\theta - \theta_0) \\
   E_{UB} & = K_{ub} (r_{ik} - r_{ub})^2

where :math:`E_a` is the angle term, :math:`E_{ba}` is a bond-angle
term, :math:`E_{UB}` is a Urey-Bradley bond term, :math:`\theta_0` is
the equilibrium angle, :math:`r_1` and :math:`r_2` are the equilibrium
bond lengths, and :math:`r_{ub}` is the equilibrium Urey-Bradley bond
length.

These formulas match how the Tinker MD code performs its angle
calculations for the AMOEBA and HIPPO force fields.  See the
:doc:`Howto amoeba <Howto_amoeba>` page for more information about
the implementation of AMOEBA and HIPPO in LAMMPS.

Note that the :math:`E_a` and :math:`E_{ba}` formulas are identical to
those used for the :doc:`angle_style class2/p6 <angle_class2>`
command, however there is no bond-bond cross term formula for
:math:`E_{bb}`.  Additionally, there is a :math:`E_{UB}` term for a
Urey-Bradley bond.  It is effectively a harmonic bond between the I
and K atoms of angle IJK, even though that bond is not enumerated in
the "Bonds" section of the data file.

There are also two ways that Tinker computes the angle :math:`\theta`
in the :math:`E_a` formula.  The first is the standard way of treating
IJK as an "in-plane" angle.  The second is an "out-of-plane" method
which Tinker may use if the center atom J in the angle is bonded to
one additional atom in addition to I and K.  In this case, all 4 atoms
are used to compute the :math:`E_a` formula, resulting in forces on
all 4 atoms.  In the Tinker PRM file, these 2 options are denoted by
*angle* versus *anglep* entries in the "Angle Bending Parameters"
section of the PRM force field file.  The *pflag* coefficient
described below selects between the 2 options.

----------


Coefficients for the :math:`E_a`, :math:`E_{bb}`, and :math:`E_{ub}`
formulas must be defined for each angle type via the :doc:`angle_coeff
<angle_coeff>` command as in the example above, or in the data file or
restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands.

These are the 8 coefficients for the :math:`E_a` formula:

* pflag = 0 or 1
* ubflag = 0 or 1
* :math:`\theta_0` (degrees)
* :math:`K_2` (energy)
* :math:`K_3` (energy)
* :math:`K_4` (energy)
* :math:`K_5` (energy)
* :math:`K_6` (energy)

A pflag value of 0 vs 1 selects between the "in-plane" and
"out-of-plane" options described above.  Ubflag is 1 if there is a
Urey-Bradley term associated with this angle type, else it is 0.
:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence the various :math:`K` values are effectively
energy per radian\^2 or radian\^3 or radian\^4 or radian\^5 or
radian\^6.

For the :math:`E_{ba}` formula, each line in a :doc:`angle_coeff
<angle_coeff>` command in the input script lists 5 coefficients, the
first of which is "ba" to indicate they are BondAngle coefficients.
In a data file, these coefficients should be listed under a "BondAngle
Coeffs" heading and you must leave out the "ba", i.e. only list 4
coefficients after the angle type.

* ba
* :math:`N_1` (energy/distance\^2)
* :math:`N_2` (energy/distance\^2)
* :math:`r_1` (distance)
* :math:`r_2` (distance)

The :math:`\theta_0` value in the :math:`E_{ba}` formula is not specified,
since it is the same value from the :math:`E_a` formula.

For the :math:`E_{ub}` formula, each line in a :doc:`angle_coeff
<angle_coeff>` command in the input script lists 3 coefficients, the
first of which is "ub" to indicate they are UreyBradley coefficients.
In a data file, these coefficients should be listed under a
"UreyBradley Coeffs" heading and you must leave out the "ub",
i.e. only list 2 coefficients after the angle type.

* ub
* :math:`K_{ub}` (energy/distance\^2)
* :math:`r_{ub}` (distance)

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the AMOEBA
package.  See the :doc:`Build package <Build_package>` doc page for
more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none
