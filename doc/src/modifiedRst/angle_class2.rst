.. index:: angle\_style class2

angle\_style class2 command
===========================

angle\_style class2/kk command
==============================

angle\_style class2/omp command
===============================

angle\_style class2/p6 command
==============================

Syntax
""""""


.. parsed-literal::

   angle_style class2

Examples
""""""""


.. parsed-literal::

   angle_style class2
   angle_coeff \* 75.0
   angle_coeff 1 bb 10.5872 1.0119 1.5228
   angle_coeff \* ba 3.6551 24.895 1.0119 1.5228

Description
"""""""""""

The *class2* angle style uses the potential

.. math::

   E & = & E_a + E_{bb} + E_{ba} \\
   E_a & = & K_2 (\theta - \theta_0)^2 + K_3 (\theta - \theta_0)^3 + K_4 (\theta - \theta_0)^4 \\
   E_{bb} & = & M (r_{ij} - r_1) (r_{jk} - r_2) \\
   E_{ba} & = & N_1 (r_{ij} - r_1) (\theta - \theta_0) + N_2 (r_{jk} - r_2) (\theta - \theta_0)


where Ea is the angle term, Ebb is a bond-bond term, and Eba is a
bond-angle term.  Theta0 is the equilibrium angle and r1 and r2 are
the equilibrium bond lengths.

See :ref:`(Sun) <angle-Sun>` for a description of the COMPASS class2 force field.

Coefficients for the Ea, Ebb, and Eba formulas must be defined for
each angle type via the :doc:`angle\_coeff <angle_coeff>` command as in
the example above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands.

These are the 4 coefficients for the Ea formula:

* theta0 (degrees)
* K2 (energy/radian\^2)
* K3 (energy/radian\^3)
* K4 (energy/radian\^4)

Theta0 is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of the various K are in per-radian.

For the Ebb formula, each line in a :doc:`angle\_coeff <angle_coeff>`
command in the input script lists 4 coefficients, the first of which
is "bb" to indicate they are BondBond coefficients.  In a data file,
these coefficients should be listed under a "BondBond Coeffs" heading
and you must leave out the "bb", i.e. only list 3 coefficients after
the angle type.

* bb
* M (energy/distance\^2)
* r1 (distance)
* r2 (distance)

For the Eba formula, each line in a :doc:`angle\_coeff <angle_coeff>`
command in the input script lists 5 coefficients, the first of which
is "ba" to indicate they are BondAngle coefficients.  In a data file,
these coefficients should be listed under a "BondAngle Coeffs" heading
and you must leave out the "ba", i.e. only list 4 coefficients after
the angle type.

* ba
* N1 (energy/distance\^2)
* N2 (energy/distance\^2)
* r1 (distance)
* r2 (distance)

The theta0 value in the Eba formula is not specified, since it is the
same value from the Ea formula.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


The *class2/p6* angle style uses the *class2* potential expanded to sixth order:

.. image:: Eqs/angle_class2_p6.jpg
   :align: center

In this expanded term 6 coefficients for the Ea formula need to be set:

* theta0 (degrees)
* K2 (energy/radian\^2)
* K3 (energy/radian\^3)
* K4 (energy/radian\^4)
* K5 (energy/radian\^5)
* K6 (energy/radian\^6)

The bond-bond and bond-angle terms remain unchanged.


----------


Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the CLASS2
package.  For the *class2/p6* style LAMMPS needs to be built with the
USER-MOFFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle\_coeff <angle_coeff>`

**Default:** none


----------


.. _angle-Sun:



**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
