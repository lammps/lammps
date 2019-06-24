.. index:: angle\_style mm3

angle\_style mm3 command
========================

Syntax
""""""


.. parsed-literal::

   angle_style mm3

Examples
""""""""


.. parsed-literal::

   angle_style mm3
   angle_coeff 1 100.0 107.0

Description
"""""""""""

The *mm3* angle style uses the potential that is anharmonic in the angle
as defined in :ref:`(Allinger) <mm3-allinger1989>`

.. math::

   E = K (\theta - \theta_0)^2 \left[ 1 - 0.014(\theta - \theta_0) + 5.6(10)^{-5} (\theta - \theta_0)^2 - 7.0(10)^{-7} (\theta - \theta_0)^3 + 9(10)^{-10} (\theta - \theta_0)^4 \right]


where theta0 is the equilibrium value of the angle, and K is a
prefactor. The anharmonic prefactors have units deg\^(-n), for example
-0.014 deg\^(-1), 5.6(10)\^(-5) deg\^(-2), ...

The following coefficients must be defined for each angle type via the
:doc:`angle\_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* K (energy/radian\^2)
* theta0 (degrees)

Theta0 is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of K are in energy/radian\^2.

Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
USER\_YAFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle\_coeff <angle_coeff>`

**Default:** none


----------



.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
