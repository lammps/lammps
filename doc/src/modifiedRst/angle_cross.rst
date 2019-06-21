.. index:: angle\_style cross

angle\_style cross command
==========================

Syntax
""""""


.. parsed-literal::

   angle_style cross

Examples
""""""""


.. parsed-literal::

   angle_style cross
   angle_coeff 1 200.0 100.0 100.0 1.25 1.25 107.0

Description
"""""""""""

The *cross* angle style uses a potential that couples the bond stretches of
a bend with the angle stretch of that bend:

.. math::

  E = K_{SS} \left(r_{12}-r_{12,0}\right)\left(r_{32}-r_{32,0}\right) + K_{BS0}\left(r_{12}-r_{12,0}\right)\left(\theta-\theta_0\right) + K_{BS1}\left(r_{32}-r_{32,0}\right)\left(\theta-\theta_0\right)


where r12,0 is the rest value of the bond length between atom 1 and 2,
r32,0 is the rest value of the bond length between atom 2 and 2,
and theta0 is the rest value of the angle. KSS is the force constant of
the bond stretch-bond stretch term and KBS0 and KBS1 are the force constants
of the bond stretch-angle stretch terms.

The following coefficients must be defined for each angle type via the
:doc:`angle\_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* KSS (energy/distance\^2)
* KBS0 (energy/distance/rad)
* KBS1 (energy/distance/rad)
* r12,0 (distance)
* r32,0 (distance)
* theta0 (degrees)

Theta0 is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of KBS0 and KBS1 are in energy/distance/radian.

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
