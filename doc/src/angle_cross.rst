.. index:: angle_style cross

angle_style cross command
==========================

Syntax
""""""


.. code-block:: LAMMPS

   angle_style cross

Examples
""""""""


.. code-block:: LAMMPS

   angle_style cross
   angle_coeff 1 200.0 100.0 100.0 1.25 1.25 107.0

Description
"""""""""""

The *cross* angle style uses a potential that couples the bond stretches of
a bend with the angle stretch of that bend:

.. math::

   E = K_{SS} \left(r_{12}-r_{12,0}\right)\left(r_{32}-r_{32,0}\right) + K_{BS0}\left(r_{12}-r_{12,0}\right)\left(\theta-\theta_0\right) + K_{BS1}\left(r_{32}-r_{32,0}\right)\left(\theta-\theta_0\right)

where :math:`r_{12,0}` is the rest value of the bond length between atom 1 and 2,
:math:`r_{32,0}` is the rest value of the bond length between atom 3 and 2,
and :math:`\theta_0` is the rest value of the angle. :math:`K_{SS}` is the force constant of
the bond stretch-bond stretch term and :math:`K_{BS0}` and :math:`K_{BS1}` are the force constants
of the bond stretch-angle stretch terms.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K_{SS}` (energy/distance\^2)
* :math:`K_{BS0}` (energy/distance/rad)
* :math:`K_{BS1}` (energy/distance/rad)
* :math:`r_{12,0}` (distance)
* :math:`r_{32,0}` (distance)
* :math:`\theta_0` (degrees)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of :math:`K_{BS0}` and :math:`K_{BS1}` are in energy/distance/radian.

Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
USER\_YAFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

**Default:** none
