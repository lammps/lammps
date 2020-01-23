.. index:: angle_style mm3

angle_style mm3 command
=======================

Syntax
""""""


.. code-block:: LAMMPS

   angle_style mm3

Examples
""""""""


.. code-block:: LAMMPS

   angle_style mm3
   angle_coeff 1 100.0 107.0

Description
"""""""""""

The *mm3* angle style uses the potential that is anharmonic in the angle
as defined in :ref:`(Allinger) <mm3-allinger1989>`

.. math::

   E = K (\theta - \theta_0)^2 \left[ 1 - 0.014(\theta - \theta_0) + 5.6(10)^{-5} (\theta - \theta_0)^2 - 7.0(10)^{-7} (\theta - \theta_0)^3 + 9(10)^{-10} (\theta - \theta_0)^4 \right]


where :math:`\theta_0` is the equilibrium value of the angle, and :math:`K` is a
prefactor. The anharmonic prefactors have units :math:`\deg^{-n}`, for example
:math:`-0.014 \deg^{-1}`, :math:`5.6 \cdot 10^{-5} \deg^{-2}`, ...

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/radian\^2)
* :math:`\theta_0` (degrees)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of :math:`K` are in energy/radian\^2.

Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
USER\_YAFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

**Default:** none
