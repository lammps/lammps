.. index:: angle_style quartic
.. index:: angle_style quartic/omp

angle_style quartic command
===========================

Accelerator Variants: *quartic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style quartic

Examples
""""""""

.. code-block:: LAMMPS

   angle_style quartic
   angle_coeff 1 129.1948 56.8726 -25.9442 -14.2221

Description
"""""""""""

The *quartic* angle style uses the potential

.. math::

   E = K_2 (\theta - \theta_0)^2 + K_3 (\theta - \theta_0)^3 + K_4 (\theta - \theta_0)^4

where :math:`\theta_0` is the equilibrium value of the angle, and :math:`K` is a
prefactor.  Note that the usual 1/2 factor is included in :math:`K`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`\theta_0` (degrees)
* :math:`K_2` (energy)
* :math:`K_3` (energy)
* :math:`K_4` (energy)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence the various :math:`K` are effectively energy
per radian\^2 or radian\^3 or radian\^4.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none
