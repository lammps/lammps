.. index:: angle_style harmonic
.. index:: angle_style harmonic/intel
.. index:: angle_style harmonic/kk
.. index:: angle_style harmonic/omp

angle_style harmonic command
============================

Accelerator Variants: *harmonic/intel*, *harmonic/kk*, *harmonic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style harmonic

Examples
""""""""

.. code-block:: LAMMPS

   angle_style harmonic
   angle_coeff 1 300.0 107.0

Description
"""""""""""

The *harmonic* angle style uses the potential

.. math::

   E = K (\theta - \theta_0)^2

where :math:`\theta_0` is the equilibrium value of the angle, and :math:`K` is a
prefactor.  Note that the usual 1/2 factor is included in :math:`K`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`\theta_0` (degrees)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence :math:`K` is effectively energy per
radian\^2.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none
