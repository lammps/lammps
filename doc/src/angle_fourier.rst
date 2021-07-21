.. index:: angle_style fourier
.. index:: angle_style fourier/omp

angle_style fourier command
===========================

Accelerator Variants: *fourier/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style fourier

Examples
""""""""

.. code-block:: LAMMPS

   angle_style fourier
   angle_coeff 75.0 1.0 1.0 1.0

Description
"""""""""""

The *fourier* angle style uses the potential

.. math::

   E = K [C_0 + C_1 \cos ( \theta) + C_2 \cos( 2 \theta) ]

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`C_0` (real)
* :math:`C_1` (real)
* :math:`C_2` (real)

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
