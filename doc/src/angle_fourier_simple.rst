.. index:: angle_style fourier/simple
.. index:: angle_style fourier/simple/omp

angle_style fourier/simple command
==================================

Accelerator Variants: *fourier/simple/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style fourier/simple

Examples
""""""""

.. code-block:: LAMMPS

   angle_style fourier/simple
   angle_coeff 100.0 -1.0 1.0

Description
"""""""""""

The *fourier/simple* angle style uses the potential

.. math::

   E = K [ 1.0 + c \cos ( n \theta) ]

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`c` (real)
* :math:`n` (real)

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
