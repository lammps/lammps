.. index:: improper_style fourier
.. index:: improper_style fourier/omp

improper_style fourier command
==============================

Accelerator Variants: *fourier/omp*

Syntax
""""""

.. code-block:: LAMMPS

   improper_style fourier

Examples
""""""""

.. code-block:: LAMMPS

   improper_style fourier
   improper_coeff 1 100.0 0.0 1.0 0.5 1

Description
"""""""""""

The *fourier* improper style uses the following potential:

.. math::

   E = K [C_0 + C_1 \cos ( \omega) + C_2 \cos( 2 \omega) ]

where K is the force constant, C0, C1, C2 are dimensionless coefficients,
and omega is the angle between the IL axis and the IJK plane:

.. image:: JPG/umbrella.jpg
   :align: center

If all parameter (see below) is not zero, the all the three possible angles will taken in account.

The following coefficients must be defined for each improper type via
the :doc:`improper_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)
* :math:`C_0` (unitless)
* :math:`C_1` (unitless)
* :math:`C_2` (unitless)
* all  (0 or 1, optional)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none
