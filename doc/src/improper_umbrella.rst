.. index:: improper_style umbrella
.. index:: improper_style umbrella/omp

improper_style umbrella command
===============================

Accelerator Variants: *umbrella/omp*

Syntax
""""""

.. code-block:: LAMMPS

   improper_style umbrella

Examples
""""""""

.. code-block:: LAMMPS

   improper_style umbrella
   improper_coeff 1 100.0 180.0

Description
"""""""""""

The *umbrella* improper style uses the following potential, which is
commonly referred to as a classic inversion and used in the
:doc:`DREIDING <Howto_bioFF>` force field:

.. math::

   E = & \frac{1}{2}K\left( \frac{1}{\sin\omega_0}\right) ^2 \left( \cos\omega - \cos\omega_0\right) ^2 \qquad \omega_0 \neq 0^o \\
   E = & K\left( 1-cos\omega\right)  \qquad \omega_0 = 0^o

where :math:`K` is the force constant and :math:`\omega` is the angle between the IL
axis and the IJK plane:

.. image:: JPG/umbrella.jpg
   :align: center

If :math:`\omega_0 = 0` the potential term has a minimum for the planar
structure.  Otherwise it has two minima at :math:`\omega +/- \omega_0`,
with a barrier in between.

See :ref:`(Mayo) <umbrella-Mayo>` for a description of the DREIDING force field.

The following coefficients must be defined for each improper type via
the :doc:`improper_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)
* :math:`\omega_0` (degrees)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none

----------

.. _umbrella-Mayo:

**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990),
