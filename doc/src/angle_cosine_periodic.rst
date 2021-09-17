.. index:: angle_style cosine/periodic
.. index:: angle_style cosine/periodic/omp

angle_style cosine/periodic command
===================================

Accelerator Variants: *cosine/periodic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine/periodic

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine/periodic
   angle_coeff * 75.0 1 6

Description
"""""""""""

The *cosine/periodic* angle style uses the following potential, which may be
particularly used for organometallic systems where :math:`n` = 4 might be used
for an octahedral complex and :math:`n` = 3 might be used for a trigonal
center:

.. math::

   E = \frac{2.0}{n^2} * C \left[ 1 - B(-1)^n\cos\left( n\theta\right) \right]

where :math:`C`, :math:`B` and :math:`n` are coefficients defined for each angle type.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`C` (energy)
* :math:`B` = 1 or -1
* :math:`n` = 1, 2, 3, 4, 5 or 6 for periodicity

Note that the prefactor :math:`C` is specified as coefficient and not the overall force
constant :math:`K = \frac{2 C}{n^2}`.  When :math:`B = 1`, it leads to a minimum for the
linear geometry.  When :math:`B = -1`, it leads to a maximum for the linear geometry.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none

