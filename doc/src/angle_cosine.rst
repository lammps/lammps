.. index:: angle_style cosine
.. index:: angle_style cosine/omp
.. index:: angle_style cosine/kk

angle_style cosine command
==========================

Accelerator Variants: *cosine/omp*, *cosine/kk*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine
   angle_coeff * 75.0

Description
"""""""""""

The *cosine* angle style uses the potential

.. math::

   E = K [1 + \cos(\theta)]

where :math:`K` is defined for each angle type.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none
