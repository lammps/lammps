.. index:: angle_style cosine/delta
.. index:: angle_style cosine/delta/omp

angle_style cosine/delta command
================================

Accelerator Variants: *cosine/delta/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine/delta

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine/delta
   angle_coeff 2*4 75.0 100.0

Description
"""""""""""

The *cosine/delta* angle style uses the potential

.. math::

   E = K [1 - \cos(\theta - \theta_0)]

where :math:`\theta_0` is the equilibrium value of the angle, and :math:`K` is a
prefactor.  Note that the usual 1/2 factor is included in :math:`K`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`\theta_0` (degrees)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to radians
internally.

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

:doc:`angle_coeff <angle_coeff>`, :doc:`angle_style cosine/squared <angle_cosine_squared>`

Default
"""""""

none
