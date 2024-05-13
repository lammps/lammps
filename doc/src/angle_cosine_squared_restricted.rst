.. index:: angle_style cosine/squared/restricted
.. index:: angle_style cosine/squared/restricted/omp

angle_style cosine/squared/restricted command
=============================================

Accelerator Variants: *cosine/squared/restricted/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine/squared/restricted

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine/squared/restricted
   angle_coeff 2*4 75.0 100.0

Description
"""""""""""

.. versionadded:: 17Apr2024

The *cosine/squared/restricted* angle style uses the potential

.. math::

   E = K [\cos(\theta) - \cos(\theta_0)]^2 / \sin^2(\theta)

, which is commonly used in the MARTINI force field,
where :math:`\theta_0` is the equilibrium value of the angle, and :math:`K`
is a prefactor.  Note that the usual 1/2 factor is included in :math:`K`.

See :ref:`(Bulacu) <restricted-Bulacu>` for a description of the restricted angle for the MARTINI force field.

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

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none

----------

.. _restricted-Bulacu:

**(Bulacu)** Bulacu, Goga, Zhao, Rossi, Monticelli, Periole, Tieleman, Marrink, J Chem Theory Comput, 9, 3282-3292
(2013).
