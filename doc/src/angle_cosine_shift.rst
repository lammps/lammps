.. index:: angle_style cosine/shift
.. index:: angle_style cosine/shift/omp

angle_style cosine/shift command
================================

Accelerator Variants: *cosine/shift/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine/shift

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine/shift
   angle_coeff * 10.0 45.0

Description
"""""""""""

The *cosine/shift* angle style uses the potential

.. math::

   E = -\frac{U_{\text{min}}}{2} \left[ 1 + \cos(\theta-\theta_0) \right]

where :math:`\theta_0` is the equilibrium angle. The potential is bounded
between :math:`-U_{\text{min}}` and zero. In the neighborhood of the minimum
:math:`E = - U_{\text{min}} + U_{\text{min}}/4(\theta - \theta_0)^2` hence
the spring constant is :math:`\frac{U_{\text{min}}}{2}`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`U_{\text{min}}` (energy)
* :math:`\theta` (angle)

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

:doc:`angle_coeff <angle_coeff>`,
:doc:`angle_style cosine/shift/exp <angle_cosine_shift_exp>`

Default
"""""""

none
