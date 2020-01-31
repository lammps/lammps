.. index:: angle_style cosine/shift

angle_style cosine/shift command
=================================

angle_style cosine/shift/omp command
====================================

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


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
USER-MISC package.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`,
:doc:`angle_style cosine/shift/exp <angle_cosine_shift_exp>`

**Default:** none
