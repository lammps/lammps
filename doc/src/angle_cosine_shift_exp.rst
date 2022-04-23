.. index:: angle_style cosine/shift/exp
.. index:: angle_style cosine/shift/exp/omp

angle_style cosine/shift/exp command
====================================

Accelerator Variants: *cosine/shift/exp/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine/shift/exp

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine/shift/exp
   angle_coeff * 10.0 45.0 2.0

Description
"""""""""""

The *cosine/shift/exp* angle style uses the potential

.. math::

   E = -U_{\text{min}} \frac{e^{-a U(\theta,\theta_0)}-1}{e^a-1} \quad \text{with} \quad U(\theta,\theta_0) = -0.5 \left(1+\cos(\theta-\theta_0) \right)

where :math:`U_{\text{min}}`, :math:`\theta`, and :math:`a` are defined for each angle type.

The potential is bounded between :math:`[-U_{\text{min}}, 0]` and the minimum is
located at the angle :math:`\theta_0`. The a parameter can be both positive or
negative and is used to control the spring constant at the
equilibrium.

The spring constant is given by :math:`k = A \exp(A) U_{\text{min}} / [2 (\exp(a)-1)]`.
For :math:`a > 3`, :math:`\frac{k}{U_{\text{min}}} = \frac{a}{2}` to better than 5% relative error. For negative
values of the :math:`a` parameter, the spring constant is essentially zero,
and anharmonic terms takes over. The potential is furthermore well
behaved in the limit :math:`a \rightarrow 0`, where it has been implemented to linear
order in :math:`a` for :math:`a < 0.001`. In this limit the potential reduces to the
cosineshifted potential.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`U_min` (energy)
* :math:`\theta` (angle)
* :math:`A` (real number)

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
:doc:`angle_style cosine/shift <angle_cosine_shift>`,
:doc:`dihedral_style cosine/shift/exp <dihedral_cosine_shift_exp>`

Default
"""""""

none
