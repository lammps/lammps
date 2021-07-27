.. index:: dihedral_style cosine/shift/exp
.. index:: dihedral_style cosine/shift/exp/omp

dihedral_style cosine/shift/exp command
=======================================

Accelerator Variants: *cosine/shift/exp/omp*

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style cosine/shift/exp

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style cosine/shift/exp
   dihedral_coeff 1 10.0 45.0 2.0

Description
"""""""""""

The *cosine/shift/exp* dihedral style uses the potential

.. math::

   E = -U_{min}\frac{e^{-a U(\theta,\theta_0)}-1}{e^a-1} \quad\mbox{with}\quad U(\theta,\theta_0)=-0.5 \left(1+\cos(\theta-\theta_0) \right)

where :math:`U_{min}`, :math:`\theta`, and :math:`a` are defined for
each dihedral type.

The potential is bounded between :math:`\left[-U_{min}:0\right]` and the minimum is located
at the angle :math:`\theta_0`. The a parameter can be both positive or negative
and is used to control the spring constant at the equilibrium.

The spring constant is given by :math:`k=a e^a \frac{U_{min}}{2 \left(e^a-1\right)}`.
For :math:`a>3` and  :math:`\frac{k}{U_{min}} = \frac{a}{2}` to better than 5% relative error. For negative
values of the a parameter, the spring constant is essentially zero,
and anharmonic terms takes over. The potential is furthermore well
behaved in the limit :math:`a \rightarrow 0`, where it has been implemented to linear
order in :math:`a` for :math:`a < 0.001`.

The following coefficients must be defined for each dihedral type via
the :doc:`dihedral_coeff <dihedral_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`U_{min}` (energy)
* :math:`\theta` (angle)
* :math:`a` (real number)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This dihedral style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`,
:doc:`angle_style cosine/shift/exp <angle_cosine_shift_exp>`

Default
"""""""

none
