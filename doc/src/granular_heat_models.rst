Models for heat transport in granular interactions
======================================================

.. _radius_heat_model:

*radius* heat model
--------------------

*Parameters:* :math:`k_s`

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 heat radius 0.1

For *heat* *radius*, the heat
:math:`Q` conducted between two particles is given by

.. math::

   Q = 2 k_{s} a \Delta T

where :math:`\Delta T` is the difference in the two particles' temperature,
:math:`k_{s}` is a non-negative numeric value for the conductivity (in units
of power/(length*temperature)), and :math:`a` is the radius of the contact and
depends on the normal force model. This is the model proposed by
:ref:`Vargas and McCarthy <VargasMcCarthy2001>`.

.. _area_heat_model:

*area* heat model
------------------

*Parameters:* :math:`h_s`

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 heat area 0.1

For *heat* *area*, the heat
:math:`Q` conducted between two particles is given by

.. math::

   Q = h_{s} A \Delta T


where :math:`\Delta T` is the difference in the two particles' temperature,
:math:`h_{s}` is a non-negative numeric value for the heat transfer
coefficient (in units of power/(area*temperature)), and :math:`A=\pi a^2` is
the area of the contact and depends on the normal force model.
