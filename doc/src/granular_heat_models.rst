Models for heat transport in granular interactions
======================================================

.. _radius_heat_model:

*radius* heat model
--------------------
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

For *heat* *area*, the heat
:math:`Q` conducted between two particles is given by

.. math::

   Q = h_{s} A \Delta T


where :math:`\Delta T` is the difference in the two particles' temperature,
:math:`h_{s}` is a non-negative numeric value for the heat transfer
coefficient (in units of power/(area*temperature)), and :math:`A=\pi a^2` is
the area of the contact and depends on the normal force model.