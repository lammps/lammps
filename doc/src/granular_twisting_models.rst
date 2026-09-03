Models for twisting friction in granular interactions
======================================================

A history-dependent spring-dashpot-slider is used to compute the 
twisting torque. Because twisting displacement is a scalar, there 
is no need to adjust for changes in the frame of reference due to 
rotations of the particle pair. The formulation in :ref:`Marshall <Marshall2009>` 
therefore provides the most straightforward treatment:

.. math::

   \tau_{twist,0} = -k_{twist}\xi_{twist} - \gamma_{twist}\Omega_{twist}

Here :math:`\xi_{twist} = \int_{t_0}^t \Omega_{twist} (\tau) \mathrm{d}\tau` is the twisting angular displacement, and
:math:`\Omega_{twist} = (\mathbf{\Omega}_i - \mathbf{\Omega}_j) \cdot \mathbf{n}` is the relative twisting angular velocity. The torque
is then truncated according to:

.. math::

   \tau_{twist} = \min(\mu_{twist} F_{n,0}, \tau_{twist,0})

Similar to the sliding and rolling displacement, the angular
displacement is rescaled so that it corresponds to the critical value
if the twisting torque exceeds this critical value:

.. math::

   \xi_{twist} = \frac{1}{k_{twist}} (\mu_{twist} F_{n,0}sgn(\Omega_{twist}) - \gamma_{twist}\Omega_{twist})

The twisting torque on each particle is then given by:

.. math::

   \mathbf{\tau}_{twist,i} = \tau_{twist}\mathbf{n}

.. math::

   \mathbf{\tau}_{twist,j} = -\mathbf{\tau}_{twist,i}

.. _sds_twisting_model:

*sds* twisting friction model
-----------------------------

For *twisting sds*, the coefficients :math:`k_{twist}, \gamma_{twist}`
and :math:`\mu_{twist}` are simply the user input parameters that follow
the *twisting sds* keywords in the *pair_coeff* command.

.. _marshall_twisting_model:

*marshall* twisting history model
---------------------------------

For *twisting marshall*, the coefficients :math:`k_{twist}, \gamma_{twist}`
and :math:`\mu_{twist}` are expressed in terms of the
tangential friction coefficients :math:`k_t, \eta_t, \mu_t`, as discussed in
:ref:`Marshall <Marshall2009>` (see equations 32 and 33 of that work, and
the doc page for :doc:`granular_tangential_models`):

.. math::

   k_{twist} = 0.5k_ta^2

.. math::

   \eta_{twist} = 0.5\eta_ta^2

.. math::

   \mu_{twist} = \frac{2}{3}a\mu_t


**(Marshall, 2009)** Marshall, J. S. (2009). Discrete-element modeling
of particulate aerosol flows.  Journal of Computational Physics,
228(5), 1541-1561.

