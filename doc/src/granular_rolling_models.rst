Models for rolling friction in granular interactions
======================================================

.. _sds_rolling_model:

*sds* rolling model
-------------------

For *rolling sds*, rolling friction is computed via a
spring-dashpot-slider, using a 'pseudo-force' formulation, as described
by :ref:`Luding <Luding2008>`. Unlike the formulation in
:ref:`Marshall <Marshall2009>`, this allows for the required adjustment of
rolling displacement due to changes in the frame of reference of the
contacting pair.  The rolling pseudo-force is computed analogously to
the tangential force:

.. math::

   \mathbf{F}_{roll,0} =  k_{roll} \mathbf{\xi}_{roll}  - \gamma_{roll} \mathbf{v}_{roll}

Here, :math:`\mathbf{v}_{roll} = -R(\boldsymbol{\Omega}_i - \boldsymbol{\Omega}_j)
\times \mathbf{n}` is the relative rolling velocity, as given in
:ref:`Wang et al <Wang2015>` and :ref:`Luding <Luding2008>`. This differs from the
expressions given by :ref:`Kuhn and Bagi <Kuhn2004>` and used in :ref:`Marshall <Marshall2009>`;
see :ref:`Wang et al <Wang2015>` for details. The rolling displacement is given by:

.. math::

   \mathbf{\xi}_{roll} = \int_{t_0}^t \mathbf{v}_{roll} (\tau) \mathrm{d} \tau

A Coulomb friction criterion truncates the rolling pseudo-force if it
exceeds a critical value:

.. math::

   \mathbf{F}_{roll} =  \min(\mu_{roll} F_{n,0}, \|\mathbf{F}_{roll,0}\|)\mathbf{k}

Here, :math:`\mathbf{k} = \mathbf{v}_{roll}/\|\mathbf{v}_{roll}\|` is the direction of
the pseudo-force.  As with tangential displacement, the rolling
displacement is rescaled when the critical force is exceeded, so that
the spring length corresponds the critical force. Additionally, the
displacement is adjusted to account for rotations of the frame of
reference of the two contacting particles in a manner analogous to the
tangential displacement.

The rolling pseudo-force does not contribute to the total force on
either particle (hence 'pseudo'), but acts only to induce an equal and
opposite torque on each particle, according to:

.. math::

   \tau_{roll,i} =  R \mathbf{n} \times \mathbf{F}_{roll}

.. math::

   \tau_{roll,j} =  -\tau_{roll,i}

References
""""""""""

.. _Kuhn2004:

**(Kuhn and Bagi, 2005)** Kuhn, M. R., & Bagi, K. (2004). Contact
rolling and deformation in granular media.  International journal of
solids and structures, 41(21), 5793-5820.

.. _Wang2015:

**(Wang et al, 2015)** Wang, Y., Alonso-Marroquin, F., & Guo,
W. W. (2015).  Rolling and sliding in 3-D discrete element
models. Particuology, 23, 49-55.