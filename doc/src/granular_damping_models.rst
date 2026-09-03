Models for normal damping in granular interactions
==================================================

The normal force is augmented by a damping term of the following 
general form:

.. math::

   \mathbf{F}_{n,damp} = -\eta_n \mathbf{v}_{n,rel}

Here, :math:`\mathbf{v}_{n,rel} = (\mathbf{v}_j - \mathbf{v}_i) \cdot \mathbf{n}\ \mathbf{n}` 
is the component of relative velocity along :math:`\mathbf{n}`, where
:math:`\mathbf{n}` is the unit vector along the direction connecting 
the two particle centers, discussed :ref:`here <normal_models_preamble>`. .

Different damping models result in different expressions for :math:`\eta_n`:

.. _velocity_damping_model:

`velocity` damping model
------------------------

For *damping velocity*, the normal damping is simply equal to the
user-specified damping coefficient in the *normal* model:

.. math::

   \eta_n = \eta_{n0}

Here, :math:`\eta_{n0}` is the damping coefficient specified for the normal
contact model, in units of *mass*\ /\ *time*\ .

.. _mass_velocity_damping_model:

`mass_velocity` damping model
-----------------------------

For *damping mass_velocity*, the normal damping is given by:

.. math::

   \eta_n = \eta_{n0} m_{eff}

Here, :math:`\eta_{n0}` is the damping coefficient specified for the normal
contact model, in units of 1/\ *time* and
:math:`m_{eff} = m_i m_j/(m_i + m_j)` is the effective mass.
Use *damping mass_velocity* to reproduce the damping behavior of
*pair gran/hooke/\**.

.. _viscoelastic_damping_model:

`viscoelastic` damping model
----------------------------

The *damping viscoelastic* model is based on the viscoelastic
treatment of :ref:`(Brilliantov et al) <Brill1996>`, where the normal
damping is given by:

.. math::

   \eta_n = \eta_{n0}\ a m_{eff}

Here, *a* is the contact radius, given by :math:`a =\sqrt{R\delta}`
for all models except *jkr*, for which it is given implicitly according
to :math:`\delta = a^2/R - 2\sqrt{\pi \gamma a/E}`.  For *damping viscoelastic*,
:math:`\eta_{n0}` is in units of 1/(\ *time*\ \*\ *distance*\ ).


.. _tsuji_damping_model:

`tsuji` damping model
---------------------


The *tsuji* model is based on the work of :ref:`(Tsuji et al) <Tsuji1992>`.
Here, the damping coefficient specified as part of the normal model is interpreted
as a restitution coefficient :math:`e`. The damping constant :math:`\eta_n` is
given by:

.. math::

   \eta_n = \alpha (m_{eff}k_n)^{1/2}

For normal contact models based on material parameters, :math:`k_n = 4/3Ea`. This
damping model is not compatible with cohesive normal models such as *JKR* or *DMT*.
The parameter :math:`\alpha` is related to the restitution coefficient *e*
according to:

.. math::

   \alpha = 1.2728-4.2783e+11.087e^2-22.348e^3+27.467e^4-18.022e^5+4.8218e^6

The dimensionless coefficient of restitution :math:`e` specified as part
of the normal contact model parameters should be between 0 and 1, but
no error check is performed on this.

.. _coeff_restitution_damping_model:

`coeff_restitution`` damping model
-----------------------------------

The *coeff_restitution* model is useful when a specific normal coefficient of
restitution :math:`e` is required. It operates much like the *Tsuji* model
but, the normal coefficient of restitution :math:`e` is specified as an input
in place of the usual :math:`\eta_{n0}` value in the normal model. Following
the approach of :ref:`(Brilliantov et al) <Brill1996>`, when using the *hooke*
normal model, *coeff_restitution* then calculates the damping coefficient as:

.. math::

   \eta_n = \sqrt{\frac{4m_{eff}k_{nd}}{1+\left( \frac{\pi}{\log(e)}\right)^2}} ,

where :math:`k_{nd}` is the same stiffness defined in the above *Tsuji* model.
For any other normal model, e.g. the *hertz* and *hertz/material* models, the damping
coefficient is:

.. math::

   \eta_n = -2\sqrt{\frac{5}{6}}\frac{\log(e)}{\sqrt{\pi^2+(\log(e))^2}}\sqrt{\frac{3}{2}k_{nd} m_{eff}} ,

Since *coeff_restitution* accounts for the effective mass, effective radius,
and pairwise overlaps (except when used with the *hooke* normal model) when calculating
the damping coefficient, it accurately reproduces the specified coefficient of
restitution for both monodisperse and polydisperse particle pairs.  This damping
model is not compatible with cohesive normal models such as *JKR* or *DMT*.

.. _mdr_damping_model:

mdr class of damping models
---------------------------

The *mdr* damping class contains multiple damping models that can be toggled between
by specifying different integer values for the :math:`d_{type}` input parameter. This
damping option is only compatible with the normal *mdr* contact model.

Setting :math:`d_{type} = 1` is the suggested damping option. This specifies a damping
model that takes into account the contact stiffness :math:`k_{mdr}` calculated
by the normal *mdr* contact model to determine the damping coefficient:

.. math::

   \eta_n = \eta_{n0} (m_{eff}k_{mdr})^{1/2},

where :math:`k_{mdr}` is proportional to contact radius :math:`a_{mdr}` tracked by the
normal *mdr* contact model:

.. math::

   k_{mdr} = 2 E_{eff} a_{mdr}.

In this case, :math:`\eta_{n0}` is simply a dimensionless coefficient that scales the
the overall damping coefficient.

The other supported option is :math:`d_{type} = 2`, which defines a simple damping model
similar to the *velocity* option

.. math::

   \eta_n = \eta_{n0},

but has additional checks to avoid non-physical damping after plastic deformation.

The total normal force is computed as the sum of the elastic and
damping components:

.. math::

   \mathbf{F}_n = \mathbf{F}_{n,e} + \mathbf{F}_{n,damp}

--------------

References
""""""""""

.. _Brill1996:

**(Brilliantov et al, 1996)** Brilliantov, N. V., Spahn, F., Hertzsch,
J. M., & Poschel, T. (1996).  Model for collisions in granular
gases. Physical review E, 53(5), 5382.

.. _Tsuji1992:

**(Tsuji et al, 1992)** Tsuji, Y., Tanaka, T., & Ishida,
T. (1992). Lagrangian numerical simulation of plug flow of
cohesionless particles in a horizontal pipe. Powder technology, 71(3),
239-250.

