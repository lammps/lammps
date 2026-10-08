Models for normal damping in granular interactions
==================================================

The normal force is augmented by a damping term of the following
general form:

.. math::

   \mathbf{F}_{n,damp} = -\eta_n \mathbf{v}_{n,rel}

Here, :math:`\mathbf{v}_{n,rel} = (\mathbf{v}_j - \mathbf{v}_i) \cdot \mathbf{n}\ \mathbf{n}`
is the component of relative velocity along :math:`\mathbf{n}`, where
:math:`\mathbf{n}` is the unit vector along the direction connecting
the two particle centers, discussed :ref:`here <normal_models_preamble>`.

Different damping models result in different expressions for :math:`\eta_n`:

.. _velocity_damping_model:

`velocity` damping model
------------------------

*Parameters:* none beyond the normal model damping coefficient
:math:`\eta_{n0}`.

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hooke 1000.0 50.0 tangential linear_nohistory 1.0 0.4 damping velocity

For *damping velocity*, the normal damping is simply equal to the
user-specified damping coefficient in the *normal* model:

.. math::

   \eta_n = \eta_{n0}

Here, :math:`\eta_{n0}` is the damping coefficient specified for the normal
contact model, in units of *mass*\ /\ *time*\ .

.. _mass_velocity_damping_model:

`mass_velocity` damping model
-----------------------------

*Parameters:* none beyond the normal model damping coefficient
:math:`\eta_{n0}`.

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hooke 1000.0 50.0 tangential linear_history 500.0 1.0 0.4 damping mass_velocity


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

*Parameters:* none beyond the normal-model damping coefficient
:math:`\eta_{n0}`.

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 damping viscoelastic



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

*Parameters:* none beyond the normal-model restitution coefficient
:math:`e`.

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hertz/material 1.0e8 0.3 0.3 tangential mindlin NULL 1.0 0.4 damping tsuji

The *tsuji* model is based on the work of :ref:`(Tsuji et al)
<Tsuji1992>`.  Here, the damping coefficient specified as part of the
normal model is interpreted as a restitution coefficient :math:`e`, assuming the
normal force is Hertzian.  The damping constant :math:`\eta_n` is given by:

.. math::

   \eta_n = \alpha (m_{eff}k_{nd})^{1/2}

where :math:`k_{nd}` is an effective harmonic stiffness equal to the
ratio of the normal force to the overlap.  For example, :math:`k_{nd} =
4/3Ea` for a Hertz contact model based on material parameters with
:math:`a` being the contact radius of :math:`\sqrt{\delta R}`.  For
Hooke, :math:`k_{nd}` is simply the spring constant or :math:`k_{n}`.
This damping model is not compatible with cohesive normal models such as
*JKR* or *DMT*.  The parameter :math:`\alpha` is related to the
restitution coefficient *e* according to:

.. math::

   \alpha / \sqrt{2} = 1.2728-4.2783e+11.087e^2-22.348e^3+27.467e^4-18.022e^5+4.8218e^6

The dimensionless coefficient of restitution :math:`e` specified as part
of the normal contact model parameters should be between 0 and 1, but no
error check is performed on this.
Using this damping model with normal contact models other than Hertz is
possible, but the resulting coefficient of restitution is not likely to
accurately match the specified value.

.. versionchanged:: 2Sep2026

This numerical solution is from :ref:`(Marshall, 2009) <Marshall2009_1>`
where the factor of :math:`\sqrt{2}` arises from a difference in convention
from Tsuji when defining :math:`\alpha` using either the mass vs. effective
mass. This factor was missing in earlier versions of LAMMPS.

.. _coeff_restitution_damping_model:

`coeff_restitution` damping model
-----------------------------------

*Parameters:* none beyond the normal-model restitution coefficient
:math:`e`.

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hertz/material 1.0e8 0.3 0.3 tangential mindlin NULL 1.0 0.4 damping coeff_restitution


The *coeff_restitution* model is useful when a specific normal coefficient of
restitution :math:`e` is required. It operates much like the *Tsuji* model,
but the normal coefficient of restitution :math:`e` is calculated using a different
expression (below) and is designed to work with Hooke and Hertz. Following
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

*Parameters:* :math:`d_{type}`

Example:

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * mdr 5.0e6 0.4 1.9e5 2.0 0.5 0.5 tangential linear_history 940.0 1.0 0.7 damping mdr 1


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

In this case, :math:`\eta_{n0}` is simply a dimensionless coefficient that scales
the overall damping coefficient.

The other supported option is :math:`d_{type} = 2`, which defines a simple damping model
similar to the *velocity* option:

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

