Models for tangential contact in granular interactions
======================================================

.. _linear_nohistory:

linear_nohistory tangential model
---------------------------------

For *tangential linear_nohistory*, a simple velocity-dependent Coulomb
friction criterion is used, which mimics the behavior of the *pair
gran/hooke* style. The tangential force :math:`\mathbf{F}_t` is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|\mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}

The tangential damping force :math:`\mathbf{F}_\mathrm{t,damp}` is given by:

.. math::

   \mathbf{F}_\mathrm{t,damp} = -\eta_t \mathbf{v}_{t,rel}

The tangential damping prefactor :math:`\eta_t` is calculated by scaling
the normal damping :math:`\eta_n` (see above):

.. math::

   \eta_t = -x_{\gamma,t} \eta_n

The normal damping prefactor :math:`\eta_n` is determined by the choice
of the *damping* keyword, as discussed above.  Thus, the *damping*
keyword also affects the tangential damping.  The parameter
:math:`x_{\gamma,t}` is a scaling coefficient. Several works in the
literature use :math:`x_{\gamma,t} = 1` (:ref:`Marshall <Marshall2009>`,
:ref:`Tsuji et al <Tsuji1992>`, :ref:`Silbert et al <Silbert2001>`).  The relative
tangential velocity at the point of contact is given by
:math:`\mathbf{v}_{t, rel} = \mathbf{v}_{t} - (R_i\mathbf{\Omega}_i + R_j\mathbf{\Omega}_j) \times \mathbf{n}`, where :math:`\mathbf{v}_{t} = \mathbf{v}_r - \mathbf{v}_r\cdot\mathbf{n}\ \mathbf{n}`,
:math:`\mathbf{v}_r = \mathbf{v}_j - \mathbf{v}_i` .
The direction of the applied force is :math:`\mathbf{t} = \mathbf{v_{t,rel}}/\|\mathbf{v_{t,rel}}\|` .

The normal force value :math:`F_{n0}` used to compute the critical force
depends on the form of the contact model. For non-cohesive models
(\ *hertz*, *hertz/material*, *hooke*\ ), it is given by the magnitude of
the normal force:

.. math::

   F_{n0} = \|\mathbf{F}_n\|

For cohesive models such as *jkr* and *dmt*, the critical force is
adjusted so that the critical tangential force approaches :math:`\mu_t F_{pulloff}`, see :ref:`Marshall <Marshall2009>`, equation 43, and
:ref:`Thornton <Thornton1991>`.  For both models, :math:`F_{n0}` takes the
form:

.. math::

   F_{n0} = \|\mathbf{F}_{ne} + 2 F_{pulloff}\|

Where :math:`F_{pulloff} = 3\pi \gamma R` for *jkr*, and
:math:`F_{pulloff} = 4\pi \gamma R` for *dmt*\ .

The remaining tangential options all use accumulated tangential
displacement (i.e. contact history), except for the options
*mindlin/force* and *mindlin_rescale/force*, that use accumulated
tangential force instead, and are discussed further below.
The accumulated tangential displacement is discussed in details below
in the context of the *linear_history* option. The same treatment of
the accumulated displacement applies to the other options as well.

For *tangential linear_history*, the tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t\mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}

Here, :math:`\mathbf{\xi}` is the tangential displacement accumulated
during the entire duration of the contact:

.. math::

   \mathbf{\xi} = \int_{t0}^t \mathbf{v}_{t,rel}(\tau) \mathrm{d}\tau

This accumulated tangential displacement must be adjusted to account
for changes in the frame of reference of the contacting pair of
particles during contact. This occurs due to the overall motion of the
contacting particles in a rigid-body-like fashion during the duration
of the contact. There are two modes of motion that are relevant: the
'tumbling' rotation of the contacting pair, which changes the
orientation of the plane in which tangential displacement occurs; and
'spinning' rotation of the contacting pair about the vector connecting
their centers of mass (:math:`\mathbf{n}`).  Corrections due to the
former mode of motion are made by rotating the accumulated
displacement into the plane that is tangential to the contact vector
at each step, or equivalently removing any component of the tangential
displacement that lies along :math:`\mathbf{n}`, and rescaling to
preserve the magnitude.  This follows the discussion in
:ref:`Luding <Luding2008>`, see equation 17 and relevant discussion in that
work:

.. math::

   \mathbf{\xi} = \left(\mathbf{\xi'} - (\mathbf{n} \cdot \mathbf{\xi'})\mathbf{n}\right) \frac{\|\mathbf{\xi'}\|}{\|\mathbf{\xi'} - (\mathbf{n}\cdot\mathbf{\xi'})\mathbf{n}\|}

Here, :math:`\mathbf{\xi'}` is the accumulated displacement prior to the
current time step and :math:`\mathbf{\xi}` is the corrected
displacement. Corrections to the displacement due to the second mode
of motion described above (rotations about :math:`\mathbf{n}`) are not
currently implemented, but are expected to be minor for most
simulations.

Furthermore, when the tangential force exceeds the critical force, the
tangential displacement is re-scaled to match the value for the
critical force (see :ref:`Luding <Luding2008>`, equation 20 and related
discussion):

.. math::

   \mathbf{\xi} = -\frac{1}{k_t}\left(\mu_t F_{n0}\mathbf{t} - \mathbf{F}_{t,damp}\right)

The tangential force is added to the total normal force (elastic plus
damping) to produce the total force on the particle. The tangential
force also acts at the contact point (defined as the center of the
overlap region) to induce a torque on each particle according to:

.. math::

   \mathbf{\tau}_i = -(R_i - 0.5 \delta) \mathbf{n} \times \mathbf{F}_t

.. math::

   \mathbf{\tau}_j = -(R_j - 0.5 \delta) \mathbf{n} \times \mathbf{F}_t

For *tangential mindlin*, the :ref:`Mindlin <Mindlin1949>` no-slip solution
is used which differs from the *linear_history* option by an additional factor
of :math:`a`, the radius of the contact region. The tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t a \mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}


Here, :math:`a` is the radius of the contact region, given by :math:`a =\sqrt{R\delta}`
for all normal contact models, except for *jkr*, where it is given
implicitly by :math:`\delta = a^2/R - 2\sqrt{\pi \gamma a/E}`, see
discussion above. To match the Mindlin solution, one should set
:math:`k_t = 8G_{eff}`, where :math:`G_{eff}` is the effective shear modulus given by:

.. math::

   G_{eff} = \left(\frac{2-\nu_i}{G_i} + \frac{2-\nu_j}{G_j}\right)^{-1}

where :math:`G` is the shear modulus, related to Young's modulus :math:`E`
and Poisson's ratio :math:`\nu` by :math:`G = E/(2(1+\nu))`. This can also be
achieved by specifying *NULL* for :math:`k_t`, in which case a
normal contact model that specifies material parameters :math:`E` and
:math:`\nu` is required (e.g. *hertz/material*, *dmt* or *jkr*\ ). In this
case, mixing of the shear modulus for different particle types *i* and
*j* is done according to the formula above.

.. note::

   The radius of the contact region :math:`a` depends on the normal overlap.
   As a result, the tangential force for *mindlin* can change due to
   a variation in normal overlap, even with no change in tangential displacement.

For *tangential mindlin/force*, the accumulated elastic tangential force
characterizes the contact history, instead of the accumulated tangential
displacement. This prevents the dependence of the tangential force on the
normal overlap as noted above. The tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|\mathbf{F}_{te} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}

The increment of the elastic component of the tangential force
:math:`\mathbf{F}_{te}` is given by:

.. math::

   \mathrm{d}\mathbf{F}_{te} = -k_t a \mathbf{v}_{t,rel} \mathrm{d}\tau

The changes in frame of reference of the contacting pair of particles during
contact are accounted for by the same formula as above, replacing the
accumulated tangential displacement :math:`\xi`, by the accumulated tangential
elastic force :math:`F_{te}`. When the tangential force exceeds the critical
force, the tangential force is directly re-scaled to match the value for
the critical force:

.. math::

   \mathbf{F}_{te} = - \mu_t F_{n0}\mathbf{t} + \mathbf{F}_{t,damp}

The same rules as those described for *mindlin* apply regarding the tangential
stiffness and mixing of the shear modulus for different particle types.

The *mindlin_rescale* option uses the same form as *mindlin*, but the
magnitude of the tangential displacement is re-scaled as the contact
unloads, i.e. if :math:`a < a_{t_{n-1}}`:

.. math::

   \mathbf{\xi} = \mathbf{\xi_{t_{n-1}}} \frac{a}{a_{t_{n-1}}}

Here, :math:`t_{n-1}` indicates the value at the previous time
step. This rescaling accounts for the fact that a decrease in the
contact area upon unloading leads to the contact being unable to
support the previous tangential loading, and spurious energy is
created without the rescaling above (:ref:`Walton <WaltonPC>` ).

.. note::

   For *mindlin*, a decrease in the tangential force already occurs as the
   contact unloads, due to the dependence of the tangential force on the normal
   force described above. By re-scaling :math:`\xi`, *mindlin_rescale*
   effectively re-scales the tangential force twice, i.e., proportionally to
   :math:`a^2`. This peculiar behavior results from use of the accumulated
   tangential displacement to characterize the contact history. Although
   *mindlin_rescale* remains available for historic reasons and backward
   compatibility purposes, it should be avoided in favor of *mindlin_rescale/force*.

The *mindlin_rescale/force* option uses the same form as *mindlin/force*,
but the magnitude of the tangential elastic force is re-scaled as the contact
unloads, i.e. if :math:`a < a_{t_{n-1}}`:

.. math::

   \mathbf{F}_{te} = \mathbf{F}_{te, t_{n-1}} \frac{a}{a_{t_{n-1}}}

This approach provides a better approximation of the :ref:`Mindlin-Deresiewicz <Mindlin1953>`
laws and is more consistent than *mindlin_rescale*. See discussions in
:ref:`Thornton et al, 2013 <Thornton2013>`, particularly equation 18(b) of that
work and associated discussion, and :ref:`Agnolin and Roux, 2007 <AgnolinRoux2007>`,
particularly Appendix A.
