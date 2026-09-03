Models for tangential friction in granular interactions
========================================================

.. _tangential_models_preamble:

Tangential models apply a force :math:`\mathbf{F}_t = \mathbf{F}_{i,t} = -\mathbf{F}_{j,t}` 
on each particle :math:`i,j` in a direction tangential
to the surfaces of the particles, and induce torques 
:math:`\tau_i,\tau_j`. The tangential force is added 
to the total normal force (elastic plus
damping) to produce the total force on each particle. The torques 
act at the contact point, defined as the center of the
overlap region:

.. math::
   
   \mathbf{\tau}_i = -(R_i - 0.5 \delta) \mathbf{n} \times \mathbf{F}_t

.. math::

   \mathbf{\tau}_j = -(R_j - 0.5 \delta) \mathbf{n} \times \mathbf{F}_t

The direction of the tangential force :math:`\mathbf{t}` is based on 
the relative tangential velocity :math:`\mathbf{v}_{t,rel}` at the 
point of contact:

.. math::

   \mathbf{v}_{t, rel} = \mathbf{v}_{t} - (R_i\mathbf{\Omega}_i + R_j\mathbf{\Omega}_j) \times \mathbf{n}

Here, :math:`\mathbf{v}_{t} = \mathbf{v}_r - \mathbf{v}_r\cdot\mathbf{n}\ \mathbf{n}`,
:math:`\mathbf{v}_r = \mathbf{v}_j - \mathbf{v}_i`, and :math:`\mathbf{n}` is 
the unit vector along the direction connecting the two 
particle centers, discussed :ref:`here <normal_models_preamble>`. 

The direction of the tangential force is then:

.. math::

   \mathbf{t} = \mathbf{v_{t,rel}}/\|\mathbf{v_{t,rel}}\|


.. _tangential_damping_discussion:

Tangential damping
--------------------------------

The tangential damping force :math:`\mathbf{F}_\mathrm{t,damp}` used by
all tangential models has the generic form:

.. math::

   \mathbf{F}_\mathrm{t,damp} = -\eta_t \mathbf{v}_{t,rel}

The tangential damping prefactor :math:`\eta_t` is calculated by scaling
the normal damping :math:`\eta_n` (see :doc:`granular_damping_models`):

.. math::

   \eta_t = -x_{\gamma,t} \eta_n

The normal damping prefactor :math:`\eta_n` is determined by the choice
of the *damping* keyword, as discussed :doc:`here <granular_damping_models>`.  
Thus, the *damping* keyword to the pair coeff command also 
affects the tangential damping.  The parameter :math:`x_{\gamma,t}` is a 
scaling coefficient. Several works in the literature use 
:math:`x_{\gamma,t} = 1` (:ref:`Marshall <Marshall2009>`,
:ref:`Tsuji et al <Tsuji1992>`, :ref:`Silbert et al <Silbert2001>`).  

.. _tangential_displacement_definition:

Tangential displacement
-------------------------------------

History-dependent tangential models (*linear_history\**, *mindlin\**) 
rely on accumulated tangential displacement (i.e. contact history), 
except for the options *mindlin/force* and *mindlin_rescale/force*, 
which use accumulated tangential force instead, and are discussed further below.
The tangential displacement accumulated during the entire duration of the contact
is denoted :math:`\mathbf{\xi}`, and is given by:

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

Furthermore, when the tangential force exceeds a critical force 
:math:`F_{t,crit}` (see :ref:`<critical_force_definition>`),
the magnitude of the tangential force is truncated to the critical force, i.e., 
:math:`F_t = F_{t,crit}`. Additionally, the tangential displacement is 
re-scaled to match the value for the critical force 
(see :ref:`Luding <Luding2008>`, equation 20 and related
discussion):

.. math::

   \mathbf{\xi} = -\frac{1}{k_t}\left(\mu_t F_{n0}\mathbf{t} - \mathbf{F}_{t,damp}\right)


.. _critical_force_definition:

Critical force
----------------------------

All friction models use a critical force to determine if slip occurs,
generically calculated as:

.. math::

   F_{t,crit} = \mu_t F_{n0}

Here, :math:`mu_t` is the tangential friction coefficient.
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



.. _linear_nohistory_tangential_model:

*linear_nohistory* tangential model
------------------------------------

For the *linear_nohistory* tangential model, a simple velocity-dependent Coulomb
friction criterion is used, which mimics the behavior of the *pair
gran/hooke* style. The tangential force :math:`\mathbf{F}_t` is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|\mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}

.. _linear_history_tangential_model:

*linear_history* tangential model
----------------------------------

For *tangential linear_history*, the tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t\mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}


.. _mindlin_tangential_model:

*mindlin* tangential model
--------------------------------


For *tangential mindlin*, the :ref:`Mindlin <Mindlin1949>` no-slip solution
is used which differs from the *linear_history* option by an additional factor
of :math:`a`, the radius of the contact region. The tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t a \mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}


Here, :math:`a` is the radius of the contact region, given by :math:`a =\sqrt{R\delta}`
for all normal contact models, except for *jkr*, where it is given
implicitly by :math:`\delta = a^2/R - 2\sqrt{\pi \gamma a/E}`, see
:ref:`jkr <jkr_normal_model>`. To match the Mindlin solution, one should set
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


.. _mindlin_force_tangential_model:

*mindlin/force* tangential model
--------------------------------

For *tangential mindlin/force*, the accumulated elastic tangential force
characterizes the contact history, instead of the accumulated tangential
displacement. This prevents the dependence of the tangential force on the
normal overlap as noted for the :ref:`mindlin <mindlin_tangential_model>`
model. The tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|\mathbf{F}_{te} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}

The increment of the elastic component of the tangential force
:math:`\mathbf{F}_{te}` is given by:

.. math::

   \mathrm{d}\mathbf{F}_{te} = -k_t a \mathbf{v}_{t,rel} \mathrm{d}\tau

The changes in frame of reference of the contacting pair of particles during
contact are accounted for by the same formula as :ref:`above <tangential_displacement_definition>`
, replacing the accumulated tangential displacement :math:`\xi`, by the accumulated tangential
elastic force :math:`F_{te}`. When the tangential force exceeds the critical
force, the tangential force is directly re-scaled to match the value for
the critical force:

.. math::

   \mathbf{F}_{te} = - \mu_t F_{n0}\mathbf{t} + \mathbf{F}_{t,damp}

The same rules as those described for *mindlin* apply regarding the tangential
stiffness and mixing of the shear modulus for different particle types.

.. _mindlin_rescale_tangential_model:

*mindlin/rescale* tangential model
----------------------------------

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

.. _mindlin_rescale_force_tangential_model:

*mindlin_rescale/force* tangential model
-----------------------------------------

The *mindlin_rescale/force* option uses the same form as :ref:`mindlin/force<mindlin_force_tangential_model>`,
but the magnitude of the tangential elastic force is re-scaled as the contact
unloads, i.e. if :math:`a < a_{t_{n-1}}`:

.. math::

   \mathbf{F}_{te} = \mathbf{F}_{te, t_{n-1}} \frac{a}{a_{t_{n-1}}}

This approach provides a better approximation of the :ref:`Mindlin-Deresiewicz <Mindlin1953>`
laws and is more consistent than *mindlin_rescale*. See discussions in
:ref:`Thornton et al, 2013 <Thornton2013>`, particularly equation 18(b) of that
work and associated discussion, and :ref:`Agnolin and Roux, 2007 <AgnolinRoux2007>`,
particularly Appendix A.


.. _linear_history_static_tangential_model:

*linear_history/static* tangential model
-----------------------------------------

The *linear_history/static* option mirrors the :ref:`linear_history <linear_history_tangential_model>`
model, but allows different static and dynamic friction coefficients, 
:math:`\mu_{t,static}` and :math:`\mu_{t,dynamic}`. This allows simulation
of distinct static and dynamic friction.

The static and dynamic critical friction forces are defined as:

.. math::

   F_{t,crit,static} = \mu_{t,static} F_{n0}


and 

.. math::

   F_{t,crit,dynamic} = \mu_{t,dynamic} F_{n0}


The state of a contact is designated as *static* or *dynamic*, and tracked
throughout a simulation. The initial state of a contact is initialized to 
*static* at the start of a simulation (or read from a restart file, similar 
to other history entries). 

While the contact is designated as *static*, the critical force for purposes
of determining slip and truncating the tangential displacement 
(see :ref:`tangential_displacement_definition`) is given by 
:math:`F_{t,crit} = F_{t,crit,static}`.
If the tangential force exceeds the static critical force while the contact
is *static*, i.e. :math:`F_t > F_{t,crit,static}`, the contact state is changed
to *dynamic*. 

Conversely, while the contact is designated as *dynamics*, the critical force 
is given by :math:`F_{t,crit} = F_{t,crit,dynamic}`.
If the tangential force drops below the dynamic critical force while the contact
is *dynamic*, i.e. :math:`F_t < F_{t,crit,dynamic}`, the contact state is changed
back to *static*. 

.. _mindlin_static_tangential_model:

*mindlin/static* tangential model
-----------------------------------------

The *mindlin/static* option mirrors the :ref:`mindlin <mindlin_tangential_model>`
model, but allows different static and dynamic friction coefficients, 
:math:`\mu_{t,static}` and :math:`\mu_{t,dynamic}`. The 
determination of the critical force and designation of a contact as *static*
or *dynamic* are the same as described for 
:ref:`linear_history/static <linear_history_static_tangential_model>`.

References
""""""""""


.. _WaltonPC:

**(Otis R. Walton)** Walton, O.R., Personal Communication

.. _Mindlin1949:

**(Mindlin, 1949)** Mindlin, R. D. (1949). Compliance of elastic bodies
in contact.  J. Appl. Mech., ASME 16, 259-268.

.. _Mindlin1953:

**(Mindlin and Deresiewicz, 1953)** Mindlin, R.D., & Deresiewicz, H (1953).
Elastic Spheres in Contact under Varying Oblique Force.
J. Appl. Mech., ASME 20, 327-344.

.. _AgnolinRoux2007:

**(Agnolin and Roux 2007)** Agnolin, I. & Roux, J-N. (2007).
Internal states of model isotropic granular packings.
I. Assembling process, geometry, and contact networks. Phys. Rev. E, 76, 061302.

.. _Thornton2013:

**(Thornton et al, 2013)** Thornton, C., Cummins, S. J., & Cleary,
P. W. (2013).  An investigation of the comparative behavior of
alternative contact force models during inelastic collisions. Powder
Technology, 233, 30-46.