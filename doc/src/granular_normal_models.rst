Models for normal contact in granular interactions
==================================================

The normal force acts along the vector connecting the center of two particles,
i.e. normal to the plane of contact between particles.
In all cases, the normal forces is modeled as a function of overlap :math:`delta`. 
The following quantities are common to all the normal models currently available:

:math:`\delta_{ij} = R_i + R_j - \|\mathbf{r}_{ij}\|` is the particle overlap, 
:math:`R_i, R_j` are the particle radii, 
:math:`\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j` is the vector separating the two
particle centers (note the i-j ordering so that the force is
positive for repulsion), and 
:math:`\mathbf{n} = \frac{\mathbf{r}_{ij}}{\|\mathbf{r}_{ij}\|}`

.. _hooke_normal_model:

`hooke` model
---------------

For the *hooke* model, the normal component of force acting
on particle *i* due to contact with particle *j* is given by:

.. math::

   \mathbf{F}_{ne, Hooke} = k_n \delta_{ij} \mathbf{n}

The units of the spring constant :math:`k_n` are
*force*\ /\ *distance*, or equivalently *mass*\ /*time\^2*.

.. _hertz_normal_model:

`hertz` model
-------------

For the *hertz* model, the normal force is given by:

.. math::

   \mathbf{F}_{ne, Hertz} = k_n R_{eff}^{1/2}\delta_{ij}^{3/2} \mathbf{n}

The units of the spring constant :math:`k_n` are *force*\ /\ *length*\ \^2, or
equivalently *pressure*\ .

.. _hertz_material_normal_model:

`hertz/material` model
----------------------

For the *hertz/material* model, the normal force is given by:

.. math::

   \mathbf{F}_{ne, Hertz/material} = \frac{4}{3} E_{eff} R_{eff}^{1/2}\delta_{ij}^{3/2} \mathbf{n}

Here, :math:`E_{eff} = E = \left(\frac{1-\nu_i^2}{E_i} + \frac{1-\nu_j^2}{E_j}\right)^{-1}` is the effective Young's
modulus, with :math:`\nu_i, \nu_j` the Poisson ratios of the particles of
types *i* and *j*\ . Note that if the elastic modulus and the shear
modulus of the two particles are the same, the *hertz/material* model
is equivalent to the *hertz* model with :math:`k_n = 4/3 E_{eff}`

.. _dmt_normal_model:

`dmt` model
----------------------

The *dmt* model corresponds to the
:ref:`(Derjaguin-Muller-Toporov) <DMT1975>` cohesive model, where the force
is simply Hertz with an additional attractive cohesion term:

.. math::

   \mathbf{F}_{ne, dmt} = \left(\frac{4}{3} E R^{1/2}\delta_{ij}^{3/2} - 4\pi\gamma R\right)\mathbf{n}

.. _jkr_normal_model:

`jkr` model
----------------------

The *jkr* model is the :ref:`(Johnson-Kendall-Roberts) <JKR1971>` model,
where the force is computed as:

.. math::

   \mathbf{F}_{ne, jkr} = \left(\frac{4Ea^3}{3R} - 2\pi a^2\sqrt{\frac{4\gamma E}{\pi a}}\right)\mathbf{n}

Here, :math:`a` is the radius of the contact zone, related to the overlap
:math:`\delta` according to:

.. math::

   \delta = a^2/R - 2\sqrt{\pi \gamma a/E}

LAMMPS internally inverts the equation above to solve for *a* in terms
of :math:`\delta`, then solves for the force in the previous
equation. Additionally, note that the JKR model allows for a tensile
force beyond contact (i.e. for :math:`\delta < 0`), up to a maximum of
:math:`3\pi\gamma R` (also known as the 'pull-off' force).  Note that this
is a hysteretic effect, where particles that are not contacting
initially will not experience force until they come into contact
:math:`\delta \geq 0`; as they move apart and (:math:`\delta < 0`), they
experience a tensile force up to :math:`3\pi\gamma R`, at which point they
lose contact.


.. _epa_linear_normal_model:

`epa_linear` model
----------------------

The *epa_linear* model is the linear elastic-plastic-adhesive model proposed
by :ref:`(Luding) <Luding2008>` model, where the force is computed according to:

.. math::

   \mathbf{F}_{ne, luding} = \left(\frac{4Ea^3}{3R} - 2\pi a^2\sqrt{\frac{4\gamma E}{\pi a}}\right)\mathbf{n}

Here, :math:`a` is the radius of the contact zone, related to the overlap
:math:`\delta` according to:

.. math::

   \delta = a^2/R - 2\sqrt{\pi \gamma a/E}


.. _epa_nonlinear_normal_model:

`epa_nonlinear` model
----------------------

The *epa_nonlinear* model is very similar to the nonlinear elastic-plastic-adhesive model proposed
by :ref:`(Thankur et al) <Thakur2014>` model, where the force is computed according to:

.. math::

   \mathbf{F}_{ne, luding} = \left(\frac{4Ea^3}{3R} - 2\pi a^2\sqrt{\frac{4\gamma E}{\pi a}}\right)\mathbf{n}

Here, :math:`a` is the radius of the contact zone, related to the overlap
:math:`\delta` according to:

.. math::

   \delta = a^2/R - 2\sqrt{\pi \gamma a/E}

References
""""""""""

.. _JKR1971:

**(Johnson et al, 1971)** Johnson, K. L., Kendall, K., & Roberts,
A. D. (1971).  Surface energy and the contact of elastic
solids. Proc. R. Soc. Lond. A, 324(1558), 301-313.

.. _DMT1975:

**(Derjaguin et al, 1975)** Derjaguin, B. V., Muller, V. M., & Toporov,
Y. P. (1975). Effect of contact deformations on the adhesion of
particles. Journal of Colloid and interface science, 53(2), 314-326.

.. _Luding2008:

**(Luding, 2008)** Luding, S. (2008). Cohesive, frictional powders:
contact models for tension. Granular matter, 10(4), 235.

**(Thakur et al, 2014)** Thakur, Subhash C., et al. (2014). 
Micromechanical analysis of cohesive granular materials using 
the discrete element method with an adhesive  elasto-plastic contact 
model. Granular Matter 16, 383-400.