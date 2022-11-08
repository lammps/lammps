.. index:: fix brownian
.. index:: fix brownian/sphere
.. index:: fix brownian/asphere

fix brownian command
===========================

fix brownian/sphere command
===========================

fix brownian/asphere command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style_name temp seed keyword args

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *brownian* or *brownian/sphere* or *brownian/asphere*
* temp = temperature
* seed = random number generator seed
* one or more keyword/value pairs may be appended
* keyword = *rng* or *dipole* or *gamma_r_eigen* or *gamma_t_eigen* or *gamma_r* or *gamma_t* or *rotation_temp* or *planar_rotation*

  .. parsed-literal::

     *rng* value = *uniform* or *gaussian* or *none*
       *uniform* = use uniform random number generator
       *gaussian* = use gaussian random number generator
       *none* = turn off noise
     *dipole* value = *mux* and *muy* and *muz* for *brownian/asphere*
       *mux*, *muy*, and *muz* = update orientation of dipole having direction (*mux*,*muy*,*muz*) in body frame of rigid body
     *gamma_r_eigen* values = *gr1* and *gr2* and *gr3* for *brownian/asphere*
       *gr1*, *gr2*, and *gr3* = diagonal entries of body frame rotational friction tensor
     *gamma_r* values = *gr* for *brownian/sphere*
       *gr* = magnitude of the (isotropic) rotational friction tensor
     *gamma_t_eigen* values = *gt1* and *gt2* and *gt3* for *brownian/asphere*
       *gt1*, *gt2*, and *gt3* = diagonal entries of body frame translational friction tensor
     *gamma_t* values = *gt* for *brownian* and *brownian/sphere*
        *gt* = magnitude of the (isotropic) translational friction tensor
     *rotation_temp* values = *T* for *brownian/sphere* and *brownian/asphere*
        *T* = rotation temperature, which can be different then *temp* when out of equilibrium
     *planar_rotation* values = none (constrains rotational diffusion to be in xy plane if in 3D)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all brownian 1.0 12908410 gamma_t 1.0
   fix 1 all brownian 1.0 12908410 gamma_t 3.0 rng gaussian
   fix 1 all brownian/sphere 1.0 1294019 gamma_t 3.0 gamma_r 1.0
   fix 1 all brownian/sphere 1.0 19581092 gamma_t 1.0 gamma_r 0.3  rng none
   fix 1 all brownian/asphere 1.0 1294019 gamma_t_eigen 1.0 2.0 3.0 gamma_r_eigen 4.0 7.0 8.0 rng gaussian
   fix 1 all brownian/asphere 1.0 1294019 gamma_t_eigen 1.0 2.0 3.0 gamma_r_eigen 4.0 7.0 8.0 dipole 1.0 0.0 0.0


Description
"""""""""""

Perform Brownian Dynamics time integration to update position, velocity,
dipole orientation (for spheres) and quaternion orientation (for
ellipsoids, with optional dipole update as well) of all particles in the
fix group in each timestep.  Brownian Dynamics uses Newton's laws of
motion in the limit that inertial forces are negligible compared to
viscous forces. The stochastic equation of motion for the center of mass
positions is

.. math::

   d\mathbf{r} = \boldsymbol{\gamma}_t^{-1}\mathbf{F}dt
   + \sqrt{2k_B T}\boldsymbol{\gamma}_t^{-1/2}d\mathbf{W}_t,

in the lab-frame (i.e., :math:`\boldsymbol{\gamma}_t` is not diagonal, but
only depends on orientation and so the noise is still additive).

The rotational motion for the spherical and ellipsoidal particles is not
as simple an expression, but is chosen to replicate the Boltzmann
distribution for the case of conservative torques (see :ref:`(Ilie)
<Ilie1>` or :ref:`(Delong) <Delong1>`).

For the style *brownian*, only the positions of the particles are
updated. This is therefore suitable for point particle simulations.

For the style *brownian/sphere*, the positions of the particles are
updated, and a dipole slaved to the spherical orientation is also
updated. This style therefore requires the hybrid atom style
:doc:`atom_style dipole <atom_style>` and :doc:`atom_style sphere
<atom_style>`. The equation of motion for the dipole is

.. math::

   \boldsymbol{\mu}(t+dt) = \frac{\boldsymbol{\mu}(t) + \boldsymbol{\omega} \times \boldsymbol{\mu}dt
   }{|\boldsymbol{\mu}(t) + \boldsymbol{\omega} \times \boldsymbol{\mu}|}

which correctly reproduces a Boltzmann distribution of orientations and
rotational diffusion moments (see :ref:`(Ilie) <Ilie1>`) when

.. math::

   \boldsymbol{\omega} = \frac{\mathbf{T}}{\gamma_r} + \sqrt{\frac{2 k_B T_{rot}}{\gamma_r}\frac{d\mathbf{W}}{dt}},

with :math:`d\mathbf{W}` being a random number with zero mean and variance :math:`dt`
and :math:`T_{rot}` is *rotation_temp*.

For the style *brownian/asphere*, the center of mass positions and the
quaternions of ellipsoidal particles are updated. This fix style is
suitable for equations of motion where the rotational and translational
friction tensors can be diagonalized in a certain (body) reference
frame. In this case, the rotational equation of motion is updated via
the quaternion

.. math::

   \mathbf{q}(t+dt) = \frac{\mathbf{q}(t) + d\mathbf{q}}{\lVert\mathbf{q}(t) + d\mathbf{q}\rVert}

which correctly reproduces a Boltzmann distribution of orientations and rotational
diffusion moments [see :ref:`(Ilie) <Ilie1>`] when the quaternion step is given by

.. math::

   d\mathbf{q} = \boldsymbol{\Psi}\boldsymbol{\omega}dt

where :math:`\boldsymbol{\Psi}` has rows :math:`(-q_1,-q_2,-q_3)`,
:math:`(q_0,-q_3,q_2)`, :math:`(q_3,q_0,-q_1)`, and :math:`(-q_2,q_1,q_0)`.
:math:`\boldsymbol{\omega}` is evaluated in the body frame of reference where the
friction tensor is diagonal.  See :ref:`(Delong) <Delong1>` for more details of
a similar algorithm.


---------

.. note::

   This integrator does not by default assume a relationship between the
   rotational and translational friction tensors, though such a
   relationship should exist in the case of no-slip boundary conditions
   between the particles and the surrounding (implicit) solvent. For example,
   in the case of spherical particles, the condition
   :math:`\gamma_t=3\gamma_r/\sigma^2` must be explicitly accounted for
   by setting *gamma_t* to 3x and *gamma_r* to x (where :math:`\sigma`
   is the sphere's diameter). A similar (though more complex)
   relationship holds for ellipsoids and rod-like particles. The
   translational diffusion and rotational diffusion are given by
   *temp/gamma_t* and *rotation_temp/gamma_r*.

---------

.. note::

   Temperature computation using the :doc:`compute temp <compute_temp>`
   will not correctly compute the temperature of these overdamped dynamics
   since we are explicitly neglecting inertial effects.  Furthermore,
   this time integrator does not add the stochastic terms or viscous
   terms to the force and/or torques.  Rather, they are just added in to
   the equations of motion to update the degrees of freedom.

---------


If the *rng* keyword is used with the *uniform* value, then the noise
is generated from a uniform distribution (see
:ref:`(Dunweg) <Dunweg7>` for why this works). This is the same method
of noise generation as used in :doc:`fix_langevin <fix_langevin>`.

If the *rng* keyword is used with the *gaussian* value, then the noise
is generated from a Gaussian distribution. Typically this added
complexity is unnecessary, and one should be fine using the *uniform*
value for reasons argued in :ref:`(Dunweg) <Dunweg7>`.

If the *rng* keyword is used with the *none* value, then the noise
terms are set to zero.

The *gamma_t* keyword sets the (isotropic) translational viscous damping.
Required for (and only compatible with) *brownian* and *brownian/sphere*.
The units of *gamma_t* are mass/time.

The *gamma_r* keyword sets the (isotropic) rotational viscous damping.
Required for (and only compatible with) *brownian/sphere*.
The units of *gamma_r* are mass*length**2/time.

The *gamma_r_eigen*, and *gamma_t_eigen* keywords are the eigenvalues of
the rotational and viscous damping tensors (having the same units as
their isotropic counterparts). Required for (and only compatible with)
*brownian/asphere*. For a 2D system, the first two values of
*gamma_r_eigen* must be *inf* (only rotation in *x*\ --\ *y* plane), and the third
value of *gamma_t_eigen* must be *inf* (only diffusion in the *x*\ --\ *y* plane).

If the *dipole* keyword is used, then the dipole moments of the particles
are updated as described above. Only compatible with *brownian/asphere*
(as *brownian/sphere* updates dipoles automatically).

If the *rotation_temp* keyword is used, then the rotational diffusion
will be occur at this prescribed temperature instead of *temp*. Only
compatible with *brownian/sphere* and *brownian/asphere*.

If the *planar_rotation* keyword is used, then rotation is constrained
to the *x*\ -- *y* plane in a 3D simulation. Only compatible with
*brownian/sphere* and *brownian/asphere* in 3D.

----------

.. note::
   For style *brownian/asphere*, the components *gamma_t_eigen* = (x,x,x) and
   *gamma_r_eigen* = (y,y,y), the dynamics will replicate those of the
   *brownian/sphere* style with *gamma_t* = x and *gamma_r* = y.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The style *brownian/sphere* fix requires that atoms store torque and
angular velocity (omega) as defined by the :doc:`atom_style sphere
<atom_style>` command.  The style *brownian/asphere* fix requires that
atoms store torque and quaternions as defined by the :doc:`atom_style
ellipsoid <atom_style>` command.  If the *dipole* keyword is used, they
must also store a dipole moment as defined by the :doc:`atom_style
dipole <atom_style>` command.

This fix is part of the BROWNIAN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix propel/self <fix_propel_self>`,
:doc:`fix langevin <fix_langevin>`, :doc:`fix nve/sphere <fix_nve_sphere>`,

Default
"""""""

The default for *rng* is *uniform*. The default for the rotational and
translational friction tensors are the identity tensor.

----------

.. _Ilie1:

**(Ilie)** Ilie, Briels, den Otter, Journal of Chemical Physics, 142, 114103 (2015).

.. _Delong1:

**(Delong)** Delong, Usabiaga, Donev, Journal of Chemical Physics. 143, 144107 (2015)

.. _Dunweg7:

**(Dunweg)** Dunweg and Paul, Int J of Modern Physics C, 2, 817-27 (1991).
