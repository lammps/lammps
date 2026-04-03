.. index:: pair_style granular/superellipsoid

pair_style granular/superellipsoid command
==========================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style granular/superellipsoid cutoff no_bounding_box curvature_gaussian

   Optional settings, see discussion below.
   * cutoff = global cutoff value
   * no_bounding_box = skip oriented bounding box check
   * curvature_gaussian = gaussian curvature coeff approximation for contact patch

Examples
""""""""

.. code-block:: LAMMPS

   pair_style granular/superellipsoid
   pair_coeff * * hooke 1000.0 50.0 tangential linear_history 1000.0 1.0 0.5 damping mass_velocity

   pair_style granular/superellipsoid 10.0 curvature_gaussian
   pair_coeff 1 1 hertz 1000.0 50.0 tangential linear_history 500.0 1.0 0.4 damping viscoelastic
   pair_coeff 2 2 hertz 500.0 50.0 tangential linear_history 250.0 1.0 0.1 damping viscoelastic

Description
"""""""""""

.. versionadded:: 30Mar2026

The *granular/superellipsoid* style calculates granular contact forces
between superellipsoidal particles (see :doc:`atom style ellipsoid
<atom_style>`).  Similar to the :doc:`granular pairstyle
<pair_granular>` which is designed for spherical particles, various
normal, damping, and tangential contact models are available (rolling
and twisting may be added later).  The total computed forces and torques
are the sum of various models selected.

All model choices and parameters are entered in the :doc:`pair_coeff
<pair_coeff>` command, as described below.  Coefficient values are not
global, but can be set to different values for different combinations of
particle types, as determined by the :doc:`pair_coeff <pair_coeff>`
command.  If the contact model choice is the same for two particle
types, the mixing for the cross-coefficients can be carried out
automatically.  This is shown in the last example, where model choices
are the same for type 1 - type 1 as for type 2 - type2 interactions, but
coefficients are different.  In this case, the mixed coefficients for
type 1 - type 2 interactions can be determined from mixing rules
discussed below.  For additional flexibility, coefficients as well as
model forms can vary between particle types.

----------

This pair_style allows granular contact between two superellipsoid
particles whose surface is implicitly defined as:

.. math::

    f(\mathbf{x}) = \left(
    \left|\frac{x}{a}\right|^{n_2} + \left|\frac{y}{b}\right|^{n_2}
    \right)^{n_1 / n_2}
    + \left|\frac{z}{c}\right|^{n_1} - 1 = 0

for a point :math:`\mathbf{x} = (x, y, z)` where the coordinates are
given in the reference of the principal directions of inertia of the
particle.  The half-diameters :math:`a`, :math:`b`, and :math:`c`
correspond to the *shape* property, and the exponents :math:`n_1` and
:math:`n_2` to the *block* property of the ellipsoid atom.  See the doc
page for the :doc:`set <set>` command for more details.

.. note::

    The contact solver strictly requires convex particle shapes to
    ensure a mathematically unique point of deepest penetration.
    Therefore, the blockiness parameters must be :math:`n_1 \ge 2.0` and
    :math:`n_2 \ge 2.0`.  Attempting to simulate concave or "pointy"
    particles (:math:`n < 2.0`) will result in an error.

.. note::

    For particles with high blockiness exponents (:math:`n > 4.0`)
    involved in edge-to-edge or corner-to-corner contacts, the surface
    normal vector varies rapidly over small distances.  The Newton
    solver may occasionally fail to converge to the strict gradient
    alignment tolerance (typically :math:`10^{-10}`).  You may see
    warning messages in the log indicating that the solver returned a
    sub-optimal solution, but the simulation will proceed using this
    best-effort contact point.

Contact detection for these aspherical particles uses the so-called
"midway" minimization approach from :ref:`(Houlsby) <Houlsby>`.
Considering two particles with shape functions, :math:`F_i` and
:math:`F_j`, the contact point :math:`\mathbf{X}_0` in the global frame
is obtained as:

.. math::

    \mathbf{X}_0 = \underset{\mathbf{X}}{\text{argmin}}
                   \ F_i(\mathbf{X}) + F_j(\mathbf{X})
                   \text{, subject to } F_i(\mathbf{X}) = F_j(\mathbf{X})

where the shape function is given by :math:`F_i(\mathbf{X}) =
f_i(\mathbf{R}_i^T (\mathbf{X} - \mathbf{X}_i))` and where
:math:`\mathbf{X}_i` and :math:`\mathbf{R}_i` are the center of mass and
rotation matrix of the particle, respectively.  The constrained
minimization problem is solved using Lagrange multipliers and Newton's
method with a line search as described by :ref:`(Podlozhnyuk)
<Podlozhnyuk>`.

.. note::

    The shape function :math:`F` is not a signed distance function and
    does not have unit gradient :math:`\|\nabla F \| \neq 1` so that the
    so-called "midway" point is not actually located at an equal
    distance from the surface of both particles.  For contact between
    non-identical particles, the contact point tends to be closer to the
    surface of the smaller and blockier particle.

.. note::

    This formulation leads to a 4x4 system of non-linear equations.
    Tikhonov regularization and step clumping is used to ensure
    robustness of the direct solver and high convergence rate, even for
    blocky particles with near flat faces.

The particles overlap if both shape functions are negative at the
contact point.  The contact normal is obtained as:
:math:`\mathbf{n}_{ij} = \nabla F_i(\mathbf{X}_0) / \| \nabla
F_i(\mathbf{X}_0)\| = - \nabla F_j(\mathbf{X}_0) / \| \nabla
F_j(\mathbf{X}_0)\|` and the overlap :math:`\delta =
\|\mathbf{X}_j^{\mathrm{surf}} - \mathbf{X}_i^{\mathrm{surf}}\|` is
computed as the distance between the points on the particles surfaces
that are closest to the contact point in the direction of the contact
normal: :math:`F_i(\mathbf{X}_i^{\mathrm{surf}} = \mathbf{X}_0 +
\lambda_i \mathbf{n}_{ij}) = 0` and
:math:`F_j(\mathbf{X}_j^{\mathrm{surf}} = \mathbf{X}_0 + \lambda_j
\mathbf{n}_{ij}) = 0`.  Newton's method is used to solve this equation
for the scalars :math:`\lambda_i` and :math:`\lambda_j` and find the
surface points :math:`\mathbf{X}_i^{\mathrm{surf}}` and
:math:`\mathbf{X}_j^{\mathrm{surf}}`.

.. note::

   A modified representation of the particle surface is defined as
   :math:`G(\mathbf{X}) = (F(\mathbf{X})+1)^{1/n_1}-1` which is a
   quasi-radial distance function formulation.  This formulation is used
   to compute the surface points once the "midway" contact point is
   found.  This formulation is also used when the *geometric* keyword is
   specified in the pair_style command and the following optimization
   problem is solved instead for the contact point: :math:`\mathbf{X}_0
   = \underset{\mathbf{X}}{\text{argmin}} \, \left( r_i
   G_i(\mathbf{X}) + r_j G_j(\mathbf{X}) \right) \text{, subject to }
   r_i G_i(\mathbf{X}) = r_j G_j(\mathbf{X})`, where :math:`r_i` and
   :math:`r_j` are the average radii of the two particles.  The
   geometric formulation thus yields a better approximation of the
   contact point for particles with different sizes, and it is slightly
   more robust for particles with high *block* exponents, albeit more
   computationally expensive.

A hierarchical approach is used to limit the cost of contact detection.
First, intersection of the bounding spheres of the two particles of
bounding radii :math:`r_i` and :math:`r_j` is checked.  If the distance
between the particles center is more than the sum of the radii
:math:`\|\mathbf{X}_j - \mathbf{X}_j\| > r_i + r_j`, the particles do
not intersect.  Then, if the bounding spheres intersect, intersection of
the oriented bounding box is checked.  This is done following the
equations of :ref:`(Eberly) <GeometricTools>`.  This check is always
performed, unless the *no_bounding_box* keyword is used.  This is
advantageous for all particle shapes except for superellipses with
aspect ratio close to one and both blockiness indexes close to 2.

.. warning::

    The Newton-Raphson minimization used to find the midway contact
    point can fail to converge if the initial starting guess is too far
    from the true physical surface.  This typically occurs if a user
    specifies a manual global *cutoff* that is significantly larger than
    the particles **and** enables the *no_bounding_box* keyword.  Under
    these conditions, the solver attempts to resolve contacts between
    widely separated particles, which might cause the math to diverge
    and instantly crashing the simulation.  It is strongly recommended
    to keep bounding box checks enabled if a large cutoff is specified.

----------

This section provides an overview of the various normal, tangential, and
damping contact models available.  For additional context, see the
discussion in the :doc:`granular pairstyle <pair_granular>` doc page
which includes all of these options.

The first required keyword for the *pair_coeff* command is the normal
contact model.  Currently supported options for normal contact models
and their required arguments are:

1. *hooke* : :math:`k_n`, :math:`\eta_{n0}` (or :math:`e`)
2. *hertz* : :math:`k_n`, :math:`\eta_{n0}` (or :math:`e`)

Here, :math:`k_n` is spring stiffness (with units that depend on model
choice, see below); :math:`\eta_{n0}` is a damping prefactor (or, in its
place a coefficient of restitution :math:`e`, depending on the choice of
damping mode, see below).

For the *hooke* model, the normal, elastic component of force acting on
particle *i* due to contact with particle *j* is given by:

.. math::

   \mathbf{F}_{ne, Hooke} = k_n \delta_{ij} \mathbf{n}

Where :math:`\delta_{ij}` is the particle overlap, (note the i-j
ordering so that :math:`\mathbf{F}_{ne}` is positive for repulsion), and
:math:`\mathbf{n}` is the contact normal vector at the contact point.
Therefore, for *hooke*, the units of the spring constant :math:`k_n` are
*force*\ /\ *distance*, or equivalently *mass*\ /*time\^2*.

For the *hertz* model, the normal component of force is given by:

.. math::

   \mathbf{F}_{ne, Hertz} = k_n R_{eff}^{1/2}\delta_{ij}^{3/2} \mathbf{n}

Here, :math:`R_{eff} = R = \frac{R_i R_j}{R_i + R_j}` is the effective
radius, and :math:`R_i` is the equivalent radius of the i-th particle at
the surface contact point with the j-th particle.  This radius is either
the inverse of the mean curvature coefficient, :math:`R_i = 2 /
(\kappa_1 + \kappa_2)`, or the gaussian curvature coefficient :math:`R_i
= 1 / \sqrt{\kappa_1 \kappa_2}`, where :math:`\kappa_{1,2}` are the
principal curvatures of the particle surface at the contact point.  For
*hertz*, the units of the spring constant :math:`k_n` are *force*\ /\
*length*\ \^2, or equivalently *pressure*\ .

.. note::

    To ensure numerical stability and preserve physical realism, the
    computed contact radius is mathematically capped.  For highly blocky
    particles undergoing flat-on-flat contact, the theoretical curvature
    approaches zero, which would yield an infinite contact radius and
    cause a force explosion.  To prevent this, the maximum contact
    radius is capped at the physical bounding radius of the smallest
    interacting particle.  Conversely, for sharp corner contacts where
    curvature approaches infinity, the calculated radius would drop to
    zero, eliminating the repulsive force entirely.  The contact radius
    is therefore lower-bounded by a minimum fraction of the physical
    radius (:math:`10^{-4} \min(r_i, r_j)`) to prevent particles from
    unphysically interpenetrating.

In addition, the normal force is augmented by a damping term of the
following general form:

.. math::

   \mathbf{F}_{n,damp} = -\eta_n \mathbf{v}_{n,rel}

Here, :math:`\mathbf{v}_{n,rel} = (\mathbf{v}_j - \mathbf{v}_i) \cdot
\mathbf{n}\ \mathbf{n}` is the component of relative velocity along
:math:`\mathbf{n}`.

The optional *damping* keyword to the *pair_coeff* command followed by a
keyword determines the model form of the damping factor :math:`\eta_n`,
and the interpretation of the :math:`\eta_{n0}` or :math:`e`
coefficients specified as part of the normal contact model settings.
The *damping* keyword and corresponding model form selection may be
appended anywhere in the *pair coeff* command.  Note that the choice of
damping model affects both the normal and tangential damping.  The
options for the damping model currently supported are:

1. *mass_velocity*
2. *viscoelastic*

If the *damping* keyword is not specified, the *viscoelastic* model is
used by default.

For *damping mass_velocity*, the normal damping is given by:

.. math::

   \eta_n = \eta_{n0} m_{eff}

Here, :math:`\eta_{n0}` is the damping coefficient specified for the
normal contact model, in units of 1/\ *time* and :math:`m_{eff} = m_i
m_j/(m_i + m_j)` is the effective mass.  Use *damping mass_velocity* to
reproduce the damping behavior of *pair gran/hooke/\**.

The *damping viscoelastic* model is based on the viscoelastic treatment
of :ref:`(Brilliantov et al) <Brill1996_2>`, where the normal damping is
given by:

.. math::

   \eta_n = \eta_{n0}\ a m_{eff}

Here, *a* is the contact radius, given by :math:`a =\sqrt{R\delta}` for
all models.  For *damping viscoelastic*, :math:`\eta_{n0}` is in units
of 1/(\ *time*\ \*\ *distance*\ ).

The total normal force is computed as the sum of the elastic and damping
components:

.. math::

   \mathbf{F}_n = \mathbf{F}_{ne} + \mathbf{F}_{n,damp}

----------

The *pair_coeff* command also requires specification of the tangential
contact model.  The required keyword *tangential* is expected, followed
by the model choice and associated parameters.  Currently there is only
one supported tangential model with expected parameters as follows:

1. *linear_history* : :math:`k_t`, :math:`x_{\gamma,t}`, :math:`\mu_s`

Here, :math:`x_{\gamma,t}` is a dimensionless multiplier for the normal
damping :math:`\eta_n` that determines the magnitude of the tangential
damping, :math:`\mu_t` is the tangential (or sliding) friction
coefficient, and :math:`k_t` is the tangential stiffness coefficient.

The tangential damping force :math:`\mathbf{F}_\mathrm{t,damp}` is given
by:

.. math::

   \mathbf{F}_\mathrm{t,damp} = -\eta_t \mathbf{v}_{t,rel}

The tangential damping prefactor :math:`\eta_t` is calculated by scaling
the normal damping :math:`\eta_n` (see above):

.. math::

   \eta_t = -x_{\gamma,t} \eta_n

The normal damping prefactor :math:`\eta_n` is determined by the choice
of the *damping* keyword, as discussed above.  Thus, the *damping*
keyword also affects the tangential damping.  The parameter
:math:`x_{\gamma,t}` is a scaling coefficient.  Several works in the
literature use :math:`x_{\gamma,t} = 1` (:ref:`Marshall <Marshall2009_2>`,
:ref:`Tsuji et al <Tsuji1992_2>`, :ref:`Silbert et al <Silbert2001_2>`).
The relative tangential velocity at the point of contact is given by
:math:`\mathbf{v}_{t, rel} = \mathbf{v}_{t} -
(R_i\boldsymbol{\Omega}_i + R_j\boldsymbol{\Omega}_j) \times
\mathbf{n}`, where :math:`\mathbf{v}_{t} = \mathbf{v}_r -
\mathbf{v}_r\cdot\mathbf{n}\ \mathbf{n}`, :math:`\mathbf{v}_r =
\mathbf{v}_j - \mathbf{v}_i` .  The direction of the applied force is
:math:`\mathbf{t} = \mathbf{v_{t,rel}}/\|\mathbf{v_{t,rel}}\|` .

The normal force value :math:`F_{n0}` used to compute the critical force
depends on the form of the contact model.  It is given by the magnitude
of the normal force:

.. math::

   F_{n0} = \|\mathbf{F}_n\|

The remaining tangential options all use accumulated tangential
displacement (i.e. contact history).  The accumulated tangential
displacement is discussed in details below in the context of the
*linear_history* option.  The same treatment of the accumulated
displacement will apply to other (future) options as well.

For *tangential linear_history*, the tangential force is given by:

.. math::

   \mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t\mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}

Here, :math:`\mathbf{\xi}` is the tangential displacement accumulated
during the entire duration of the contact:

.. math::

   \mathbf{\xi} = \int_{t0}^t \mathbf{v}_{t,rel}(\tau) \mathrm{d}\tau

This accumulated tangential displacement must be adjusted to account for
changes in the frame of reference of the contacting pair of particles
during contact.  This occurs due to the overall motion of the contacting
particles in a rigid-body-like fashion during the duration of the
contact.  There are two modes of motion that are relevant: the
'tumbling' rotation of the contacting pair, which changes the
orientation of the plane in which tangential displacement occurs; and
'spinning' rotation of the contacting pair about the vector connecting
their centers of mass (:math:`\mathbf{n}`).  Corrections due to the
former mode of motion are made by rotating the accumulated displacement
into the plane that is tangential to the contact vector at each step, or
equivalently removing any component of the tangential displacement that
lies along :math:`\mathbf{n}`, and rescaling to preserve the magnitude.
This follows the discussion in :ref:`Luding <Luding2008_2>`, see equation
17 and relevant discussion in that work:

.. math::

   \mathbf{\xi} = \left(\mathbf{\xi'} - (\mathbf{n} \cdot \mathbf{\xi'})\mathbf{n}\right) \frac{\|\mathbf{\xi'}\|}{\|\mathbf{\xi'} - (\mathbf{n}\cdot\mathbf{\xi'})\mathbf{n}\|}

Here, :math:`\mathbf{\xi'}` is the accumulated displacement prior to the
current time step and :math:`\mathbf{\xi}` is the corrected
displacement.  Corrections to the displacement due to the second mode of
motion described above (rotations about :math:`\mathbf{n}`) are not
currently implemented, but are expected to be minor for most
simulations.

Furthermore, when the tangential force exceeds the critical force, the
tangential displacement is re-scaled to match the value for the critical
force (see :ref:`Luding <Luding2008_2>`, equation 20 and related
discussion):

.. math::

   \mathbf{\xi} = -\frac{1}{k_t}\left(\mu_t F_{n0}\mathbf{t} - \mathbf{F}_{t,damp}\right)

The tangential force is added to the total normal force (elastic plus
damping) to produce the total force on the particle.

Unlike perfect spheres, the surface normal at the contact point of a
superellipsoid does not generally pass through the particle's center of
mass.  Therefore, both the normal and tangential forces act at the
contact point to induce a torque on each particle.

Using the exact contact point :math:`\mathbf{X}_0` determined by the
geometric solver, the branch vectors from the particle centers of mass
to the contact point are defined as :math:`\mathbf{r}_{ci} =
\mathbf{X}_0 - \mathbf{x}_i` and :math:`\mathbf{r}_{cj} = \mathbf{X}_0 -
\mathbf{x}_j`.  The resulting torques are calculated as:

.. math::

   \mathbf{\tau}_i = \mathbf{r}_{ci} \times \mathbf{F}_{tot}

.. math::

   \mathbf{\tau}_j = -\mathbf{r}_{cj} \times \mathbf{F}_{tot}

----------

If two particles are moving away from each other while in contact, there
is a possibility that the particles could experience an effective
attractive force due to damping.  If the optional *limit_damping*
keyword is used, this option will zero out the normal component of the
force if there is an effective attractive force.

----------

LAMMPS automatically sets pairwise cutoff values for *pair_style
granular/superellipsoid* based on particle radii.  In the vast majority
of situations, this is adequate.  However, a cutoff value can optionally
be appended to the *pair_style granular/superellipsoid* command to
specify a global cutoff (i.e.  a cutoff for all atom types).  This
option may be useful in some rare cases where the automatic cutoff
determination is not sufficient.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for granular pair styles.

Mixing of coefficients is carried out using geometric averaging for most
quantities, e.g. if friction coefficient for type 1-type 1 interactions
is set to :math:`\mu_1`, and friction coefficient for type 2-type 2
interactions is set to :math:`\mu_2`, the friction coefficient for
type1-type2 interactions is computed as :math:`\sqrt{\mu_1\mu_2}`
(unless explicitly specified to a different value by a *pair_coeff 1 2
...* command).

This pair style writes its information to :doc:`binary restart files
<restart>`, so a pair_style command does not need to be specified in an
input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

The ``single()`` function of these pair styles returns 0.0 for the
energy of a pairwise interaction, since energy is not conserved in these
dissipative potentials.  It also returns only the normal component of
the pairwise interaction force.

----------

Restrictions
""""""""""""

The *atom_style* must be set to *ellipsoid superellipsoid* to enable
superellipsoid particles' shape parameters (3 lengths and two blockiness
parameters), see :doc:`atom_style <atom_style>` for more details.

This pair style require Newton's third law be set to *off* for pair
interactions.

There are currently no versions of *fix wall/gran* or *fix
wall/gran/region* that are compatible with superellipsoid particles.

This pair style is part of the ASPHERE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires that atoms store per-particle bounding radius,
shapes, blockiness, inertia, torque, and angular momentum (omega) as
defined by the :doc:`atom_style ellipsoid superellipsoid <atom_style>`.

This pair style requires you to use the :doc:`comm_modify vel yes
<comm_modify>` command so that velocities are stored by ghost atoms.

This pair style will not restart exactly when using the
:doc:`read_restart <read_restart>` command, though it should provide
statistically similar results.  This is because the forces it computes
depend on atom velocities and the atom velocities have been propagated
half a timestep between the force computation and when the restart is
written, due to using Velocity Verlet time integration.  See the
:doc:`read_restart <read_restart>` command for more details.

Accumulated values for individual contacts are saved to restart files
but are not saved to data files.  Therefore, forces may differ
significantly when a system is reloaded using the :doc:`read_data
<read_data>` command.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`
:doc:`pair granular <pair_granular>`

Default
"""""""

For the *pair_coeff* settings: *damping viscoelastic*

References
""""""""""

.. _Brill1996_2:

**(Brilliantov et al, 1996)** Brilliantov, N. V., Spahn, F., Hertzsch,
J. M., & Poschel, T. (1996).  Model for collisions in granular
gases. Physical review E, 53(5), 5382.

.. _Tsuji1992_2:

**(Tsuji et al, 1992)** Tsuji, Y., Tanaka, T., & Ishida,
T. (1992). Lagrangian numerical simulation of plug flow of
cohesionless particles in a horizontal pipe. Powder technology, 71(3),
239-250.

.. _Luding2008_2:

**(Luding, 2008)** Luding, S. (2008). Cohesive, frictional powders:
contact models for tension. Granular matter, 10(4), 235.

.. _Marshall2009_2:

**(Marshall, 2009)** Marshall, J. S. (2009). Discrete-element modeling
of particulate aerosol flows.  Journal of Computational Physics,
228(5), 1541-1561.

.. _Silbert2001_2:

**(Silbert, 2001)** Silbert, L. E., Ertas, D., Grest, G. S., Halsey,
T. C., Levine, D., & Plimpton, S. J. (2001).  Granular flow down an
inclined plane: Bagnold scaling and rheology. Physical Review E,
64(5), 051302.


.. _Thornton1991_2:

**(Thornton, 1991)** Thornton, C. (1991). Interparticle sliding in the
presence of adhesion.  J. Phys. D: Appl. Phys. 24 1942

.. _Thornton2013_2:

**(Thornton et al, 2013)** Thornton, C., Cummins, S. J., & Cleary,
P. W. (2013).  An investigation of the comparative behavior of
alternative contact force models during inelastic collisions. Powder
Technology, 233, 30-46.

.. _WaltonPC_2:

**(Otis R. Walton)** Walton, O.R., Personal Communication

.. _Podlozhnyuk:

**(Podlozhnyuk)** Podlozhnyuk, Pirker, Kloss, Comp. Part. Mech., 4:101-118 (2017).

.. _Houlsby:

**(Houlsby)** Houlsby, Computers and Geotechnics, 36, 953-959 (2009).

.. _GeometricTools:

**(Eberly)** Eberly, Geometric Tools: Dynamic Collision Detection Using Oriented Bounding Boxes (2008).

