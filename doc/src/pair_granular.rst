.. index:: pair_style granular

pair_style granular command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style granular cutoff

* cutoff = global cutoff (optional).  See discussion below.

Examples
""""""""

.. code-block:: LAMMPS

   pair_style granular
   pair_coeff * * hooke 1000.0 50.0 tangential linear_nohistory 1.0 0.4 damping mass_velocity

   pair_style granular
   pair_coeff * * hooke 1000.0 50.0 tangential linear_history 500.0 1.0 0.4 damping mass_velocity

   pair_style granular
   pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 limit_damping

   pair_style granular
   pair_coeff * * hertz/material 1e8 0.3 0.3 tangential mindlin_rescale NULL 1.0 0.4 damping tsuji

   pair_style granular
   pair_coeff * * hertz/material 1e8 0.3 0.3 tangential mindlin_rescale NULL 1.0 0.4 damping coeff_restitution synchronized_verlet

   pair_style granular
   pair_coeff 1 * jkr 1000.0 500.0 0.3 10 tangential mindlin 800.0 1.0 0.5 rolling sds 500.0 200.0 0.5 twisting marshall
   pair_coeff 2 2 hertz 200.0 100.0 tangential linear_history 300.0 1.0 0.1 rolling sds 200.0 100.0 0.1 twisting marshall

   pair_style granular
   pair_coeff 1 1 dmt 1000.0 50.0 0.3 0.0 tangential mindlin NULL 0.5 0.5 rolling sds 500.0 200.0 0.5 twisting marshall
   pair_coeff 2 2 dmt 1000.0 50.0 0.3 10.0 tangential mindlin NULL 0.5 0.1 rolling sds 500.0 200.0 0.1 twisting marshall

   pair_style granular
   pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 heat area 0.1

   pair_style granular
   pair_coeff * * mdr 5e6 0.4 1.9e5 2.0 0.5 0.5 tangential linear_history 940.0 1.0 0.7 rolling sds 2.7e5 0.0 0.6 damping mdr 1

Description
"""""""""""

The *granular* styles support a variety of options for the normal,
tangential, rolling and twisting forces resulting from contact between
two granular particles. This expands on the options offered by the
:doc:`pair gran/\* <pair_gran>` pair styles. The total computed forces
and torques are the sum of various models selected for the
:doc:`normal <granular_normal_models>`, :doc:`damping <granular_damping_models>`,
:doc:`tangential <granular_tangential_models>`,
:doc:`rolling <granular_rolling_models>`, and
:doc:`twisting <granular_twisting_models>` interactions. Optionally,
:doc:`heat transport <granular_heat_models>` between particles can also
be modeled.

.. toctree::
   :maxdepth: 1
   :hidden:

   granular_normal_models
   granular_damping_models
   granular_tangential_models
   granular_rolling_models
   granular_twisting_models
   granular_heat_models

All model choices and parameters are entered in the :doc:`pair_coeff
<pair_coeff>` command, as described below.  Unlike e.g. :doc:`pair
gran/hooke <pair_gran>`, coefficient values are not global, but can be
set to different values for different combinations of particle types, as
determined by the :doc:`pair_coeff <pair_coeff>` command.  If the
contact model choice is the same for two particle types, the mixing for
the cross-coefficients can be carried out automatically.  This is shown
in one of the examples, where model choices are the same for type 1 -
type 1 as for type 2 - type2 interactions, but coefficients are
different.  In this case, the mixed coefficients for type 1 - type 2
interactions can be determined from mixing rules discussed below.  For
additional flexibility, coefficients as well as model forms can vary
between particle types, as shown in the sixth example: type 1 - type 1
interactions are based on a Johnson-Kendall-Roberts normal contact model
and 2-2 interactions are based on a DMT cohesive model (see below).  In
that example, 1-1 and 2-2 interactions have different model forms, in
which case mixing of coefficients cannot be determined, so 1-2
interactions must be explicitly defined via the *pair_coeff 1 \**
command, otherwise an error would result.

----------

:doc:`Normal contact models <granular_normal_models>`
-----------------------------------------------------

The first required keyword for the *pair_coeff* command is the normal
contact model, i.e. the force-displacement relation associated with
motion of particles along the normal direction of the contact plane.
Currently supported options for normal contact models
and their required arguments are:

1. :ref:`hooke <hooke_normal_model>` : :math:`k_n`, :math:`\eta_{n0}` (or :math:`e`)
2. :ref:`hertz <hertz_normal_model>` : :math:`k_n`, :math:`\eta_{n0}` (or :math:`e`)
3. :ref:`hertz/material <hertz_material_normal_model>` : :math:`E`, :math:`\eta_{n0}` (or :math:`e`), :math:`\nu`
4. :ref:`dmt <dmt_normal_model>` : :math:`E`, :math:`\eta_{n0}` (or :math:`e`), :math:`\nu`, :math:`\gamma`
5. :ref:`jkr <jkr_normal_model>` : :math:`E`, :math:`\eta_{n0}` (or :math:`e`), :math:`\nu`, :math:`\gamma`
6. :ref:`mdr <mdr_normal_model>` : :math:`E`, :math:`\nu`, :math:`Y`, :math:`\Delta\gamma`,
   :math:`\psi_b`, :math:`\eta_{n0}`
7. :ref:`epa_linear <epa_linear_normal_model>` : :math:`k_1`, :math:`\eta_{n0}` (or :math:`e`), :math:`\hat{k_2}`,
   :math:`k_c`, :math:`\phi_f`, :math:`f_0`
8. :ref:`epa_nonlinear <epa_nonlinear_normal_model>` : :math:`E`, :math:`\eta_{n0}` (or :math:`e`), :math:`\nu`,
   :math:`\lambda_p`, :math:`f_0`, :math:`k_{adh}`, :math:`m`, :math:`n`

.. note::

   Typically, neighbor lists are constructed for pair granular by testing
   whether finite sized particles overlap (using their radii). However,
   this is not the case for normal models which can interact beyond
   contact, e.g. *jkr*. Instead, the maximum radius for each particle
   type is first calculated then used to calculate a maximum per-type
   cutoff distance. For polydisperse systems, this affects the performance
   of the :doc:`multi neighbor <neigh_modify>` option where one should
   assign atoms of similar radii the same type. See the
   :doc:`pair lj/cut/sphere <pair_lj_cut_sphere>` page for a related discussion.

----------

:doc:`Damping models <granular_damping_models>`
-----------------------------------------------

In addition, the normal force is augmented by a damping term of the
following general form:

.. math::

   \mathbf{F}_{n,damp} = -\eta_n \mathbf{v}_{n,rel}

Here, :math:`\mathbf{v}_{n,rel} = (\mathbf{v}_j - \mathbf{v}_i) \cdot
\mathbf{n}\ \mathbf{n}` is the component of relative velocity along
:math:`\mathbf{n}`.

The optional *damping* keyword to the *pair_coeff* command followed by
a keyword determines the model form of the damping factor :math:`\eta_n`,
and the interpretation of the :math:`\eta_{n0}` or :math:`e` coefficients
specified as part of the normal contact model settings. The *damping*
keyword and corresponding model form selection may be appended
anywhere in the *pair coeff* command.  Note that the choice of damping
model affects both the normal and tangential damping (and depending on
other settings, potentially also the twisting damping).  The options
for the damping model currently supported are as follows:

1. :ref:`velocity <velocity_damping_model>`
2. :ref:`mass_velocity <mass_velocity_damping_model>`
3. :ref:`viscoelastic <viscoelastic_damping_model>`
4. :ref:`tsuji <tsuji_damping_model>`
5. :ref:`coeff_restitution <coeff_restitution_damping_model>`
6. :ref:`mdr (class) <mdr_damping_model>`: :math:`d_{type}`

If the *damping* keyword is not specified, the *viscoelastic* model is
used by default.

----------

:doc:`Tangential models <granular_tangential_models>`
-----------------------------------------------------

The *pair_coeff* command requires specification of the tangential
contact model. The required keyword *tangential* is required, followed
by the model choice and associated parameters. Currently supported
tangential model choices and their expected parameters are as follows:

1. :ref:`linear_nohistory <linear_nohistory_tangential_model>` : :math:`x_{\gamma,t}`, :math:`\mu_t`
2. :ref:`linear_history <linear_history_tangential_model>` : :math:`k_t`, :math:`x_{\gamma,t}`, :math:`\mu_t`
3. :ref:`mindlin <mindlin_tangential_model>` : :math:`k_t` or NULL, :math:`x_{\gamma,t}`, :math:`\mu_t`
4. :ref:`mindlin/force <mindlin_force_tangential_model>` : :math:`k_t` or NULL, :math:`x_{\gamma,t}`, :math:`\mu_t`
5. :ref:`mindlin_rescale <mindlin_rescale_tangential_model>` : :math:`k_t` or NULL, :math:`x_{\gamma,t}`, :math:`\mu_t`
6. :ref:`mindlin_rescale/force <mindlin_rescale_force_tangential_model>` : :math:`k_t` or NULL, :math:`x_{\gamma,t}`, :math:`\mu_t`
7. :ref:`linear_history/static <linear_history_static_tangential_model>` : :math:`k_t`, :math:`x_{\gamma,t}`, :math:`\mu_{t,static}`, :math:`\mu_{t,dynamic}`
8. :ref:`mindlin/static <mindlin_static_tangential_model>` : :math:`k_t` or NULL, :math:`x_{\gamma,t}`, :math:`\mu_{t,static}`, :math:`\mu_{t,dynamic}`

----------

:doc:`Rolling friction models <granular_rolling_models>`
--------------------------------------------------------

The optional *rolling* keyword enables rolling friction, which resists
pure rolling motion of particles.  The options currently supported are:

1. *none*
2. :ref:`sds <sds_rolling_model>` : :math:`k_{roll}`, :math:`\gamma_{roll}`, :math:`\mu_{roll}`

If the *rolling* keyword is not specified, the model defaults to *none*\
.

:doc:`Twisting friction models <granular_twisting_models>`
----------------------------------------------------------

The optional *twisting* keyword enables twisting friction, which resists
rotation of two contacting particles about the vector :math:`\mathbf{n}`
that connects their centers.  The options currently supported are:

1. *none*
2. :ref:`sds <sds_twisting_model>` : :math:`k_{twist}`, :math:`\gamma_{twist}`, :math:`\mu_{twist}`
3. :ref:`marshall <marshall_twisting_model>`

If the *twisting* keyword is not specified, the model defaults to
*none*\ .

:doc:`Heat conduction models <granular_heat_models>`
----------------------------------------------------

The optional *heat* keyword enables heat conduction.  The options
currently supported are:

1. *none*
2. :ref:`radius <radius_heat_model>` : :math:`k_{s}`
3. :ref:`area <area_heat_model>` : :math:`h_{s}`

If the *heat* keyword is not specified, the model defaults to
*none*.  All heat models calculate an additional pairwise quantity
accessible by the single() function (described below) which is the heat
conducted between the two particles.

Note that the option *none* must either be used in all or none of the
*pair_coeff* calls.  See :doc:`fix heat/flow <fix_heat_flow>` and
:doc:`fix property/atom <fix_property_atom>` for more information on
this option.


----------

If two particles are moving away from each other while in contact, there
is a possibility that the particles could experience an effective
attractive force due to damping.  If the optional *limit_damping*
keyword is used, this option will zero out the normal component of the
force if there is an effective attractive force.  This keyword cannot be
used with the JKR or DMT models.

----------

The standard velocity-Verlet integration scheme's half-step staggering
of position and velocity can introduce inaccuracies in frictional
tangential force calculations, resulting in unphysical kinematics in
certain systems.  These effects are particularly pronounced in
polydisperse frictional flows characterized by large-to-small size
ratios exceeding three.  The *synchronized_verlet* flag implements an
alternate Velocity-Verlet integration scheme, as detailed in :ref:`Vyas
et al <Vyas2025>`, that synchronizes position and velocity updates for
force evaluation.  By refining tangential force calculations, the
*synchronized_verlet* method ensures physically consistent results
without significantly impacting computational cost.

----------

The *granular* pair style can reproduce the behavior of the *pair
gran/\** styles with the appropriate settings (some very minor
differences can be expected due to corrections in displacement history
frame-of-reference, and the application of the torque at the center of
the contact rather than at each particle).  The first example above is
equivalent to *pair gran/hooke 1000.0 NULL 50.0 50.0 0.4 1*\ .  The
second example is equivalent to *pair gran/hooke/history 1000.0 500.0
50.0 50.0 0.4 1*\ .  The third example is equivalent to *pair
gran/hertz/history 1000.0 500.0 50.0 50.0 0.4 1 limit_damping*\ .

----------

LAMMPS automatically sets pairwise cutoff values for *pair_style
granular* based on particle radii (and in the case of *jkr* pull-off
distances).  In the vast majority of situations, this is adequate.
However, a cutoff value can optionally be appended to the *pair_style
granular* command to specify a global cutoff (i.e. a cutoff for all atom
types).  Additionally, the optional *cutoff* keyword can be passed to
the *pair_coeff* command, followed by a cutoff value.  This will set a
pairwise cutoff for the atom types in the *pair_coeff* command.  These
options may be useful in some rare cases where the automatic cutoff
determination is not sufficient, e.g.  if particle diameters are being
modified via the *fix adapt* command.  In that case, the global cutoff
specified as part of the *pair_style granular* command is applied to all
atom types, unless it is overridden for a given atom type combination by
the *cutoff* value specified in the *pair coeff* command.  If *cutoff*
is only specified in the *pair coeff* command and no global cutoff is
appended to the *pair_style granular* command, then LAMMPS will use that
cutoff for the specified atom type combination, and automatically set
pairwise cutoffs for the remaining atom types.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for any granular pair styles.

Mixing of coefficients is carried out using geometric averaging for most
quantities, e.g. if friction coefficient for type 1-type 1 interactions
is set to :math:`\mu_1`, and friction coefficient for type 2-type 2
interactions is set to :math:`\mu_2`, the friction coefficient for
type1-type2 interactions is computed as :math:`\sqrt{\mu_1\mu_2}`
(unless explicitly specified to a different value by a *pair_coeff 1 2
...* command).  The exception to this is elastic modulus, only
applicable to *hertz/material*, *dmt*, *jkr* and *epa_linear*
normal contact models. In that case, the effective elastic modulus is
computed as:

.. math::

   E_\text{eff,ij} = \left(\frac{1-\nu_i^2}{E_i} + \frac{1-\nu_j^2}{E_j}\right)^{-1}

If the *i-j* coefficients :math:`E_{ij}` and :math:`\nu_{ij}` are
explicitly specified, the effective modulus is computed as:

.. math::

   E_\text{eff,ij} = \left(\frac{1-\nu_{ij}^2}{E_{ij}} + \frac{1-\nu_{ij}^2}{E_{ij}}\right)^{-1}

or

.. math::

   E_\text{eff,ij} = \frac{E_{ij}}{2(1-\nu_{ij}^2)}

This pair style writes its information to :doc:`binary restart files
<restart>`, so a pair_style command does not need to be specified in an
input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

The single() function of this pair style returns 0.0 for the energy of a
pairwise interaction, since energy is not conserved in these dissipative
potentials.  It also returns only the normal component of the pairwise
interaction force.  However, the single() function also calculates at
least 12 extra pairwise quantities.  The first 3 are the components of
the tangential force between particles I and J, acting on particle I.
The fourth is the magnitude of this tangential force.  The next 3 (5-7)
are the components of the rolling torque acting on particle I.  The next
entry (8) is the magnitude of the rolling torque.  The next entry (9) is
the magnitude of the twisting torque acting about the vector connecting
the two particle centers.  The next 3 (10-12) are the components of the
vector connecting the centers of the two particles (x_I - x_J).  If a
granular sub-model calculates additional contact information (e.g. the
heat sub-models calculate the amount of heat exchanged), these
quantities are appended to the end of this list.  First, any extra
values from the normal sub-model are appended followed by the damping,
tangential, rolling, twisting, then heat models.  See the descriptions
of individual granular sub-models linked above for information on any extra
quantities. If two or more models are defined by pair coefficients, the
size of the array is set by the maximum number of extra quantities in a
model but the order of quantities is determined by each model's specific
set of sub-models.  Any unused quantities are zeroed.

These extra quantities can be accessed by the :doc:`compute pair/local
<compute_pair_local>` command, as *p1*, *p2*, ..., *p12*\ .

----------

Restrictions
""""""""""""

This pair style is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style requires that atoms store per-particle radius, torque,
and angular velocity (omega) as defined by the :doc:`atom_style sphere
<atom_style>`.

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
:doc:`pair gran/\* <pair_gran>`

Default
"""""""

For the *pair_coeff* settings: *damping viscoelastic*, *rolling none*,
*twisting none*\ .

References
""""""""""

.. _Marshall2009_1:

**(Marshall, 2009)** Marshall, J. S. (2009). Discrete-element modeling
of particulate aerosol flows.  Journal of Computational Physics,
228(5), 1541-1561.

.. _Silbert2001_1:

**(Silbert, 2001)** Silbert, L. E., Ertas, D., Grest, G. S., Halsey,
T. C., Levine, D., & Plimpton, S. J. (2001).  Granular flow down an
inclined plane: Bagnold scaling and rheology. Physical Review E,
64(5), 051302.

.. _Thornton1991:

**(Thornton, 1991)** Thornton, C. (1991). Interparticle sliding in the
presence of adhesion.  J. Phys. D: Appl. Phys. 24 1942

.. _VargasMcCarthy2001:

**(Vargas and McCarthy 2001)** Vargas, W.L. and McCarthy, J.J. (2001).
Heat conduction in granular materials.
AIChE Journal, 47(5), 1052-1059.

.. _Vyas2025:

**(Vyas et al, 2025)**  Vyas D. R., Ottino J. M., Lueptow R. M., & Umbanhowar P. B. (2025).
Improved Velocity-Verlet Algorithm for the Discrete Element Method.
Computer Physics Communications, 109524.
