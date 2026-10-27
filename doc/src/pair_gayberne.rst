.. index:: pair_style gayberne
.. index:: pair_style gayberne/gpu
.. index:: pair_style gayberne/intel
.. index:: pair_style gayberne/omp

pair_style gayberne command
===========================

Accelerator Variants: *gayberne/gpu*, *gayberne/intel*, *gayberne/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style gayberne gamma upsilon mu cutoff

* gamma = shift for potential minimum (typically 1)
* upsilon = exponent for eta orientation-dependent energy function
* mu = exponent for chi orientation-dependent energy function
* cutoff = global cutoff for interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style gayberne 1.0 1.0 1.0 10.0
   pair_coeff * * 1.0 1.7 1.7 3.4 3.4 1.0 1.0 1.0

Description
"""""""""""

The *gayberne* styles compute a Gay-Berne anisotropic LJ interaction
:ref:`(Berardi) <Berardi>` between pairs of ellipsoidal particles or an
ellipsoidal and spherical particle via the formulas

.. math::

   U ( \mathbf{A}_1, \mathbf{A}_2, \mathbf{r}_{12} ) = & U_r (
   \mathbf{A}_1, \mathbf{A}_2, \mathbf{r}_{12}, \gamma ) \cdot \eta_{12} (
   \mathbf{A}_1, \mathbf{A}_2, \upsilon ) \cdot \chi_{12} ( \mathbf{A}_1,
   \mathbf{A}_2, \mathbf{r}_{12}, \mu ) \\
   U_r = & 4 \epsilon ( \varrho^{12} - \varrho^6) \\
   \varrho = & \frac{\sigma}{ h_{12} + \gamma \sigma}

where :math:`\mathbf{A}_1` and :math:`\mathbf{A}_2` are the
transformation matrices from the simulation box frame to the body frame
and :math:`r_{12}` is the center to center vector between the particles.
:math:`U_r` controls the shifted distance dependent interaction based on
the distance of closest approach of the two particles (:math:`h_{12}`)
and the user-specified shift parameter :math:`\gamma`.  When both
particles are spherical, the formula reduces to the usual Lennard-Jones
interaction (see details below for when Gay-Berne treats a particle as
"spherical").

For large uniform molecules it has been shown that the energy
parameters are approximately representable in terms of local contact
curvatures :ref:`(Everaers) <Everaers2>`:

.. math::

   \epsilon_a = \sigma \cdot { \frac{a}{ b \cdot c } }; \epsilon_b =
   \sigma \cdot { \frac{b}{ a \cdot c } }; \epsilon_c = \sigma \cdot {
   \frac{c}{ a \cdot b } }

The variable names utilized as potential parameters are for the most
part taken from :ref:`(Everaers) <Everaers2>` in order to be consistent with
the :doc:`RE-squared pair potential <pair_resquared>`.  The role of
the :math:`\upsilon` and :math:`\mu` parameters and more details of the
Gay-Berne formulation as implemented in LAMMPS are given in the
:ref:`Implementation details <gayberne_details>` section below and in
the references listed at the end of this page.

Use of this pair style requires the NVE, NVT, or NPT fixes with the
*asphere* extension (e.g. :doc:`fix nve/asphere <fix_nve_asphere>`) in
order to integrate particle rotation.  Additionally, :doc:`atom_style
ellipsoid <atom_style>` should be used since it defines the rotational
state and the size and shape of each ellipsoidal particle.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below:

* :math:`\epsilon` = well depth (energy units)
* :math:`\sigma` = minimum effective particle radii (distance units)
* :math:`\epsilon_{i,a}` = relative well depth of type I for side-to-side interactions
* :math:`\epsilon_{i,b}` = relative well depth of type I for face-to-face interactions
* :math:`\epsilon_{i,c}` = relative well depth of type I for end-to-end interactions
* :math:`\epsilon_{j,a}` = relative well depth of type J for side-to-side interactions
* :math:`\epsilon_{j,b}` = relative well depth of type J for face-to-face interactions
* :math:`\epsilon_{j,c}` = relative well depth of type J for end-to-end interactions
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global
cutoff specified in the pair_style command is used.

It is typical with the Gay-Berne potential to define :math:`\sigma` as
the minimum of the 3 shape diameters of the particles involved in an I,I
interaction, though this is not required.  Note that this is a different
meaning for :math:`\sigma` than the :doc:`pair_style resquared
<pair_resquared>` potential uses.

The :math:`\epsilon_i` and :math:`\epsilon_j` coefficients are actually
defined for atom types, not for pairs of atom types.  Thus, in a series
of pair_coeff commands, they only need to be specified once for each
atom type.

Specifically, if any of :math:`\epsilon_{i,a}`, :math:`\epsilon_{i,b}`,
:math:`\epsilon_{i,c}` are non-zero, the three values are assigned to
atom type I.  If all the :math:`\epsilon_i` values are zero, they are
ignored.  If any of :math:`\epsilon_{j,a}`, :math:`\epsilon_{j,b}`,
:math:`\epsilon_{j,c}` are non-zero, the three values are assigned to
atom type J.  If all three epsilon_j values are zero, they are ignored.
Thus the typical way to define the :math:`\epsilon_i` and
:math:`\epsilon_j` coefficients is to list their values in "pair_coeff
I J" commands when I = J, but set them to 0.0 when I != J.  If you do
list them when I != J, you should ensure they are consistent with their
values in other pair_coeff commands, since only the last setting will
be in effect.

Note that if this potential is being used as a sub-style of
:doc:`pair_style hybrid <pair_hybrid>`, and there is no "pair_coeff I I"
setting made for Gay-Berne for a particular type I (because I-I
interactions are computed by another hybrid pair potential), then you
still need to ensure the :math:`\epsilon` a,b,c coefficients are assigned to
that type. e.g. in a "pair_coeff I J" command.

.. note::

   If the :math:`\epsilon_{a}` = :math:`\epsilon_{b}` = :math:`\epsilon_{c}` for an atom type, and if the shape
   of the particle itself is spherical, meaning its 3 shape parameters
   are all the same, then the particle is treated as an LJ sphere by the
   Gay-Berne potential.  This is significant because if two LJ spheres
   interact, then the simple Lennard-Jones formula is used to compute
   their interaction energy/force using the specified epsilon and sigma
   as the standard LJ parameters.  This is much cheaper to compute than
   the full Gay-Berne formula.  To treat the particle as a LJ sphere
   with sigma = D, you should normally set :math:`\epsilon_{a}` = :math:`\epsilon_{b}` = :math:`\epsilon_{c}` =
   1.0, set the pair_coeff :math:`\sigma = D`, and also set the 3 shape
   parameters for the particle to D.  The one exception is that if the 3
   shape parameters are set to 0.0, which is a valid way in LAMMPS to
   specify a point particle, then the Gay-Berne potential will treat
   that particle as a sphere with radius 1.0, i.e. the same as shape
   parameters of 2.0, since it requires finite-size particles.  This
   size only matters for interactions with ellipsoids, since the
   interaction between two such particles is computed with the
   Lennard-Jones formula and the pair_coeff :math:`\sigma`.  To be
   consistent with the rule above, you should set the pair_coeff
   :math:`\sigma` to 2.0 in this case.

----------

.. _gayberne_details:

Implementation details
""""""""""""""""""""""

The following is a summary of the Gay-Berne formulas as implemented in
LAMMPS by Mike Brown (Sandia National Labs) in 2007.  Let the shape
matrices :math:`\mathbf{S}_i = \mathrm{diag}(a_i, b_i, c_i)` be given by
the radii of particle *i*, i.e. half of the diameters set with the
:doc:`set shape <set>` command or in the data file.  Let the relative
energy matrices :math:`\mathbf{E}_i = \mathrm{diag}(\epsilon_{i,a}^{-1/\mu},
\epsilon_{i,b}^{-1/\mu}, \epsilon_{i,c}^{-1/\mu})` be given by the
relative well depths of the :doc:`pair_coeff <pair_coeff>` command.
With the rotation matrices :math:`\mathbf{A}_i` from the simulation box
frame to the body frame of particle *i*, the following matrices are
defined:

.. math::

   \mathbf{G}_i = \mathbf{A}_i^T \mathbf{S}_i^2 \mathbf{A}_i, \qquad
   \mathbf{B}_i = \mathbf{A}_i^T \mathbf{E}_i \mathbf{A}_i, \qquad
   \mathbf{G}_{12} = \mathbf{G}_1 + \mathbf{G}_2, \qquad
   \mathbf{B}_{12} = \mathbf{B}_1 + \mathbf{B}_2.

The distance of closest approach :math:`h_{12}` in the formula for
:math:`U_r` above is computed with the approximation of :ref:`(Perram)
<Perram>`:

.. math::

   h_{12} = r - \sigma_{12}, \qquad \sigma_{12} = \left[ \frac{1}{2}
   \hat{\mathbf{r}}_{12}^T \mathbf{G}_{12}^{-1} \hat{\mathbf{r}}_{12}
   \right]^{-1/2},

where :math:`\mathbf{r}_{12}` is the vector from the center of particle
1 to the center of particle 2, :math:`r = |\mathbf{r}_{12}|`, and
:math:`\hat{\mathbf{r}}_{12} = \mathbf{r}_{12}/r`.  The
orientation-dependent energy terms with the exponents :math:`\upsilon`
and :math:`\mu` from the pair_style command are

.. math::

   \eta_{12} = & \left[ \frac{2 s_1 s_2}{\det(\mathbf{G}_{12})}
   \right]^{\upsilon/2}, \qquad s_i = \left( a_i b_i + c_i^2 \right)
   \sqrt{a_i b_i}, \\
   \chi_{12} = & \left[ 2 \hat{\mathbf{r}}_{12}^T \mathbf{B}_{12}^{-1}
   \hat{\mathbf{r}}_{12} \right]^{\mu}.

For a spherical particle with equal relative well depths,
:math:`\mathbf{G}_i` and :math:`\mathbf{B}_i` are multiples of the unit
matrix, which do not depend on the orientation of the particle.

Since the analytic forces and torques were not available in the
literature when this pair style was implemented, they are listed here.
With the abbreviations

.. math::

   \boldsymbol{\kappa} = \mathbf{G}_{12}^{-1} \mathbf{r}_{12}, \qquad
   \boldsymbol{\iota} = \mathbf{B}_{12}^{-1} \mathbf{r}_{12}, \qquad
   u' = \frac{24 \epsilon}{\sigma} \left( 2 \varrho^{13} - \varrho^7
   \right) = -\frac{\partial U_r}{\partial h_{12}},

the force on particle 1 is

.. math::

   \mathbf{f}_1 = - \eta_{12} \left( U_r \frac{\partial
   \chi_{12}}{\partial \mathbf{x}_1} + \chi_{12} \frac{\partial
   U_r}{\partial \mathbf{x}_1} \right)

with the derivatives with respect to the position :math:`\mathbf{x}_1`
of particle 1 (see also :ref:`(Allen) <Allen3>`):

.. math::

   \frac{\partial U_r}{\partial \mathbf{x}_1} = & u' \,
   \hat{\mathbf{r}}_{12} + \frac{u' \sigma_{12}^3}{2 r^2} \left[
   \boldsymbol{\kappa} - \left( \boldsymbol{\kappa} \cdot
   \hat{\mathbf{r}}_{12} \right) \hat{\mathbf{r}}_{12} \right], \\
   \frac{\partial \chi_{12}}{\partial \mathbf{x}_1} = & - \frac{4
   \mu}{r^2} \chi_{12}^{(\mu - 1)/\mu} \left[ \boldsymbol{\iota} -
   \left( \boldsymbol{\iota} \cdot \hat{\mathbf{r}}_{12} \right)
   \hat{\mathbf{r}}_{12} \right].

The force on particle 2 is :math:`\mathbf{f}_2 = -\mathbf{f}_1`.  The
torque on particle *i* is

.. math::

   \boldsymbol{\tau}_i = - \left( U_r \eta_{12} \frac{\partial
   \chi_{12}}{\partial \boldsymbol{\phi}_i} + U_r \chi_{12} \frac{\partial
   \eta_{12}}{\partial \boldsymbol{\phi}_i} + \eta_{12} \chi_{12}
   \frac{\partial U_r}{\partial \boldsymbol{\phi}_i} \right),

where :math:`\partial / \partial \boldsymbol{\phi}_i` denotes the
derivative with respect to an infinitesimal rotation of particle *i*.
Let :math:`\mathbf{a}_{m,i}` be the *m*-th row of :math:`\mathbf{A}_i`,
i.e. the *m*-th axis of particle *i* in the simulation box frame, and
:math:`s_{m,i}` the corresponding radius.  Then

.. math::

   \frac{\partial U_r}{\partial \boldsymbol{\phi}_i} = & - \frac{u'
   \sigma_{12}^3}{2 r^2} \, \boldsymbol{\kappa} \times \left(
   \mathbf{G}_i \boldsymbol{\kappa} \right), \\
   \frac{\partial \chi_{12}}{\partial \boldsymbol{\phi}_i} = & \frac{4
   \mu}{r^2} \chi_{12}^{(\mu - 1)/\mu} \, \boldsymbol{\iota} \times
   \left( \mathbf{B}_i \boldsymbol{\iota} \right), \\
   \frac{\partial \eta_{12}}{\partial \boldsymbol{\phi}_i} = &
   \sum_{m=1}^3 \mathbf{a}_{m,i} \times \frac{\partial \eta_{12}}{\partial
   \mathbf{a}_{m,i}} = - \upsilon \, \eta_{12} \sum_{m=1}^3 s_{m,i}^2 \,
   \mathbf{a}_{m,i} \times \left( \mathbf{G}_{12}^{-1} \mathbf{a}_{m,i}
   \right).

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for this pair style can be mixed.  The default mix
value is *geometric*\ .  See the "pair_modify" command for details.

This pair style supports the :doc:`pair_modify <pair_modify>` shift
option for the energy of the Lennard-Jones portion of the pair
interaction, but only for sphere-sphere interactions.  There is no
shifting performed for ellipsoidal interactions due to the anisotropic
dependence of the interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *gayberne* style is part of the ASPHERE package.  It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

These pair styles require that atoms store torque and a quaternion to
represent their orientation, as defined by the
:doc:`atom_style <atom_style>`.  It also require they store a per-type
:doc:`shape <set>`.  The particles cannot store a per-particle
diameter.

This pair style requires that atoms be ellipsoids as defined by the
:doc:`atom_style ellipsoid <atom_style>` command.

Particles acted on by the potential can be finite-size aspherical or
spherical particles, or point particles.  Spherical particles have all
3 of their shape parameters equal to each other.  Point particles have
all 3 of their shape parameters equal to 0.0.

The Gay-Berne potential does not become isotropic as r increases
:ref:`(Everaers) <Everaers2>`.  The distance-of-closest-approach
approximation used by LAMMPS becomes less accurate when high-aspect
ratio ellipsoids are used.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix nve/asphere <fix_nve_asphere>`,
:doc:`compute temp/asphere <compute_temp_asphere>`, :doc:`pair_style resquared <pair_resquared>`

Default
"""""""

none

----------

.. _Everaers2:

**(Everaers)** Everaers and Ejtehadi, Phys Rev E, 67, 041710 (2003).

.. _Berardi:

**(Berardi)** Berardi, Fava, Zannoni, Chem Phys Lett, 297, 8-14 (1998).
Berardi, Muccioli, Zannoni, J Chem Phys, 128, 024905 (2008).

.. _Perram:

**(Perram)** Perram and Rasmussen, Phys Rev E, 54, 6565-6572 (1996).

.. _Allen3:

**(Allen)** Allen and Germano, Mol Phys 104, 3225-3235 (2006).
