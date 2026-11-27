.. index:: pair_style resquared
.. index:: pair_style resquared/gpu
.. index:: pair_style resquared/omp

pair_style resquared command
============================

Accelerator Variants: *resquared/gpu*, *resquared/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style resquared cutoff

* cutoff = global cutoff for interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style resquared 10.0
   pair_coeff * * 1.0 1.0 1.7 3.4 3.4 1.0 1.0 1.0

Description
"""""""""""

Style *resquared* computes the RE-squared anisotropic interaction
:ref:`(Everaers) <Everaers3>`, :ref:`(Babadi) <Babadi>` between pairs of
ellipsoidal and/or spherical Lennard-Jones particles.  For ellipsoidal
interactions, the potential considers the ellipsoid as being comprised
of small spheres of size :math:`\sigma`.  LJ particles are a single sphere of
size :math:`\sigma`.  The distinction is made to allow the pair style to make
efficient calculations of ellipsoid/solvent interactions.

Details for the equations used are given in the references below and
in the :ref:`Implementation details <resquared_details>` section below.

Use of this pair style requires the NVE, NVT, or NPT fixes with the
*asphere* extension (e.g. :doc:`fix nve/asphere <fix_nve_asphere>`) in
order to integrate particle rotation.  Additionally, :doc:`atom_style ellipsoid <atom_style>` should be used since it defines the
rotational state and the size and shape of each ellipsoidal particle.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* A12 = Energy Prefactor/Hamaker constant (energy units)
* :math:`\sigma` = atomic interaction diameter (distance units)
* :math:`\epsilon_{i,a}` = relative well depth of type I for side-to-side interactions
* :math:`\epsilon_{i,b}` = relative well depth of type I for face-to-face interactions
* :math:`\epsilon_{i,c}` = relative well depth of type I for end-to-end interactions
* :math:`\epsilon_{j,a}` = relative well depth of type J for side-to-side interactions
* :math:`\epsilon_{j,b}` = relative well depth of type J for face-to-face interactions
* :math:`\epsilon_{j,c}` = relative well depth of type J for end-to-end interactions
* cutoff (distance units)

The parameters used depend on the type of the interacting particles,
i.e. ellipsoids or LJ spheres.  The type of a particle is determined
by the diameters specified for its 3 shape parameters.  If all 3 shape
parameters = 0.0, then the particle is treated as an LJ sphere.  The
:math:`\epsilon_{i,*}` or :math:`\epsilon_{j,*}` parameters are ignored
for LJ spheres.  If
the 3 shape parameters are > 0.0, then the particle is treated as an
ellipsoid (even if the 3 parameters are equal to each other).

A12 specifies the energy prefactor which depends on the types of the
two interacting particles.

For ellipsoid/ellipsoid interactions, the interaction is computed by
the formulas in the :ref:`Implementation details <resquared_details>`
section below.  A12 is
the Hamaker constant as described in :ref:`(Everaers) <Everaers3>`. In LJ
units:

.. math::

   A_{12} = 4\pi^2\epsilon_{\mathrm{LJ}}(\rho\sigma^3)^2

where :math:`\rho` gives the number density of the spherical particles
composing the ellipsoids and :math:`\epsilon_{\mathrm{LJ}}` determines
the interaction strength of the spherical particles.

For ellipsoid/LJ sphere interactions, the interaction is also computed
by the formulas in the :ref:`Implementation details <resquared_details>`
section below.  A12 has a modified form:

.. math::

   A_{12} = 4\pi^2\epsilon_{\mathrm{LJ}}(\rho\sigma^3)

For ellipsoid/LJ sphere interactions, a correction to the distance-
of-closest approach equation has been implemented to reduce the error
from two particles of disparate sizes; see the :ref:`Implementation
details <resquared_details>` section below.

For LJ sphere/LJ sphere interactions, the interaction is computed
using the standard Lennard-Jones formula, which is much cheaper to
compute than the ellipsoidal formulas.  A12 is used as epsilon in the
standard LJ formula:

.. math::

   A_{12} = \epsilon_{\mathrm{LJ}}

and the specified :math:`\sigma` is used as the :math:`\sigma` in the
standard LJ formula.

When one of both of the interacting particles are ellipsoids, then
:math:`\sigma` specifies the diameter of the continuous distribution of
constituent particles within each ellipsoid used to model the RE-squared
potential.  Note that this is a different meaning for :math:`\sigma`
than the :doc:`pair_style gayberne <pair_gayberne>` potential uses.

The :math:`\epsilon_i` and :math:`\epsilon_j` coefficients are defined
for atom types, not for pairs of atom types.  Thus, in a series of
pair_coeff commands, they only need to be specified once for each atom
type.

Specifically, if any of :math:`\epsilon_{i,a}`, :math:`\epsilon_{i,b}`,
:math:`\epsilon_{i,c}` are non-zero, the three values are assigned to
atom type I.  If all the :math:`\epsilon_i` values are zero, they are
ignored.  If any of :math:`\epsilon_{j,a}`, :math:`\epsilon_{j,b}`,
:math:`\epsilon_{j,c}` are non-zero, the three values are assigned to
atom type J.  If all three :math:`\epsilon_i` values are zero, they are
ignored.  Thus the typical way to define the :math:`\epsilon_i` and
:math:`\epsilon_j` coefficients is to list their values in "pair_coeff
I J" commands when I = J, but set them to 0.0 when I != J.  If you do
list them when I != J, you should ensure they are consistent with their
values in other pair_coeff commands.

Note that if this potential is being used as a sub-style of
:doc:`pair_style hybrid <pair_hybrid>`, and there is no "pair_coeff I I"
setting made for RE-squared for a particular type I (because I-I
interactions are computed by another hybrid pair potential), then you
still need to ensure the epsilon a,b,c coefficients are assigned to
that type in a "pair_coeff I J" command.

For large uniform molecules it has been shown that the :math:`\epsilon_{*,*}`
energy parameters are approximately representable in terms of local
contact curvatures :ref:`(Everaers) <Everaers3>`:

.. math::

   \epsilon_a = \sigma \cdot { \frac{a}{ b \cdot c } }; \epsilon_b =
   \sigma \cdot { \frac{b}{ a \cdot c } }; \epsilon_c = \sigma \cdot {
   \frac{c}{ a \cdot b } }

where a, b, and c give the particle diameters.

The last coefficient is optional.  If not specified, the global cutoff
specified in the pair_style command is used.

----------

.. _resquared_details:

Implementation details
""""""""""""""""""""""

The following is a summary of the RE-squared formulas as implemented in
LAMMPS by Mike Brown (Sandia National Labs) in 2007.  Let the shape
matrices :math:`\mathbf{S}_i = \mathrm{diag}(a_i, b_i, c_i)` be given by
the radii of ellipsoid *i*, i.e. half of the diameters set with the
:doc:`set shape <set>` command or in the data file.  Let the relative
energy matrices :math:`\mathbf{E}_i = \mathrm{diag}(\epsilon_{i,a},
\epsilon_{i,b}, \epsilon_{i,c})` be given by the relative well depths of
the :doc:`pair_coeff <pair_coeff>` command.  Let :math:`\mathbf{A}_i` be
the rotation matrix from the simulation box frame to the body frame of
ellipsoid *i*, :math:`\mathbf{r}` the vector between the centers of the
two particles, :math:`r = |\mathbf{r}|`, and :math:`\hat{\mathbf{r}} =
\mathbf{r}/r`.

The RE-squared interaction between two ellipsoids is the sum of an
attractive and a repulsive part:

.. math::

   U = & U_A + U_R, \\
   U_\alpha = & \frac{A_{12}}{m_\alpha} \left( \frac{\sigma}{h}
   \right)^{n_\alpha} \left( 1 + o_\alpha \eta \chi \frac{\sigma}{h}
   \right) \prod_{i=1}^2 \frac{a_i b_i c_i}{(a_i + h/p_\alpha) (b_i +
   h/p_\alpha) (c_i + h/p_\alpha)}, \qquad \alpha = A, R

with the constants

.. math::

   m_A = -36, \quad n_A = 0, \quad o_A = 3, \quad p_A = 2, \qquad
   m_R = 2025, \quad n_R = 6, \quad o_R = \frac{45}{56}, \quad p_R =
   60^{1/3}.

The orientation-dependent terms are

.. math::

   \chi = & 2 \hat{\mathbf{r}}^T \mathbf{B}^{-1} \hat{\mathbf{r}},
   \qquad \mathbf{B} = \mathbf{A}_1^T \mathbf{E}_1 \mathbf{A}_1 +
   \mathbf{A}_2^T \mathbf{E}_2 \mathbf{A}_2, \\
   \eta = & \frac{\det(\mathbf{S}_1)/\sigma_1^2 +
   \det(\mathbf{S}_2)/\sigma_2^2}{\left[ \det(\mathbf{H}) / (\sigma_1 +
   \sigma_2) \right]^{1/2}}, \qquad \sigma_i = \left(
   \hat{\mathbf{r}}^T \mathbf{A}_i^T \mathbf{S}_i^{-2} \mathbf{A}_i
   \hat{\mathbf{r}} \right)^{-1/2}, \qquad \mathbf{H} =
   \frac{1}{\sigma_1} \mathbf{A}_1^T \mathbf{S}_1^2 \mathbf{A}_1 +
   \frac{1}{\sigma_2} \mathbf{A}_2^T \mathbf{S}_2^2 \mathbf{A}_2.

The distance of closest approach *h* is computed with the approximation
of :ref:`(Perram) <Perram2>`:

.. math::

   h = r - \sigma_{12}, \qquad \sigma_{12} = \left[ \frac{1}{2}
   \hat{\mathbf{r}}^T \mathbf{G}^{-1} \hat{\mathbf{r}} \right]^{-1/2},
   \qquad \mathbf{G} = \mathbf{A}_1^T \mathbf{S}_1^2 \mathbf{A}_1 +
   \mathbf{A}_2^T \mathbf{S}_2^2 \mathbf{A}_2.

The interaction between an ellipsoid and an LJ sphere is obtained from
the limit :math:`a_2 \to 0` with :math:`a_2 = b_2 = c_2` and
:math:`\frac{4}{3}\pi a_2^3 \rho = 1`:

.. math::

   U_\mathrm{elj} = & U_{A,\mathrm{elj}} + U_{R,\mathrm{elj}}, \\
   U_{\alpha,\mathrm{elj}} = & \frac{3 \sigma^3 p_\alpha^3}{4 \pi h^3}
   \frac{A_{12}}{m_\alpha} \left( \frac{\sigma}{h} \right)^{n_\alpha}
   \left( 1 + o_\alpha \chi_\mathrm{elj} \frac{\sigma}{h} \right)
   \frac{a_1 b_1 c_1}{(a_1 + h/p_\alpha) (b_1 + h/p_\alpha) (c_1 +
   h/p_\alpha)},

where :math:`A_{12}` has the modified form given above and
:math:`\chi_\mathrm{elj}` is computed as above, but with
:math:`\mathbf{B}_\mathrm{elj} = \mathbf{A}_1^T \mathbf{E}_1
\mathbf{A}_1 + \mathbf{I}`.  To reduce the error of the distance of
closest approach for two particles of very different size, the radii of
the ellipsoid are enlarged by the radius :math:`\sigma/2` of the LJ
sphere, :math:`\tilde{\mathbf{S}}_1 = \mathrm{diag}(a_1 + \sigma/2,
b_1 + \sigma/2, c_1 + \sigma/2)`, and the sphere is treated as a point:

.. math::

   h = r - \sigma_{12}, \qquad \sigma_{12} = \left[ \hat{\mathbf{r}}^T
   \left( \mathbf{A}_1^T \tilde{\mathbf{S}}_1^2 \mathbf{A}_1
   \right)^{-1} \hat{\mathbf{r}} \right]^{-1/2}.

Thus :math:`\sigma_{12}` is the distance from the center of the
ellipsoid to the surface of the enlarged ellipsoid along
:math:`\hat{\mathbf{r}}`.

The interaction between two LJ spheres is computed as

.. math::

   U_\mathrm{lj} = 4 A_{12} \left[ \left( \frac{\sigma}{r} \right)^{12}
   - \left( \frac{\sigma}{r} \right)^6 \right].

The forces and torques are computed from the analytic derivatives of
these expressions.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance can be mixed, but only for sphere pairs.  The
default mix value is *geometric*\ .  See the "pair_modify" command for
details.  Other type pairs cannot be mixed, due to the different
meanings of the energy prefactors used to calculate the interactions
and the implicit dependence of the ellipsoid-sphere interaction on the
equation for the Hamaker constant presented here.  Mixing of sigma and
epsilon followed by calculation of the energy prefactors using the
equations above is recommended.

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
*inner*, *middle*, *outer* keywords of the :doc:`run_style command <run_style>`.

----------

Restrictions
""""""""""""

This style is part of the ASPHERE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This pair style requires that atoms be ellipsoids as defined by the
:doc:`atom_style ellipsoid <atom_style>` command.

Particles acted on by the potential can be finite-size aspherical or
spherical particles, or point particles.  Spherical particles have all
3 of their shape parameters equal to each other.  Point particles have
all 3 of their shape parameters equal to 0.0.

The distance-of-closest-approach approximation used by LAMMPS becomes
less accurate when high-aspect ratio ellipsoids are used.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix nve/asphere <fix_nve_asphere>`,
:doc:`compute temp/asphere <compute_temp_asphere>`, :doc:`pair_style gayberne <pair_gayberne>`

Default
"""""""

none

----------

.. _Everaers3:

**(Everaers)** Everaers and Ejtehadi, Phys Rev E, 67, 041710 (2003).

.. _Babadi:

**(Babadi)** Babadi, Ejtehadi, Everaers, J Comp Phys, 219, 770-779 (2006).

.. _Perram2:

**(Perram)** Perram and Rasmussen, Phys Rev E, 54, 6565-6572 (1996).
