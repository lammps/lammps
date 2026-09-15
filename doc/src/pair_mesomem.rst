.. index:: pair_style mesomem


pair_style mesomem command
======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style mesomem cutoff


* cutoff = global cutoff for interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style   mesomem  2.5
   pair_coeff   *  *  1.0  1.0  15.0  1.0  2.5  2.0  5  0.0


Description
"""""""""""

.. versionadded:: 15Dec2026

The *mesomem* style computes an anisotropic interaction
between pairs of coarse-grained particles considering the relative particle orientations. This potential is  developed as a particle-based solvent-free model for biological membranes :ref:`(Sillano2026) <_Sillano>`. Unlike :doc:`pair_style ylz <pair_ylz>`, where the angular dependence multiplies the isotropic radial potential, this pair style adds an independent tilt/splay energy, weighted by a smooth radial window function, on top of a purely radial isotropic potential. The total pair energy is:

.. math::

   U ( \mathbf{r}_{ij}, \mathbf{n}_i, \mathbf{n}_j ) = U_{iso}(r) + w(r) \left[ U_{tilt} ( \mathbf{r}_{ij}, \mathbf{n}_i, \mathbf{n}_j ) + U_{splay} ( \mathbf{n}_i, \mathbf{n}_j ) \right]

where :math:`\mathbf{r}_{i}` and :math:`\mathbf{r}_{j}` are the center position vectors of particles i and j, respectively, :math:`\mathbf{r}_{ij}=\mathbf{r}_{i}-\mathbf{r}_{j}` is the inter-particle distance vector, :math:`r=\left|\mathbf{r}_{ij} \right|` and :math:`{\hat{\mathbf{r}}}_{ij}=\mathbf{r}_{ij}/r`.  The unit vectors
:math:`\mathbf{n}_{i}` and :math:`\mathbf{n}_{j}` represent the axes of symmetry of particles i and j, respectively, derived from the per-particle dipole vector.

The isotropic part :math:`U_{iso}` is purely radial:

.. math::

   U_{iso}(r) = \left\{\begin{matrix}
      \epsilon \left [ \left ( \dfrac{\sigma}{r} \right )^{4} - 2 \left ( \dfrac{\sigma}{r} \right )^{2} \right ], & r < \sigma \\[1.2ex]
      -\epsilon \; \cos^{2\zeta} \left [ \dfrac{\pi}{2} \dfrac{r - \sigma}{r_{c} - \sigma} \right ], & \sigma \le r < r_{c} \\[1.2ex]
      0, & r \ge r_{c}
   \end{matrix}\right.

where :math:`\epsilon` is the well depth, :math:`\sigma` is the distance that minimizes :math:`U_{iso}(r)`,
:math:`r_{c}` is the (per type pair) isotropic cutoff, and :math:`\zeta` controls the slope of the attractive branch and hence the diffusivity of the particles in the in-plane direction of the membrane.

The tilt and splay energies act only over a shorter-ranged,
orientation-dependent cutoff :math:`w_{c}` and are switched off smoothly as :math:`r \to w_{c}` by the weight function

.. math::

   w(r) = \left\{\begin{matrix}
      \exp \left [ - \dfrac{r^{2}}{ \left ( w_{c}/2 \right )^{2} \left ( 1 - \left ( r / w_{c} \right )^{4} \right ) } \right ], & r < w_{c} \\[1.2ex]
      0, & r \ge w_{c}
   \end{matrix}\right.

Defining the spontaneous-curvature shift :math:`s = \tfrac{1}{2} r C_{0}` and the per-particle tilt deviations
:math:`d_i = \mathbf{n}_i \cdot \hat{\mathbf{r}}_{ij} + s` and :math:`d_j = \mathbf{n}_j \cdot \hat{\mathbf{r}}_{ij} - s`, the tilt and splay energies are

.. math::

   U_{tilt} = \frac{1}{2} k_{tilt} \left ( d_i^{2} + d_j^{2} \right ) \\\\
   U_{splay} = \frac{1}{2} k_{splay} \left ( \mathbf{n}_i \cdot \mathbf{n}_j - 1 + 2 s^{2} \right )^{2}

where :math:`k_{tilt}` and :math:`k_{splay}` are the tilt and splay spring constants and :math:`C_{0}` is the spontaneous curvature parameter of the type pair.  Setting :math:`C_{0}=0` recovers the standard tilt/splay form (:math:`d_i = \mathbf{n}_i \cdot \hat{\mathbf{r}}_{ij}`, :math:`d_j = \mathbf{n}_j \cdot \hat{\mathbf{r}}_{ij}`, and :math:`U_{splay} \propto (\mathbf{n}_i \cdot \mathbf{n}_j - 1)^{2}`).

Use of this pair style requires the NVE, NVT, or NPT fixes with the *sphere* extension (e.g. :doc:`fix nve/sphere <fix_nve_sphere>`) in order to integrate particle rotation.  Additionally, :doc:`atom_style hybrid dipole sphere <atom_style>` should be used since it defines the orientation of each particle.

The following coefficients must be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\sigma` = minimum effective particle radii, :math:`r_{min}` in :math:`U_{iso}` (distance units)
* :math:`\epsilon` = well depth of the isotropic potential (energy units)
* :math:`k_{tilt}` = tilt spring constant (energy units)
* :math:`k_{splay}` = splay spring constant (energy units)
* :math:`r_{c}` = isotropic cutoff for this type pair (distance units)
* :math:`w_{c}` = orientation (tilt/splay) cutoff for this type pair (distance units)
* :math:`\zeta` = tuning parameter for the slope of the attractive branch
* :math:`C_{0}` = spontaneous curvature parameter (inverse distance units)

All eight coefficients must always be specified; none of them is optional.
:math:`w_{c}` must not be larger than either :math:`r_{c}` or the global
cutoff specified in the pair_style command.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Coefficients for all type pairs, including I == J and I != J, must be set explicitly via the :doc:`pair_coeff <pair_coeff>` command; an unset type pair is a fatal error.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need to be specified in an input script that reads a restart file.


----------

Restrictions
""""""""""""

The *mesomem* style is part of the DIPOLE package.  It is only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This pair style requires that atoms store torque and a dipole vector to represent their orientation, as defined by the :doc:`atom_style <atom_style>`.

This pair style requires that **all** atoms are hybrid dipole sphere style as defined by the :doc:`atom_style hybrid <atom_style>` command.


Related commands
""""""""""""""""
:doc:`pair_coeff <pair_coeff>`,
:doc:`pair_style ylz <pair_ylz>`,

Default
"""""""

none

----------

.. _Sillano:

**(Sillano2026)** Sillano, Marrink, Idema, Phys. Rev. E 2026.
