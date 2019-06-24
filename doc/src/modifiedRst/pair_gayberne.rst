.. index:: pair\_style gayberne

pair\_style gayberne command
============================

pair\_style gayberne/gpu command
================================

pair\_style gayberne/intel command
==================================

pair\_style gayberne/omp command
================================

Syntax
""""""


.. parsed-literal::

   pair_style gayberne gamma upsilon mu cutoff

* gamma = shift for potential minimum (typically 1)
* upsilon = exponent for eta orientation-dependent energy function
* mu = exponent for chi orientation-dependent energy function
* cutoff = global cutoff for interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style gayberne 1.0 1.0 1.0 10.0
   pair_coeff \* \* 1.0 1.7 1.7 3.4 3.4 1.0 1.0 1.0

Description
"""""""""""

The *gayberne* styles compute a Gay-Berne anisotropic LJ interaction
:ref:`(Berardi) <Berardi>` between pairs of ellipsoidal particles or an
ellipsoidal and spherical particle via the formulas

.. math::

   \mathbf{A}_1, \mathbf{A}_2, \mathbf{r}_{12}, \gamma ) \cdot \eta_{12} (
   \mathbf{A}_1, \mathbf{A}_2, \upsilon ) \cdot \chi_{12} ( \mathbf{A}_1,
   \mathbf{A}_2, \mathbf{r}_{12}, \mu ) $$
   

>>>image was here
   
   \end{document}
   :align: center

where A1 and A2 are the transformation matrices from the simulation
box frame to the body frame and r12 is the center to center vector
between the particles.  Ur controls the shifted distance dependent
interaction based on the distance of closest approach of the two
particles (h12) and the user-specified shift parameter gamma.  When
both particles are spherical, the formula reduces to the usual
Lennard-Jones interaction (see details below for when Gay-Berne treats
a particle as "spherical").

For large uniform molecules it has been shown that the energy
parameters are approximately representable in terms of local contact
curvatures :ref:`(Everaers) <Everaers2>`:

.. math::

   \sigma \cdot { \frac{b}{ a \cdot c } }; \epsilon_c = \sigma \cdot {
   \frac{c}{ a \cdot b } } $$
   
   \end{document}
   :align: center

The variable names utilized as potential parameters are for the most
part taken from :ref:`(Everaers) <Everaers2>` in order to be consistent with
the :doc:`RE-squared pair potential <pair_resquared>`.  Details on the
upsilon and mu parameters are given
`here <PDF/pair_resquared_extra.pdf>`_.

More details of the Gay-Berne formulation are given in the references
listed below and in `this supplementary document <PDF/pair_gayberne_extra.pdf>`_.

Use of this pair style requires the NVE, NVT, or NPT fixes with the
*asphere* extension (e.g. :doc:`fix nve/asphere <fix_nve_asphere>`) in
order to integrate particle rotation.  Additionally, :doc:`atom\_style ellipsoid <atom_style>` should be used since it defines the
rotational state and the size and shape of each ellipsoidal particle.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands, or by mixing as described below:

* epsilon = well depth (energy units)
* sigma = minimum effective particle radii (distance units)
* epsilon\_i_a = relative well depth of type I for side-to-side interactions
* epsilon\_i_b = relative well depth of type I for face-to-face interactions
* epsilon\_i_c = relative well depth of type I for end-to-end interactions
* epsilon\_j_a = relative well depth of type J for side-to-side interactions
* epsilon\_j_b = relative well depth of type J for face-to-face interactions
* epsilon\_j_c = relative well depth of type J for end-to-end interactions
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global
cutoff specified in the pair\_style command is used.

It is typical with the Gay-Berne potential to define *sigma* as the
minimum of the 3 shape diameters of the particles involved in an I,I
interaction, though this is not required.  Note that this is a
different meaning for *sigma* than the :doc:`pair\_style resquared <pair_resquared>` potential uses.

The epsilon\_i and epsilon\_j coefficients are actually defined for atom
types, not for pairs of atom types.  Thus, in a series of pair\_coeff
commands, they only need to be specified once for each atom type.

Specifically, if any of epsilon\_i_a, epsilon\_i_b, epsilon\_i_c are
non-zero, the three values are assigned to atom type I.  If all the
epsilon\_i values are zero, they are ignored.  If any of epsilon\_j_a,
epsilon\_j_b, epsilon\_j_c are non-zero, the three values are assigned
to atom type J.  If all three epsilon\_j values are zero, they are
ignored.  Thus the typical way to define the epsilon\_i and epsilon\_j
coefficients is to list their values in "pair\_coeff I J" commands when
I = J, but set them to 0.0 when I != J.  If you do list them when I !=
J, you should insure they are consistent with their values in other
pair\_coeff commands, since only the last setting will be in effect.

Note that if this potential is being used as a sub-style of
:doc:`pair\_style hybrid <pair_hybrid>`, and there is no "pair\_coeff I I"
setting made for Gay-Berne for a particular type I (because I-I
interactions are computed by another hybrid pair potential), then you
still need to insure the epsilon a,b,c coefficients are assigned to
that type. e.g. in a "pair\_coeff I J" command.

.. note::

   If the epsilon a = b = c for an atom type, and if the shape of
   the particle itself is spherical, meaning its 3 shape parameters are
   all the same, then the particle is treated as an LJ sphere by the
   Gay-Berne potential.  This is significant because if two LJ spheres
   interact, then the simple Lennard-Jones formula is used to compute
   their interaction energy/force using the specified epsilon and sigma
   as the standard LJ parameters.  This is much cheaper to compute than
   the full Gay-Berne formula.  To treat the particle as a LJ sphere with
   sigma = D, you should normally set epsilon a = b = c = 1.0, set the
   pair\_coeff sigma = D, and also set the 3 shape parameters for the
   particle to D.  The one exception is that if the 3 shape parameters
   are set to 0.0, which is a valid way in LAMMPS to specify a point
   particle, then the Gay-Berne potential will treat that as shape
   parameters of 1.0 (i.e. a LJ particle with sigma = 1), since it
   requires finite-size particles.  In this case you should still set the
   pair\_coeff sigma to 1.0 as well.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for this pair style can be mixed.  The default mix
value is *geometric*\ .  See the "pair\_modify" command for details.

This pair styles supports the :doc:`pair\_modify <pair_modify>` shift
option for the energy of the Lennard-Jones portion of the pair
interaction, but only for sphere-sphere interactions.  There is no
shifting performed for ellipsoidal interactions due to the anisotropic
dependence of the interaction.

The :doc:`pair\_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


The *gayberne* style is part of the ASPHERE package.  It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

These pair style require that atoms store torque and a quaternion to
represent their orientation, as defined by the
:doc:`atom\_style <atom_style>`.  It also require they store a per-type
:doc:`shape <set>`.  The particles cannot store a per-particle
diameter.

This pair style requires that atoms be ellipsoids as defined by the
:doc:`atom\_style ellipsoid <atom_style>` command.

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

:doc:`pair\_coeff <pair_coeff>`, :doc:`fix nve/asphere <fix_nve_asphere>`,
:doc:`compute temp/asphere <compute_temp_asphere>`, :doc:`pair\_style resquared <pair_resquared>`

**Default:** none


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


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
