.. index:: pair_style bondval
.. index:: pair_style bondval/kk
.. index:: pair_style bondval/vec
.. index:: pair_style bondval/vec/kk


pair_style bondval command
==========================

Accelerator Variants: *bondval/kk*

pair_style bondval/vec command
==============================

Accelerator Variants: *bondval/vec/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style bondval power cutoff
   pair_style bondval/vec power cutoff

* *power* = exponent for the bond-valence force (typically 2.0)
* *cutoff* = global cutoff for bond-valence interactions (distance units)

Examples
""""""""

A three-element perovskite BaTiO3 system:

.. code-block:: LAMMPS

   pair_style hybrid/overlay lj/cut/coul/long 8.0 8.0 bondval 2.0 8.0 bondval/vec 2.0 8.0

   # lj/cut/coul/long parameters: epsilon sigma
   pair_coeff 1 1 lj/cut/coul/long 2.0 2.44805
   pair_coeff 1 2 lj/cut/coul/long 2.0 2.32592
   pair_coeff 1 3 lj/cut/coul/long 2.0 1.98792
   pair_coeff 2 2 lj/cut/coul/long 2.0 2.73825
   pair_coeff 2 3 lj/cut/coul/long 2.0 1.37741
   pair_coeff 3 3 lj/cut/coul/long 2.0 1.99269

   # bondval parameters: r0 alpha S V0
   pair_coeff 1 1 bondval 0.0000000 5.0000000 0.5973900 2.0000000
   pair_coeff 1 2 bondval 0.0000000 5.0000000 0.0000000 0.0000000
   pair_coeff 1 3 bondval 2.2900000 8.9400000 0.0000000 0.0000000
   pair_coeff 2 2 bondval 0.0000000 5.0000000 0.1653300 4.0000000
   pair_coeff 2 3 bondval 1.7980000 5.2000000 0.0000000 0.0000000
   pair_coeff 3 3 bondval 0.0000000 5.0000000 0.9306300 2.0000000

   # bondval/vec parameters: r0 alpha D W0
   pair_coeff 1 1 bondval/vec 0.0000000 5.0000000 0.0842900 0.1156100
   pair_coeff 1 2 bondval/vec 0.0000000 5.0000000 0.0000000 0.0000000
   pair_coeff 1 3 bondval/vec 2.1430000 8.9400000 0.0000000 0.0000000
   pair_coeff 2 2 bondval/vec 0.0000000 5.0000000 0.8248400 0.3943700
   pair_coeff 2 3 bondval/vec 1.7980000 5.2000000 0.0000000 0.0000000
   pair_coeff 3 3 bondval/vec 0.0000000 5.0000000 0.2800600 0.3165100

Description
"""""""""""

.. versionadded:: 4Jul2026

The bond-valence potential is an empirical potential based on the
conservation of the bond-valence (bondval) and bond-valence vector
(bondval/vec), fitted to DFT calculations for a given bulk semiconductor.
It is typically used with :doc:`pair_style hybrid/overlay <pair_hybrid>`
to combine Coulombic, Lennard-Jones repulsion, bond-valence, and
bond-valence vector contributions.

The bond-valence for a given atom pair (:math:`V_{ij}`) is a measure
of the bonding strength calculated from the length of the bond
(:math:`r_{ij}`) by:

.. math::

   V_{ij} = \left( \frac{r0_{ij}}{r_{ij}} \right)^{\alpha_{ij}}

where :math:`r0_{ij}` and :math:`\alpha_{ij}` are Brown's empirical
parameters for bond-valence.  The bond-valence vector is a vector lying
along the bond between atom i and j:
:math:`\vec{V}_{ij} = V_{ij} \hat{R}_{ij}`, where
:math:`\hat{R}_{ij}` is the unit vector pointing from atom i toward atom j.

The total potential energy of the system consists of a Coulombic energy
(:math:`E_c`), a short-range repulsive energy (:math:`E_r`), a
bond-valence energy (:math:`E_{BV}`), a bond-valence vector energy
(:math:`E_{BVV}`), and an optional angle potential (:math:`E_a`):

.. math::

   E = E_c + E_r + E_{BV} + E_{BVV} + E_a

.. math::

   E_c = \sum_{i<j} \frac{q_i q_j}{r_{ij}}

.. math::

   E_r = \sum_{i<j} \left(\frac{B_{ij}}{r_{ij}}\right)^{12}

.. math::

   E_{BV} = \sum_i S_i \left(V_i - V_{0,i}\right)^2

.. math::

   E_{BVV} = \sum_i D_i \left(|\vec{W}_i|^2-|\vec{W}_{0,i}|^2\right)^2

.. math::

   E_a = k \sum_i^{N_{\mathrm{oxygen}}} \left(\theta_i - 180^{\circ}\right)^2

where :math:`V_i = \sum_{j \neq i} V_{ij}` is the bond-valence sum
(BVS), and :math:`\vec{W}_i = \sum_{j \neq i} \vec{V}_{ij}` is
the bond-valence vector sum (BVVS), :math:`q_i` is the ionic
charge, :math:`B_{ij}` is the short-range repulsion parameter,
:math:`S_i` and :math:`D_i` are scaling parameters (energy units),
:math:`k` is the angle spring constant, and :math:`\theta_i` is the
O-O-O angle along the common axis of two adjacent oxygen octahedra.

The pair style *bondval* computes :math:`E_{BV}`, the bond-valence
energy term.  The pair style *bondval/vec* computes :math:`E_{BVV}`, the
bond-valence vector energy term.  The remaining energy contributions
(:math:`E_c` and :math:`E_r`) are typically provided by :doc:`pair_style
lj/cut/coul/long <pair_lj_cut_coul>`, and :math:`E_a` by an
:doc:`angle_style <angle_style>`.  The *power* argument to the
pair_style command specifies the exponent used in computing the forces
and is usually set to 2.0.

The quantities :math:`r_{ij}`, :math:`V_i`, and :math:`\vec{W}_i`
are computed at each timestep from the current atom positions.  All
other parameters must be provided by the user.

The potential parameters are obtained by training against DFT energies
for a given system.  Trained and tested parameters for several bulk
perovskite oxide materials are published by the group of Andrew M. Rappe
(see references below).  The published parameters use LAMMPS *metal*
units.  Atom charges are also fitted parameters; they can be set via
:doc:`read_data <read_data>` or :doc:`set <set>`.  The angle potential,
if used, can be set with :doc:`angle_style harmonic <angle_harmonic>`
together with :doc:`angle_coeff <angle_coeff>`.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands.  When used within :doc:`pair_style hybrid/overlay <pair_hybrid>`,
the pair style name must be included as shown in the examples.

For *bondval*:

* :math:`r0_{ij}` = first of Brown's empirical bond-valence parameter (distance units)
* :math:`\alpha_{ij}` = second of Brown's empirical bond-valence exponent (dimensionless)
* :math:`S_i` = penalty for deviating from ideal bond valence (energy units)
* :math:`V_{0,i}` = ideal bond-valence for a given atom type (dimensionless)
* cutoff (distance units) -- optional

The first two parameters (:math:`r0_{ij}` and :math:`\alpha_{ij}`) are
pair atom-type dependent and contribute to the bond-valence calculation for all
atom-type pairs.  The penalty constant :math:`S_i` and ideal value
:math:`V_{0,i}` are single atom-type dependent. Thus only same-species values
(I = J) are nonzero in input file.

For *bondval/vec*:

* :math:`r0_{ij}` = first of Brown's empirical bond-valence parameter (distance units)
* :math:`\alpha_{ij}` = second of Brown's empirical bond-valence exponent (dimensionless)
* :math:`D_i` = penalty for deviating from ideal bond valence vector (energy units)
* :math:`W_{0,i}` = ideal bond-valence vector for a given atom type (dimensionless)
* cutoff (distance units) -- optional

The same distinction applies: :math:`r0_{ij}` and :math:`\alpha_{ij}`
are pair atom-type dependent and used for all pairs, while :math:`D_i` and
:math:`W_{0,i}` are atom-type dependent. Thus only same-species values
(I = J) are nonzero in input file.


The final cutoff coefficient is optional for both styles.  If not
specified, the global cutoff given in the pair_style command is used.

---------

.. include:: accel_styles.rst

---------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These pair styles do not support parameter mixing.  Coefficients must
be specified explicitly for all atom type pairs.

These pair styles do not support the :doc:`pair_modify <pair_modify>`
shift, table, or tail options.

These pair styles write their information to :doc:`binary restart
files <restart>`, so pair_style and pair_coeff commands do not need to
be specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

Restrictions
""""""""""""

These pair styles are part of the EXTRA-PAIR package.  They are only
enabled if LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

These pair styles must be used with
:doc:`pair_style hybrid/overlay <pair_hybrid>`.  They cannot be used
as standalone pair styles.

For a physically correct simulation, *bondval*, *bondval/vec*, and a
:doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>` contribution
for :math:`E_c` and :math:`E_r` must all be combined via
``hybrid/overlay``.  The published parameters for this potential are
fitted to only include :math:`r^{-12}` repulsion term
(:math:`E_r`) in the Lennard-Jones potential, while the attractive :math:`r^{-6}`
contribution is set to 0.
To run with the published parameters correctly, users must manually initialize
the internal variables ``lj2`` and ``lj4`` in the source code of
:doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>` to be zero
in order to remove the attractive :math:`r^{-6}` contribution.

Related commands
""""""""""""""""

* :doc:`pair_style hybrid/overlay <pair_hybrid>`
* :doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`
* :doc:`angle_style harmonic <angle_harmonic>`
* :doc:`kspace_style <kspace_style>`

Default
"""""""

none


----------

**(Grinberg)** I. Grinberg, V. R. Cooper, and A. M. Rappe, *Nature*, **419**, 909 (2002).

**(Shin1)** Y.-H. Shin, J.-Y. Son, B.-J. Lee, I. Grinberg, and A. M. Rappe, *J. Phys.: Condens. Matter*, **20**, 015224 (2008).

**(Shin2)** Y.-H. Shin, V. R. Cooper, I. Grinberg, and A. M. Rappe, *Phys. Rev. B*, **71**, 054104 (2005).

**(Liu1)** S. Liu, I. Grinberg, and A. M. Rappe, *J. Phys.: Condens. Matter*, **25**, 102202 (2013).

**(Liu2)** S. Liu, I. Grinberg, H. Takenaka, and A. M. Rappe, *Phys. Rev. B*, **88**, 104102 (2013).

**(Brown1)** I. D. Brown, *Chem. Rev.*, **109**, 6858 (2009).

**(Brown2)** I. Brown and R. Shannon, *Acta Crystallogr. A*, **29**, 266 (1973).

**(Brown3)** I. Brown and K. K. Wu, *Acta Crystallogr. B*, **32**, 1957 (1976).
