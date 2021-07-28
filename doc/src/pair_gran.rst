.. index:: pair_style gran/hooke
.. index:: pair_style gran/hooke/omp
.. index:: pair_style gran/hooke/history
.. index:: pair_style gran/hooke/history/omp
.. index:: pair_style gran/hooke/history/kk
.. index:: pair_style gran/hertz/history
.. index:: pair_style gran/hertz/history/omp

pair_style gran/hooke command
=============================

Accelerator Variants: *gran/hooke/omp*

pair_style gran/hooke/history command
=====================================

Accelerator Variants: *gran/hooke/history/omp*, *gran/hooke/history/kk*

pair_style gran/hertz/history command
=====================================

Accelerator Variants: *gran/hertz/history/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style Kn Kt gamma_n gamma_t xmu dampflag keyword

* style = *gran/hooke* or *gran/hooke/history* or *gran/hertz/history*
* Kn = elastic constant for normal particle repulsion (force/distance units or pressure units - see discussion below)
* Kt = elastic constant for tangential contact (force/distance units or pressure units - see discussion below)
* gamma_n = damping coefficient for collisions in normal direction (1/time units or 1/time-distance units - see discussion below)
* gamma_t = damping coefficient for collisions in tangential direction (1/time units or 1/time-distance units - see discussion below)
* xmu = static yield criterion (unitless value between 0.0 and 1.0e4)
* dampflag = 0 or 1 if tangential damping force is excluded or included

* keyword = *limit_damping*

  .. parsed-literal::

      *limit_damping* value = none
         limit damping to prevent attractive interaction

.. note::

   Versions of LAMMPS before 9Jan09 had different style names for
   granular force fields.  This is to emphasize the fact that the
   Hertzian equation has changed to model polydispersity more accurately.
   A side effect of the change is that the Kn, Kt, gamma_n, and gamma_t
   coefficients in the pair_style command must be specified with
   different values in order to reproduce calculations made with earlier
   versions of LAMMPS, even for monodisperse systems.  See the NOTE below
   for details.

Examples
""""""""

.. code-block:: LAMMPS

   pair_style gran/hooke/history 200000.0 NULL 50.0 NULL 0.5 1
   pair_style gran/hooke 200000.0 70000.0 50.0 30.0 0.5 0
   pair_style gran/hooke 200000.0 70000.0 50.0 30.0 0.5 0 limit_damping


Description
"""""""""""

The *gran* styles use the following formulas for the frictional force
between two granular particles, as described in
:ref:`(Brilliantov) <Brilliantov>`, :ref:`(Silbert) <Silbert>`, and
:ref:`(Zhang) <Zhang3>`, when the distance r between two particles of radii
Ri and Rj is less than their contact distance d = Ri + Rj.  There is
no force between the particles when r > d.

The two Hookean styles use this formula:

.. math::

   F_{hk} = (k_n \delta \mathbf{n}_{ij} -
   m_{eff} \gamma_n\mathbf{ v}_n) -
   (k_t \mathbf{ \Delta s}_t +
   m_{eff} \gamma_t \mathbf{v}_t)

The Hertzian style uses this formula:

.. math::

   F_{hz} = \sqrt{\delta} \sqrt{\frac{R_i R_j}{R_i + R_j}} F_{hk} =
     \sqrt{\delta} \sqrt{\frac{R_i R_j}{R_i + R_j}}
     \Big[ (k_n \delta \mathbf{n}_{ij} -
       m_{eff} \: \gamma_n \mathbf{ v}_n) -
       (k_t \mathbf{ \Delta s}_t +
       m_{eff} \: \gamma_t \mathbf{v}_t) \Big]

In both equations the first parenthesized term is the normal force
between the two particles and the second parenthesized term is the
tangential force.  The normal force has 2 terms, a contact force and a
damping force.  The tangential force also has 2 terms: a shear force
and a damping force.  The shear force is a "history" effect that
accounts for the tangential displacement between the particles for the
duration of the time they are in contact.  This term is included in
pair styles *hooke/history* and *hertz/history*, but is not included
in pair style *hooke*\ .  The tangential damping force term is included
in all three pair styles if *dampflag* is set to 1; it is not included
if *dampflag* is set to 0.

The other quantities in the equations are as follows:

* :math:`\delta` = d - r = overlap distance of 2 particles
* :math:`K_n` = elastic constant for normal contact
* :math:`K_t` = elastic constant for tangential contact
* :math:`\gamma_n` = viscoelastic damping constant for normal contact
* :math:`\gamma_t` = viscoelastic damping constant for tangential contact
* :math:`m_{eff} = M_i M_j / (M_i + M_j) =` effective mass of 2 particles of mass M_i and M_j
* :math:`\mathbf{\Delta s}_t =` tangential displacement vector between 2 particles       which is truncated to satisfy a frictional yield criterion
* :math:`n_{ij} =` unit vector along the line connecting the centers of the 2 particles
* :math:`V_n =` normal component of the relative velocity of the 2 particles
* :math:`V_t =` tangential component of the relative velocity of the 2 particles

The :math:`K_n`, :math:`K_t`, :math:`\gamma_n`, and :math:`\gamma_t`
coefficients are specified as parameters to the pair_style command.  If
a NULL is used for :math:`K_t`, then a default value is used where
:math:`K_t = 2/7 K_n`.  If a NULL is used for :math:`\gamma_t`, then a
default value is used where :math:`\gamma_t = 1/2 \gamma_n`.

The interpretation and units for these 4 coefficients are different in
the Hookean versus Hertzian equations.

The Hookean model is one where the normal push-back force for two
overlapping particles is a linear function of the overlap distance.
Thus the specified :math:`K_n` is in units of (force/distance).  Note
that this push-back force is independent of absolute particle size (in
the monodisperse case) and of the relative sizes of the two particles
(in the polydisperse case).  This model also applies to the other terms
in the force equation so that the specified :math:`\gamma_n` is in units
of (1/time), :math:`K_t` is in units of (force/distance), and
:math:`\gamma_t` is in units of (1/time).

The Hertzian model is one where the normal push-back force for two
overlapping particles is proportional to the area of overlap of the
two particles, and is thus a non-linear function of overlap distance.
Thus Kn has units of force per area and is thus specified in units of
(pressure).  The effects of absolute particle size (monodispersity)
and relative size (polydispersity) are captured in the radii-dependent
pre-factors.  When these pre-factors are carried through to the other
terms in the force equation it means that the specified :math:`\gamma_n` is in
units of (1/(time\*distance)), :math:`K_t` is in units of (pressure), and
:math:`\gamma_t` is in units of (1/(time\*distance)).

Note that in the Hookean case, :math:`K_n` can be thought of as a linear
spring constant with units of force/distance.  In the Hertzian case,
:math:`K_n` is like a non-linear spring constant with units of
force/area or pressure, and as shown in the :ref:`(Zhang) <Zhang3>`
paper, :math:`K_n = 4G / (3(1-\nu))` where :math:`\nu =` the Poisson ratio,
G = shear modulus = :math:`E / (2(1+\nu))`, and E = Young's modulus.  Similarly,
:math:`K_t = 4G / (2-\nu)`.  (NOTE: in an earlier version of the manual, we incorrectly
stated that :math:`K_t = 8G / (2-\nu)`.)

Thus in the Hertzian case :math:`K_n` and :math:`K_t` can be set to
values that corresponds to properties of the material being modeled.
This is also true in the Hookean case, except that a spring constant
must be chosen that is appropriate for the absolute size of particles in
the model.  Since relative particle sizes are not accounted for, the
Hookean styles may not be a suitable model for polydisperse systems.

.. note::

   In versions of LAMMPS before 9Jan09, the equation for Hertzian
   interactions did not include the :math:`\sqrt{r_i r_j / (r_i + r_j)}`
   term and thus was not as accurate for polydisperse systems.  For
   monodisperse systems, :math:`\sqrt{ r_i r_j /(r_i+r_j)}` is a
   constant factor that effectively scales all 4 coefficients:
   :math:`K_n, K_t, \gamma_n, \gamma_t`.  Thus you can set the values of
   these 4 coefficients appropriately in the current code to reproduce
   the results of a previous Hertzian monodisperse calculation.  For
   example, for the common case of a monodisperse system with particles
   of diameter 1, all 4 of these coefficients should now be set 2x
   larger than they were previously.

Xmu is also specified in the pair_style command and is the upper limit
of the tangential force through the Coulomb criterion Ft = xmu\*Fn,
where Ft and Fn are the total tangential and normal force components
in the formulas above.  Thus in the Hookean case, the tangential force
between 2 particles grows according to a tangential spring and
dash-pot model until Ft/Fn = xmu and is then held at Ft = Fn\*xmu until
the particles lose contact.  In the Hertzian case, a similar analogy
holds, though the spring is no longer linear.

.. note::

   Normally, xmu should be specified as a fractional value between
   0.0 and 1.0, however LAMMPS allows large values (up to 1.0e4) to allow
   for modeling of systems which can sustain very large tangential
   forces.

The effective mass *m_eff* is given by the formula above for two
isolated particles.  If either particle is part of a rigid body, its
mass is replaced by the mass of the rigid body in the formula above.
This is determined by searching for a :doc:`fix rigid <fix_rigid>`
command (or its variants).

For granular styles there are no additional coefficients to set for
each pair of atom types via the :doc:`pair_coeff <pair_coeff>` command.
All settings are global and are made via the pair_style command.
However you must still use the :doc:`pair_coeff <pair_coeff>` for all
pairs of granular atom types.  For example the command

.. code-block:: LAMMPS

   pair_coeff * *

should be used if all atoms in the simulation interact via a granular
potential (i.e. one of the pair styles above is used).  If a granular
potential is used as a sub-style of :doc:`pair_style hybrid <pair_hybrid>`, then specific atom types can be used in the
pair_coeff command to determine which atoms interact via a granular
potential.

If two particles are moving away from each other while in contact, there
is a possibility that the particles could experience an effective attractive
force due to damping. If the *limit_damping* keyword is used, this option
will zero out the normal component of the force if there is an effective
attractive force.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for granular pair styles.

These pair styles write their information to :doc:`binary restart files <restart>`, so a pair_style command does not need to be
specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

The single() function of these pair styles returns 0.0 for the energy
of a pairwise interaction, since energy is not conserved in these
dissipative potentials.  It also returns only the normal component of
the pairwise interaction force.  However, the single() function also
calculates 10 extra pairwise quantities.  The first 3 are the
components of the tangential force between particles I and J, acting
on particle I.  The fourth is the magnitude of this tangential force.
The next 3 (5-7) are the components of the relative velocity in the
normal direction (along the line joining the 2 sphere centers).  The
last 3 (8-10) the components of the relative velocity in the
tangential direction.

These extra quantities can be accessed by the :doc:`compute pair/local <compute_pair_local>` command, as *p1*, *p2*, ...,
*p10*\ .

----------

Restrictions
""""""""""""

All the granular pair styles are part of the GRANULAR package.  It is
only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

These pair styles require that atoms store torque and angular velocity
(omega) as defined by the :doc:`atom_style <atom_style>`.  They also
require a per-particle radius is stored.  The *sphere* atom style does
all of this.

This pair style requires you to use the :doc:`comm_modify vel yes <comm_modify>` command so that velocities are stored by ghost
atoms.

These pair styles will not restart exactly when using the
:doc:`read_restart <read_restart>` command, though they should provide
statistically similar results.  This is because the forces they
compute depend on atom velocities.  See the
:doc:`read_restart <read_restart>` command for more details.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none

----------

.. _Brilliantov:

**(Brilliantov)** Brilliantov, Spahn, Hertzsch, Poschel, Phys Rev E, 53,
p 5382-5392 (1996).

.. _Silbert:

**(Silbert)** Silbert, Ertas, Grest, Halsey, Levine, Plimpton, Phys Rev
E, 64, p 051302 (2001).

.. _Zhang3:

**(Zhang)** Zhang and Makse, Phys Rev E, 72, p 011301 (2005).
