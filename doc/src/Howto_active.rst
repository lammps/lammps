Active matter
=============

Active matter consists of particles that consume energy to propel
themselves, for example swimming bacteria, self-propelled colloids
("Janus particles"), or flocks of birds :ref:`(Marchetti2013)
<Marchetti2013>`, :ref:`(Bechinger2016) <Bechinger2016>`.  Since the
propulsion drives the system out of thermal equilibrium, active matter
shows behavior without an equilibrium counterpart, such as
motility-induced phase separation or collective motion ("flocking").
This page gives an overview of how common active matter models are set
up in LAMMPS with the commands of the :ref:`BROWNIAN <PKG-BROWNIAN>` and
:ref:`DIPOLE <PKG-DIPOLE>` packages.

Rather than providing a monolithic command for each model, LAMMPS
provides building blocks that are combined in the input script:

1. a per-particle **orientation** vector :math:`\mathbf{e}_i` (a unit
   vector), stored as the dipole moment of :doc:`atom_style dipole
   <atom_style>` or given by the quaternion of :doc:`atom_style
   ellipsoid <atom_style>`,
2. a **self-propulsion force** along the orientation (:doc:`fix
   propel/self <fix_propel_self>`), or a fluctuating active force that
   does not require an orientation (:doc:`fix propel/ou <fix_propel_ou>`),
3. the **dynamics of the orientation**: rotational diffusion (:doc:`fix
   brownian/sphere <fix_brownian>`, or rotational friction and noise
   from :doc:`fix langevin <fix_langevin>` with the *omega* keyword),
   random reorientations (:doc:`fix tumble <fix_tumble>`), and alignment
   torques (:doc:`fix align/self <fix_align_self>`, :doc:`fix
   align/neighbor <fix_align_neighbor>`),
4. the **translational dynamics**: overdamped (Brownian) dynamics with
   :doc:`fix brownian <fix_brownian>` and its variants, or inertial
   (Langevin) dynamics with :doc:`fix nve/sphere <fix_nve_sphere>` and
   :doc:`fix langevin <fix_langevin>`,
5. optionally, **interactions** between the particles through a
   :doc:`pair style <pair_style>`, typically a purely repulsive
   potential such as the WCA potential (:doc:`pair_style lj/cut
   <pair_lj>` with a cutoff of :math:`2^{1/6}\sigma` and
   :doc:`pair_modify shift yes <pair_modify>`).

The following sections describe how to combine these building blocks
into the most common models.  Complete input scripts are in the
``examples/PACKAGES/brownian`` directory.

Orientation vectors
-------------------

The dipole vector of :doc:`atom_style dipole <atom_style>` is the most
convenient way to represent the orientation of an active particle: it
is read from data files, stored in restart files, can be written to dump
files (as *mux*, *muy*, *muz*), and it is automatically communicated to
neighboring processors so that fixes and pair styles can use the
orientations of neighboring particles.  The fixes that apply torques
(:doc:`fix align/self <fix_align_self>`, :doc:`fix align/neighbor
<fix_align_neighbor>`) and :doc:`fix brownian/sphere <fix_brownian>`
require that the atom style also stores torques, which is the case for
:doc:`atom_style hybrid sphere dipole <atom_style>`.  The
self-propulsion fixes use the dipole vector as is, so it should be set
to a unit vector, e.g. with random orientations:

.. code-block:: LAMMPS

   atom_style      hybrid sphere dipole
   ...
   set             group all dipole/random 12345 1.0

In 2d simulations, the *dipole/random* keyword creates orientations in
the xy plane, and the fixes described here keep them in that plane.

Active Brownian particles
-------------------------

Active Brownian particles (ABPs) :ref:`(Romanczuk2012) <Romanczuk2012>`
move with a constant self-propulsion speed :math:`v_0` along their
orientation, which undergoes rotational diffusion with the rotational
diffusion coefficient :math:`D_r`.  In the overdamped limit,

.. math::

   \frac{d\mathbf{r}_i}{dt} = v_0 \mathbf{e}_i + \frac{1}{\gamma_t} \mathbf{F}_i + \sqrt{2 D_t}\, \mathbf{\xi}_i(t),
   \qquad
   \frac{d\mathbf{e}_i}{dt} = \sqrt{2 D_r}\, \mathbf{\eta}_i(t) \times \mathbf{e}_i

with :math:`\mathbf{F}_i` the force from interactions with other
particles, :math:`\gamma_t` the translational friction coefficient,
:math:`D_t = k_B T/\gamma_t` and :math:`D_r = k_B T/\gamma_r` the
translational and rotational diffusion coefficients, and
:math:`\mathbf{\xi}_i` and :math:`\mathbf{\eta}_i` Gaussian white noises.
This is obtained with :doc:`fix brownian/sphere <fix_brownian>` for the
overdamped translational and rotational dynamics and :doc:`fix
propel/self <fix_propel_self>` for the self-propulsion force :math:`f_P
= \gamma_t v_0`:

.. code-block:: LAMMPS

   units           lj
   dimension       2
   atom_style      hybrid sphere dipole
   ...
   set             group all dipole/random 12345 1.0

   pair_style      lj/cut 1.122462
   pair_coeff      * * 1.0 1.0
   pair_modify     shift yes

   fix             1 all brownian/sphere 1.0 12345 gamma_t 1.0 gamma_r 3.0
   fix             2 all propel/self dipole 40.0
   fix             3 all enforce2d

The same input works for 3d simulations after removing the
:doc:`dimension <dimension>` and :doc:`fix enforce2d <fix_enforce2d>`
commands.  The persistence length of the motion is :math:`v_0/D_r` and
the persistence time :math:`1/D_r` (in 2d); their ratios to the particle
size and the time to diffuse over a particle size define the Peclet
number, which is the key parameter of ABPs.

Alternatively, the inertial version of the model uses :doc:`fix
nve/sphere <fix_nve_sphere>` with the *update dipole* keyword to rotate
the orientations with the angular velocity of the particles, and
:doc:`fix langevin <fix_langevin>` with the *omega* keyword to provide
the translational and rotational friction and noise:

.. code-block:: LAMMPS

   fix             1 all nve/sphere update dipole
   fix             2 all langevin 1.0 1.0 0.1 12345 omega yes
   fix             3 all propel/self dipole 40.0
   fix             4 all enforce2d

Run-and-tumble particles
------------------------

Run-and-tumble particles model the motion of swimming bacteria: they
move at constant speed along their orientation (the "run") and at
random times, with an average rate :math:`\lambda`, pick a new random
orientation (the "tumble").  This is the ABP setup with :doc:`fix
tumble <fix_tumble>` in place of, or in addition to, the rotational
diffusion.  Without rotational diffusion, :doc:`fix brownian
<fix_brownian>` (which does not require torques, so :doc:`atom_style
dipole <atom_style>` is sufficient) can be used for the translational
motion:

.. code-block:: LAMMPS

   fix             1 all brownian 1.0 12345 gamma_t 1.0
   fix             2 all propel/self dipole 40.0
   fix             3 all tumble dipole 1.0 12345

For non-interacting particles, run-and-tumble particles with the tumble
rate :math:`\lambda = (d-1) D_r` have the same long-time diffusive
behavior as ABPs with the rotational diffusion coefficient :math:`D_r`
in :math:`d` dimensions :ref:`(Cates2013) <Cates2013>`.

Active Ornstein-Uhlenbeck particles
-----------------------------------

Active Ornstein-Uhlenbeck particles (AOUPs) are driven by an active
force whose Cartesian components are independent Gaussian processes with
an exponentially decaying time correlation :ref:`(Martin2021)
<Martin2021>`.  No orientation vector is needed, so this model works
with any atom style.  :doc:`fix propel/ou <fix_propel_ou>` adds the
active force with the given root mean square value per component and
correlation time:

.. code-block:: LAMMPS

   atom_style      atomic
   ...
   fix             1 all brownian 1.0 12345 gamma_t 1.0
   fix             2 all propel/ou 5.0 2.0 12345

Collective motion and alignment
-------------------------------

Collective motion arises when the orientations of neighboring particles
align.  :doc:`fix align/neighbor <fix_align_neighbor>` applies a torque
that rotates the orientation of a particle toward the orientations of
the other particles within a cutoff distance, either in parallel (polar
alignment, as in the Vicsek model) or in parallel or antiparallel
(nematic alignment, as for self-propelled rods).  The torques are turned
into rotations by :doc:`fix brownian/sphere <fix_brownian>` with the
rotational friction coefficient *gamma_r*, so that the alignment rate is
the ratio of the torque magnitude and *gamma_r*, and the noise of the
Vicsek model corresponds to the rotational diffusion:

.. code-block:: LAMMPS

   comm_modify     cutoff 2.0
   fix             1 all brownian/sphere 1.0 12345 gamma_t 1.0 gamma_r 1.0
   fix             2 all propel/self dipole 1.0
   fix             3 all align/neighbor dipole 2.0 1.5
   fix             4 all enforce2d

The alignment cutoff may be larger than the cutoff of the pair style,
in which case the communication cutoff has to be increased with the
:doc:`comm_modify cutoff <comm_modify>` command as shown above.

:doc:`fix align/self <fix_align_self>` provides an alternative mechanism
for collective motion, a torque that aligns the orientation of a
particle with its own velocity.  In dense systems this leads to
collective motion through the interactions between the particles.

Order parameters
----------------

The degree of collective motion is measured by the polar order
parameter, the length of the average orientation vector, which is 1 for
perfectly aligned particles and close to 0 for random orientations:

.. code-block:: LAMMPS

   compute         mu all property/atom mux muy muz
   compute         pol all reduce ave c_mu[1] c_mu[2] c_mu[3]
   variable        order equal sqrt(c_pol[1]^2+c_pol[2]^2+c_pol[3]^2)
   thermo_style    custom step v_order

In 2d, the corresponding nematic order parameter is

.. code-block:: LAMMPS

   variable        c2 atom 2.0*c_mu[1]^2-1.0
   variable        s2 atom 2.0*c_mu[1]*c_mu[2]
   compute         nem all reduce ave v_c2 v_s2
   variable        nematic equal sqrt(c_nem[1]^2+c_nem[2]^2)

The diffusive behavior of the particles can be measured with
:doc:`compute msd <compute_msd>`.  Note that a temperature is not well
defined for active systems and the kinetic energy of overdamped
particles is meaningless, so :doc:`compute pressure <compute_pressure>`
should be used with the *virial* keyword and without a temperature
compute.  Also see the documentation of :doc:`fix propel/self
<fix_propel_self>` for the "active pressure" contribution of the
self-propulsion forces to the virial.

----------

.. _Marchetti2013:

**(Marchetti2013)** M. C. Marchetti, J. F. Joanny, S. Ramaswamy, T. B. Liverpool, J. Prost, M. Rao, and R. A. Simha, Hydrodynamics of soft active matter, Rev. Mod. Phys. 85, 1143 (2013).

.. _Bechinger2016:

**(Bechinger2016)** C. Bechinger, R. Di Leonardo, H. Loewen, C. Reichhardt, G. Volpe, and G. Volpe, Active particles in complex and crowded environments, Rev. Mod. Phys. 88, 045006 (2016).

.. _Romanczuk2012:

**(Romanczuk2012)** P. Romanczuk, M. Baer, W. Ebeling, B. Lindner, and L. Schimansky-Geier, Active Brownian particles: From individual to collective stochastic dynamics, Eur. Phys. J. Special Topics 202, 1 (2012).
