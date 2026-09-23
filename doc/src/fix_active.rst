.. index:: fix active

fix active command
==================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID active model keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* active = style name of this fix command
* model = *abp* or *aoup* or *chiral* or *rtp* or *vicsek* or *rod* or *spinner* or *iabp* or *mips* or *nematic*
* zero or more keyword/value pairs may be appended
* keyword = *v0* or *Dr* or *Dt* or *gamma* or *kT* or *seed* or *dimension* or *omega* or *tumble_rate* or *Rcut* or *alignment* or *aspect* or *omega_spin* or *eta_odd* or *tau* or *Da* or *activity* or *K* or *rho_max*

  .. parsed-literal::

       *v0* value = self-propulsion speed (distance/time units)
       *Dr* value = rotational diffusion coefficient (1/time units)
       *Dt* value = translational diffusion coefficient (distance\ :sup:`2`/time units), only used by model *abp*
       *gamma* value = friction coefficient that scales the active force (mass/time units)
       *kT* value = thermal energy (energy units), only used by models *rod* and *iabp*
       *seed* value = integer seed for the random number generator
       *dimension* value = 2 or 3, overrides the dimension of the simulation for the initial orientations and the density normalization of model *mips*
       *omega* value = intrinsic angular velocity (1/time units), only used by model *chiral*
       *tumble_rate* value = tumbling rate (1/time units), only used by model *rtp*
       *Rcut* value = cutoff for alignment or density estimation (distance units), used by models *vicsek*, *mips*, and *nematic*
       *alignment* value = weight of the neighbor orientations (unitless), only used by model *vicsek*
       *aspect* value = ratio of perpendicular to parallel friction (unitless), only used by model *rod*
       *omega_spin* value = spinning frequency (1/time units), only used by model *spinner*
       *eta_odd* value = odd viscosity coefficient (mass/time units), only used by model *spinner*
       *tau* value = persistence time (time units), only used by model *aoup*
       *Da* value = variance of the active velocity (distance\ :sup:`2`/time\ :sup:`2` units), only used by model *aoup*
       *activity* value = propulsion speed (distance/time units), only used by model *nematic*
       *K* value = alignment stiffness (1/time units), only used by model *nematic*
       *rho_max* value = density at which the propulsion speed vanishes (1/distance\ :sup:`2` units), only used by model *mips*

Examples
""""""""

.. code-block:: LAMMPS

   # overdamped active Brownian particles
   fix 1 all nve
   fix 2 all langevin 0.01 0.01 1.0 48279
   fix 3 all active abp v0 1.0 Dr 0.1 gamma 1.0 seed 54321
   fix 4 all enforce2d

   # Vicsek alignment, neighbor list with a 1.5 distance cutoff
   pair_style zero 1.5
   pair_coeff * *
   fix 1 all nve
   fix 2 all langevin 0.0 0.0 1.0 48279
   fix 3 all active vicsek v0 1.0 Dr 0.5 Rcut 1.0 alignment 1.0 seed 54321
   fix 4 all enforce2d

   # inertial active Brownian particles (self-contained drag and noise)
   fix 1 all nve
   fix 2 all active iabp v0 1.0 Dr 0.1 gamma 10.0 kT 0.1 seed 54321
   fix 3 all enforce2d

Description
"""""""""""

Add a self-propulsion force to each atom in the group. Ten different
models of active matter can be selected with the *model* argument, so
that different microscopic mechanisms of self-propulsion can be compared
within the same input script by changing a single keyword. The
:doc:`fix propel/self <fix_propel_self>` command provides a
self-propulsion force along a dipole, velocity, or quaternion direction
for particles that already carry that vector; this fix instead stores
its own orientation for each atom and evolves it according to the
selected model. A detailed description of all models, together with
validation against analytical results and further benchmarks, is given
in :ref:`(Chand) <active-Chand>`.

Each atom in the group carries an orientation angle :math:`\theta_i`,
which defines the unit vector :math:`\hat{e}_i = (\cos\theta_i,
\sin\theta_i)` along which the atom is propelled. The angles are
initialized to random values uniformly distributed between 0 and
:math:`2\pi` when the fix is defined. Except for model *spinner*, which
has no orientation, the orientation is updated once per time step with
the explicit Euler-Maruyama scheme. The force is added in the
post-force stage of the time step and is therefore added to the forces
of the pair styles and of all other fixes that were already applied.
The position and velocity of the atoms are not changed by this fix but
by the time integrator that is used together with it.

The parameter :math:`\gamma` (keyword *gamma*) is the friction
coefficient that converts the self-propulsion speed into a force,
:math:`F = \gamma v_0`. For the models that only add a force, it should
be set to the friction coefficient of the thermostat that supplies the
drag, so that the terminal velocity of a free particle is :math:`v_0`.
For :doc:`fix langevin <fix_langevin>`, this friction coefficient is
:math:`\gamma = m/T_{damp}`, where :math:`m` is the mass of the atom and
:math:`T_{damp}` is the damping parameter of that command.

The models are described in the following. All random numbers
:math:`\xi` are drawn from a Gaussian distribution with zero mean and
unit variance, and :math:`\Delta t` is the time step.

----------

Model *abp* implements active Brownian particles :ref:`(Howse) <active-Howse>`.
The orientation performs rotational diffusion

.. math::

   \theta_i(t+\Delta t) = \theta_i(t) + \sqrt{2 D_r \Delta t}\, \xi_i

and the force :math:`\vec{F}_i = \gamma v_0 \hat{e}_i` is added. If
*Dt* is larger than zero, a Gaussian white-noise force with the
amplitude :math:`\gamma \sqrt{2 D_t/\Delta t}` is added to each
Cartesian component as well. This should be left at zero when the
translational noise is already provided by the thermostat, for example
by :doc:`fix langevin <fix_langevin>` with the temperature
:math:`k_B T` and the friction coefficient :math:`\gamma`, which
corresponds to :math:`D_t = k_B T/\gamma`. In two dimensions the
long-time effective diffusion coefficient is
:math:`D_{eff} = D_t + v_0^2/(2 D_r)`.

Model *aoup* implements active Ornstein-Uhlenbeck particles
:ref:`(Fodor) <active-Fodor>`. Instead of an orientation, each atom
carries an active velocity :math:`\vec{v}^a_i` that is updated with the
exact solution of the Ornstein-Uhlenbeck process

.. math::

   \vec{v}^a_i(t+\Delta t) = e^{-\Delta t/\tau}\, \vec{v}^a_i(t) + \sqrt{D_a \left(1 - e^{-2\Delta t/\tau}\right)}\, \vec{\xi}_i

and the force :math:`\vec{F}_i = \gamma \vec{v}^a_i` is added. The
stationary variance of each component of the active velocity is
:math:`D_a`, and the initial active velocities are drawn from that
distribution. The update is exact for any time step. In two dimensions
the long-time effective diffusion coefficient is :math:`D_{eff} = D_t +
D_a \tau`.

Model *chiral* implements chiral active Brownian particles
:ref:`(Liebchen) <active-Liebchen>`, which move on circles because their
orientation rotates with the constant angular velocity :math:`\omega`
in addition to the rotational diffusion

.. math::

   \theta_i(t+\Delta t) = \theta_i(t) + \omega \Delta t + \sqrt{2 D_r \Delta t}\, \xi_i

and the force :math:`\vec{F}_i = \gamma v_0 \hat{e}_i` is added.

Model *rtp* implements run-and-tumble particles
:ref:`(Tailleur) <active-Tailleur>`. The orientation is constant except
for tumbles. In every time step, an atom tumbles with the probability
:math:`1 - \exp(-\alpha \Delta t)`, where :math:`\alpha` is the value of
the *tumble_rate* keyword, and its new orientation is drawn uniformly
between 0 and :math:`2\pi`. The keyword *Dr* has no effect on this
model. The force :math:`\vec{F}_i = \gamma v_0 \hat{e}_i` is added.

Model *vicsek* implements Vicsek-type alignment
:ref:`(Vicsek) <active-Vicsek>`. The new orientation of an atom is the
direction of the vector

.. math::

   (1-a)\, \hat{e}_i + a \sum_{j \in B_i} \hat{e}_j

with the addition of the angular noise :math:`\sqrt{2 D_r \Delta t}\,
\xi_i`. Here :math:`a` is the value of the *alignment* keyword and
:math:`B_i` is the set that contains atom *i* itself and all its
neighbors closer than *Rcut*. For :math:`a = 1` this is the direction of
the circular mean of the orientations in the neighborhood, which is the
original Vicsek rule. All atoms are updated synchronously, using the
orientations from the start of the time step. The force
:math:`\vec{F}_i = \gamma v_0 \hat{e}_i` is added.

Model *rod* implements active rods with an anisotropic friction
:ref:`(Wensink) <active-Wensink>`. The friction coefficients along and
perpendicular to the rod axis are :math:`\gamma_\parallel =
\gamma/\text{aspect}` and :math:`\gamma_\perp = \gamma`. The
orientation is rotated by the torque :math:`\tau_z = \hat{e}_i \times
\vec{F}^{ext}_i` of the force :math:`\vec{F}^{ext}_i` that is acting on
the atom when this fix is applied, i.e. the forces from the pair styles
and from the fixes that were applied before

.. math::

   \theta_i(t+\Delta t) = \theta_i(t) + \frac{\tau_z}{\gamma_\perp} \Delta t + \sqrt{2 D_r \Delta t}\, \xi_i

Then the force

.. math::

   \vec{F}_i = -\gamma_\parallel v_\parallel \hat{e}_i - \gamma_\perp \vec{v}_\perp + \gamma_\parallel v_0 \hat{e}_i + \sqrt{2 \gamma_\parallel k_B T/\Delta t}\, \xi_\parallel \hat{e}_i + \sqrt{2 \gamma_\perp k_B T/\Delta t}\, \xi_\perp \hat{e}_\perp

is added, where :math:`v_\parallel` and :math:`\vec{v}_\perp` are the
components of the current velocity of the atom along and perpendicular
to :math:`\hat{e}_i`. This model supplies its own drag and thermal noise.

Model *spinner* implements spinning particles with an odd (Hall)
viscosity :ref:`(Banerjee) <active-Banerjee>`. This model has no
orientation and does not use random numbers. It adds a force that is
perpendicular to the velocity of the atom

.. math::

   \vec{F}_i = \eta_{odd}\, \omega_{spin}\, (\hat{z} \times \vec{v}_i) = \eta_{odd}\, \omega_{spin}\, (-v_{y,i},\, v_{x,i})

The friction is not part of this model and has to be supplied by a
thermostat. If a constant external force is applied in the presence of a
friction coefficient :math:`\gamma`, the velocity of the atom is rotated
away from the direction of that force by the Hall angle :math:`\phi`
with :math:`\tan\phi = \eta_{odd}\, \omega_{spin}/\gamma`.

Model *iabp* implements inertial active Brownian particles
:ref:`(Loewen) <active-Loewen>`, for which the momentum of the particle
does not relax instantaneously. The orientation performs rotational
diffusion as in model *abp*, and the force

.. math::

   \vec{F}_i = -\gamma \vec{v}_i + \gamma v_0 \hat{e}_i + \sqrt{2 \gamma k_B T/\Delta t}\, \vec{\xi}_i

is added. This model supplies its own drag and thermal noise. The Stokes
number :math:`m D_r/\gamma` measures the importance of the inertia.

Model *mips* implements a quorum-sensing model of motility-induced phase
separation :ref:`(Cates) <active-Cates>`, in which the propulsion speed
decreases with the local density :math:`\rho_i`

.. math::

   v(\rho_i) = v_0 \max\left(1 - \rho_i/\rho_{max},\, 0\right)

The local density is estimated by counting atom *i* and all its
neighbors closer than *Rcut*, divided by the area :math:`\pi R_{cut}^2`
(or by the volume :math:`4\pi R_{cut}^3/3` if the *dimension* keyword is
set to 3). The orientation performs rotational diffusion as in model
*abp* and the force :math:`\vec{F}_i = \gamma v(\rho_i) \hat{e}_i` is
added. Steric repulsion between the atoms is not part of this model and
can be added with any pair style. Note that motility-induced phase
separation can also occur for the other models when they are combined
with a repulsive pair style at a sufficiently large Peclet number and
density.

Model *nematic* implements a dry active nematic
:ref:`(Doostmohammadi) <active-Doostmohammadi>`. The self-propulsion of
each atom is polar, but the alignment interaction is apolar, i.e.
invariant under :math:`\theta_i \rightarrow \theta_i + \pi`. In the
neighborhood :math:`B_i` (atom *i* and all its neighbors closer than
*Rcut*) the local nematic tensor is computed

.. math::

   Q_{xx} = \frac{1}{N_i} \sum_{j \in B_i} \left(\cos^2\theta_j - \tfrac{1}{2}\right), \qquad Q_{xy} = \frac{1}{N_i} \sum_{j \in B_i} \cos\theta_j \sin\theta_j

with :math:`N_i` the number of atoms in :math:`B_i`. The local director
angle is :math:`\phi_i = \tfrac{1}{2}\, \mathrm{atan2}(2 Q_{xy}, 2 Q_{xx})`
and the orientation is updated according to

.. math::

   \theta_i(t+\Delta t) = \theta_i(t) - K \sin\left(2(\theta_i - \phi_i)\right) \Delta t + \sqrt{2 D_r \Delta t}\, \xi_i

Then the force :math:`\vec{F}_i = \gamma\, \text{activity}\, \hat{e}_i`
is added. All atoms are updated synchronously, using the orientations
from the start of the time step.

----------

The following table lists the keywords that are used by each model.
Keywords that are valid but are not used by the selected model are
accepted and ignored, which allows one to run the same input script with
different models. An unknown keyword is an error.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - model
     - keywords that are used
   * - *abp*
     - *v0*, *Dr*, *Dt*, *gamma*, *seed*
   * - *aoup*
     - *Da*, *tau*, *gamma*, *seed*
   * - *chiral*
     - *v0*, *Dr*, *omega*, *gamma*, *seed*
   * - *rtp*
     - *v0*, *tumble_rate*, *gamma*, *seed*
   * - *vicsek*
     - *v0*, *Dr*, *Rcut*, *alignment*, *gamma*, *seed*
   * - *rod*
     - *v0*, *Dr*, *aspect*, *kT*, *gamma*, *seed*
   * - *spinner*
     - *omega_spin*, *eta_odd*
   * - *iabp*
     - *v0*, *Dr*, *gamma*, *kT*, *seed*
   * - *mips*
     - *v0*, *Dr*, *Rcut*, *rho_max*, *gamma*, *dimension*, *seed*
   * - *nematic*
     - *activity*, *K*, *Rcut*, *Dr*, *gamma*, *seed*

Combination with other fixes
""""""""""""""""""""""""""""

This fix only adds forces. The positions and velocities of the atoms are
integrated by another fix.

For the models *abp*, *aoup*, *chiral*, *rtp*, *vicsek*, *spinner*,
*mips*, and *nematic*, use :doc:`fix nve <fix_nve>` together with
:doc:`fix langevin <fix_langevin>`, which provides the friction and the
thermal noise. A small damping parameter of :doc:`fix langevin
<fix_langevin>` compared to the orientational persistence time gives the
overdamped limit that is described by the equations of motion of
Brownian particles.

.. note::

   The models *rod* and *iabp* add their own friction and thermal noise.
   They must be used with :doc:`fix nve <fix_nve>` only. If they are
   combined with :doc:`fix langevin <fix_langevin>` or :doc:`fix
   brownian <fix_brownian>`, the friction and the noise are applied
   twice and the dynamics of the system are wrong. No warning is issued.

For simulations in two dimensions, :doc:`fix enforce2d <fix_enforce2d>`
should be used as well.

Parallel simulations
""""""""""""""""""""

The orientations (and the active velocities of model *aoup*) are stored
per atom and are moved together with the atoms when they migrate between
processors or when LAMMPS sorts the atoms. For the models *vicsek* and
*nematic*, the orientations of the ghost atoms are communicated at every
time step, so that the alignment with neighbors on other processors is
correct for any domain decomposition.

Each processor uses its own random number generator, which is seeded
with the value of the *seed* keyword plus the rank of the processor. The
random numbers are consumed in the order of the atoms on each processor.
The trajectories of a simulation are therefore reproducible when the same
input script is run on the same number of processors, but they are **not**
independent of the number of processors: the same input script run on a
different number of processors gives a different random sequence and thus
different trajectories, even though all statistical properties, such as
the mean-square displacement or the order parameter, are the same. This
is also the case for :doc:`fix langevin <fix_langevin>`.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`. The orientations are not restarted; they are initialized
randomly whenever the fix is defined.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix. No global or per-atom quantities are stored by this fix for access
by various :doc:`output commands <Howto_output>`. No parameter of this
fix can be used with the *start/stop* keywords of the :doc:`run <run>`
command. This fix is not invoked during :doc:`energy minimization
<minimize>`. The active force is not included in the virial and thus
not in the pressure.

Restrictions
""""""""""""

This fix is part of the BROWNIAN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

The active force is only applied in the xy-plane. If the fix is used in
a three-dimensional simulation, a warning is printed, the force has no
z-component, and the initial orientation of all atoms is along the
x-axis.

The models *vicsek*, *mips*, and *nematic* use a full neighbor list. Only
atoms that are closer than the neighbor cutoff (cutoff of the pair style
plus the skin distance) are found. If *Rcut* is larger than that cutoff,
a warning is printed. When the atoms do not interact through a pair
style, :doc:`pair_style zero <pair_zero>` can be used to set the
neighbor cutoff.

Related commands
""""""""""""""""

:doc:`fix propel/self <fix_propel_self>`,
:doc:`fix align/self <fix_align_self>`,
:doc:`fix brownian <fix_brownian>`,
:doc:`fix langevin <fix_langevin>`,
:doc:`fix nve <fix_nve>`,
:doc:`fix enforce2d <fix_enforce2d>`

Default
"""""""

The keyword defaults are *v0* = 1.0, *Dr* = 0.1, *Dt* = 0.0, *gamma* =
1.0, *kT* = 1.0, *seed* = 12345, *dimension* = the dimension of the
simulation, *omega* = 0.0, *tumble_rate* = 0.1, *Rcut* = 1.0,
*alignment* = 1.0, *aspect* = 2.0, *omega_spin* = 1.0, *eta_odd* = 0.0,
*tau* = 1.0, *Da* = 1.0, *activity* = 1.0, *K* = 1.0, and *rho_max* =
1.0.

----------

.. _active-Howse:

**(Howse)** J. R. Howse, R. A. L. Jones, A. J. Ryan, T. Gough,
R. Vafabakhsh, and R. Golestanian, Phys. Rev. Lett. 99, 048102 (2007).

.. _active-Fodor:

**(Fodor)** E. Fodor, C. Nardini, M. E. Cates, J. Tailleur, P. Visco, and
F. van Wijland, Phys. Rev. Lett. 117, 038103 (2016).

.. _active-Liebchen:

**(Liebchen)** B. Liebchen and H. Loewen, J. Chem. Phys. 157, 090901 (2022).

.. _active-Tailleur:

**(Tailleur)** J. Tailleur and M. E. Cates, Phys. Rev. Lett. 100, 218103 (2008).

.. _active-Vicsek:

**(Vicsek)** T. Vicsek, A. Czirok, E. Ben-Jacob, I. Cohen, and O. Shochet,
Phys. Rev. Lett. 75, 1226 (1995).

.. _active-Wensink:

**(Wensink)** H. H. Wensink, J. Dunkel, S. Heidenreich, K. Drescher,
R. E. Goldstein, H. Loewen, and J. M. Yeomans, Proc. Natl. Acad. Sci.
U.S.A. 109, 14308 (2012).

.. _active-Banerjee:

**(Banerjee)** D. Banerjee, A. Souslov, A. G. Abanov, and V. Vitelli,
Nat. Commun. 8, 1573 (2017).

.. _active-Loewen:

**(Loewen)** H. Loewen, J. Chem. Phys. 152, 040901 (2020).

.. _active-Cates:

**(Cates)** M. E. Cates and J. Tailleur, Annu. Rev. Condens. Matter
Phys. 6, 219 (2015).

.. _active-Doostmohammadi:

**(Doostmohammadi)** A. Doostmohammadi, J. Ignes-Mullol, J. M. Yeomans,
and F. Sagues, Nat. Commun. 9, 3246 (2018).

.. _active-Chand:

**(Chand)** R. Chand and S. Bukhari, Modelling Simul. Mater. Sci. Eng.
(2026), https://doi.org/10.1088/1361-651X/aeaa0a
