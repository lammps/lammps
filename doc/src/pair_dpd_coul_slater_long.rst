.. index:: pair_style dpd/coul/slater/long
.. index:: pair_style dpd/coul/slater/long/gpu

pair_style dpd/coul/slater/long command
=======================================

Accelerator Variants: *dpd/coul/slater/long/gpu*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style dpd/coul/slater/long T cutoff_DPD seed lambda cutoff_coul

* T = temperature (temperature units)
* cutoff_DPD = global cutoff for DPD interactions (distance units)
* seed = random # seed (positive integer)
* lambda = decay length of the charge (distance units)
* cutoff_coul = global cutoff for Coulombic interactions (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style dpd/coul/slater/long 1.0 2.5 34387 0.25 3.0

   pair_coeff 1 1 78.0 4.5          # not charged by default
   pair_coeff 2 2 78.0 4.5 yes


Description
"""""""""""

.. versionadded:: 27June2024

Style *dpd/coul/slater/long* computes a force field for dissipative
particle dynamics (DPD) following the exposition in :ref:`(Groot)
<Groot5>`.  It also allows for the use of charged particles in the
model by adding a long-range Coulombic term to the DPD interactions.
The short-range portion of the Coulombics is calculated by this pair
style.  The long-range Coulombics are computed by use of the
:doc:`kspace_style <kspace_style>` command, e.g. using the Ewald or
PPPM styles.

Coulombic forces in mesoscopic models such as DPD employ potentials
without explicit excluded-volume interactions.  The goal is to prevent
artificial ionic pair formation by including a charge distribution in
the Coulomb potential, following the formulation in :ref:`(Melchor1)
<Melchor1>`.

.. note::

   This pair style is effectively the combination of the
   :doc:`pair_style dpd <pair_dpd>` and :doc:`pair_style
   coul/slater/long <pair_coul_slater>` commands, but should be more
   efficient (especially on GPUs) than using :doc:`pair_style
   hybrid/overlay dpd coul/slater/long <pair_hybrid>`. That is
   particularly true for the GPU package version of the pair style since
   this version is compatible with computing neighbor lists on the GPU
   instead of the CPU as is required for hybrid styles.

In the charged DPD model, the force on bead I due to bead J is given
as a sum of 4 terms:

.. math::

   \vec{f}  = & (F^C + F^D + F^R + F^E) \hat{r_{ij}} \\
   F^C      = & A w(r) \qquad \qquad \qquad \qquad \qquad r < r_{DPD} \\
   F^D      = & - \gamma w^2(r) (\hat{r_{ij}} \bullet \vec{v}_{ij}) \qquad \qquad r < r_{DPD} \\
   F^R      = & \sigma w(r) \alpha (\Delta t)^{-1/2} \qquad \qquad \qquad r < r_{DPD} \\
   w(r)     = & 1 - \frac{r}{r_{DPD}} \\
   F^E      = & \frac{C q_iq_j}{\epsilon r^2} \left( 1- exp\left( \frac{2r_{ij}}{\lambda} \right) \left( 1 + \frac{2r_{ij}}{\lambda} \left( 1 + \frac{r_{ij}}{\lambda} \right)\right) \right)

where :math:`F^C` is a conservative force, :math:`F^D` is a
dissipative force, :math:`F^R` is a random force, and :math:`F^E` is
an electrostatic force.  :math:`\hat{r_{ij}}` is a unit vector in the
direction :math:`r_i - r_j`, :math:`\vec{v}_{ij}` is the vector
difference in velocities of the two atoms :math:`\vec{v}_i -
\vec{v}_j`, :math:`\alpha` is a Gaussian random number with zero mean
and unit variance, *dt* is the timestep size, and :math:`w(r)` is a
weighting factor that varies between 0 and 1.

:math:`\sigma` is set equal to :math:`\sqrt{2 k_B T \gamma}`, where
:math:`k_B` is the Boltzmann constant and *T* is the temperature
parameter in the pair_style command.

:math:`r_{DPD}` is the pairwise cutoff for the first 3 DPD terms in
the formula as specified by *cutoff_DPD*.  For the :math:`F^E` term,
pairwise interactions within the specified *cutoff_coul* distance are
computed directly; interactions beyond that distance are computed in
reciprocal space.  *C* is the same Coulomb conversion factor used in
the Coulombic formulas described on the :doc:`pair_coul <pair_coul>`
doc page.

The following parameters must be defined for each pair of atoms types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* A (force units)
* :math:`\gamma` (force/velocity units)
* is_charged (optional boolean, default = no)

The *is_charged* parameter is optional and can be specified as *yes* or
*no*.  *Yes* should be used for interactions between two types of
charged particles.  *No* is the default and should be used for
interactions between two types of particles when one or both are
uncharged.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.  Note that the
user-specified random number seed is stored in the restart file, so
when a simulation is restarted, each processor will re-initialize its
random number generator the same way it did initially.  This means the
random forces will be random, but will not be the same as they would
have been if the original simulation had continued past the restart
time.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This style is part of the DPD-BASIC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

The default frequency for rebuilding neighbor lists is every 10 steps
(see the :doc:`neigh_modify <neigh_modify>` command). This may be too
infrequent since particles move rapidly and can overlap by large
amounts.  If this setting yields a non-zero number of "dangerous"
reneighborings (printed at the end of a simulation), you should
experiment with forcing reneighboring more often and see if system
energies/trajectories change.

This pair style requires use of the :doc:`comm_modify vel yes
<comm_modify>` command so that velocities are stored by ghost atoms.

This pair style also requires use of a long-range solvers from the
KSPACE package.

This pair style will not restart exactly when using the
:doc:`read_restart <read_restart>` command, though they should provide
statistically similar results.  This is because the forces they compute
depend on atom velocities.  See the :doc:`read_restart <read_restart>`
command for more details.

Related commands
""""""""""""""""

:doc:`pair_style dpd <pair_dpd>`, :doc:`pair_style coul/slater/long <pair_coul_slater>`,

Default
"""""""

For the pair_coeff command, the default is is_charged = no.

----------

.. _Groot5:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _Melchor1:

**(Melchor)** Gonzalez-Melchor, Mayoral, Velazquez, and Alejandre, J Chem Phys, 125, 224107 (2006).
