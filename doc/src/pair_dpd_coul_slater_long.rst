.. index:: pair_style dpd/coul/slater/long
.. index:: pair_style dpd/coul/slater/long/gpu

pair_style dpd/coul/slater/long command
======================

Accelerator Variants: *dpd/coul/slater/long/gpu*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style dpd/coul/slater/long T cutoff_DPD seed lambda cutoff_coul

   pair_coeff I J a_IJ Gamma is_charged

* T = temperature (temperature units) (dpd only)
* cutoff_DPD = global cutoff for DPD interactions (distance units)
* seed = random # seed (positive integer)
* lambda = decay length of the charge (distance units)
* cutoff_coul = real part cutoff for Coulombic interactions (distance units)
* I,J = numeric atom types, or type labels
* Gamma = DPD Gamma coefficient
* is_charged (boolean) set to yes if I and J are charged beads

Examples
""""""""

.. code-block:: LAMMPS

   pair_style dpd/coul/slater/long 1.0 2.5 34387 0.25 3.0
   pair_coeff 1 1 78.0 4.5 # not charged by default
   pair_coeff 2 2 78.0 4.5 yes


Description
"""""""""""

Style *dpd/coul/slater/long* computes a force field for dissipative particle dynamics
(DPD) following the exposition in :ref:`(Groot) <Groot1>` with the addition of
electrostatic interactions. The coulombic forces in mesoscopic models
employ potentials without explicit excluded-volume interactions.
The goal is to prevent artificial ionic pair formation by including a charge
distribution in the Coulomb potential, following the formulation of
:ref:`(Melchor) <Melchor>`:

The force on bead I due to bead J is given as a sum
of 4 terms

.. math::

   \vec{f}  = & (F^C + F^D + F^R + F^E) \hat{r_{ij}} \\
   F^C      = & A w(r) \qquad \qquad \qquad \qquad \qquad r < r_c \\
   F^D      = & - \gamma w^2(r) (\hat{r_{ij}} \bullet \vec{v}_{ij}) \qquad \qquad r < r_c \\
   F^R      = & \sigma w(r) \alpha (\Delta t)^{-1/2} \qquad \qquad \qquad r < r_c \\
   w(r)     = & 1 - \frac{r}{r_c} \\
   F^E      = & \frac{Cq_iq_j}{\epsilon r^2} \left( 1- exp\left( \frac{2r_{ij}}{\lambda} \right) \left( 1 + \frac{2r_{ij}}{\lambda} \left( 1 + \frac{r_{ij}}{\lambda} \right)\right) \right)

where :math:`F^C` is a conservative force, :math:`F^D` is a dissipative
force, :math:`F^R` is a random force, and :math:`F^E` is an electrostatic force.
:math:`\hat{r_{ij}}` is a unit vector in the direction
:math:`r_i - r_j`, :math:`\vec{v}_{ij}` is
the vector difference in velocities of the two atoms :math:`\vec{v}_i -
\vec{v}_j`, :math:`\alpha` is a Gaussian random number with zero mean
and unit variance, *dt* is the timestep size, and :math:`w(r)` is a
weighting factor that varies between 0 and 1.  :math:`r_c` is the
pairwise cutoff.  :math:`\sigma` is set equal to :math:`\sqrt{2 k_B T
\gamma}`, where :math:`k_B` is the Boltzmann constant and *T* is the
temperature parameter in the pair_style command.
C is the same Coulomb conversion factor as in the pair_styles
coul/cut and coul/long. In this way the Coulomb
interaction between ions is corrected at small distances r, and
the long-range interactions are computed either by the Ewald or the PPPM technique.


The following parameters must be defined for each
pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in
the examples above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* A (force units)
* :math:`\gamma` (force/velocity units)
* is_charged (boolean)


.. note::

   This style is the combination of :doc:`pair_style dpd <pair_dpd>` and :doc:`pair_style coul/slater/long <pair_coul_slater>`.

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
user-specified random number seed is stored in the restart file, so when
a simulation is restarted, each processor will re-initialize its random
number generator the same way it did initially.  This means the random
forces will be random, but will not be the same as they would have been
if the original simulation had continued past the restart time.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.


----------

Restrictions
""""""""""""

This style is part of the DPD-BASIC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

The default frequency for rebuilding neighbor lists is every 10 steps
(see the :doc:`neigh_modify <neigh_modify>` command). This may be too
infrequent since particles move rapidly and
can overlap by large amounts.  If this setting yields a non-zero number
of "dangerous" reneighborings (printed at the end of a simulation), you
should experiment with forcing reneighboring more often and see if
system energies/trajectories change.

This pair style requires you to use the :doc:`comm_modify vel yes
<comm_modify>` command so that velocities are stored by ghost atoms.

This pair style also requires the long-range solvers included in the KSPACE package.


This pair style will not restart exactly when using the
:doc:`read_restart <read_restart>` command, though they should provide
statistically similar results.  This is because the forces they compute
depend on atom velocities.  See the :doc:`read_restart <read_restart>`
command for more details.

Related commands
""""""""""""""""

:doc:`pair_style dpd <pair_dpd>`, :doc:`pair_style coul/slater/long <pair_coul_slater>`,
:doc:`pair_coeff <pair_coeff>`, :doc:`fix nvt <fix_nh>`, :doc:`fix langevin <fix_langevin>`,
:doc:`pair_style srp <pair_srp>`, :doc:`fix mvv/dpd <fix_mvv_dpd>`.

Default
"""""""

is_charged = no

----------

.. _Groot1:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _Melchor:

**(Melchor)** Gonzalez-Melchor, Mayoral, Velazquez, and Alejandre, J Chem Phys, 125, 224107 (2006).
