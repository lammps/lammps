.. index:: pair_style dpd/ext
.. index:: pair_style dpd/ext/tstat

pair_style dpd/ext command
==========================

pair_style dpd/ext/tstat command
================================

Syntax
""""""


.. code-block:: LAMMPS

   pair_style dpd/ext T cutoff seed
   pair_style dpd/ext/tstat Tstart Tstop cutoff seed

* T = temperature (temperature units)
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* cutoff = global cutoff for DPD interactions (distance units)
* seed = random # seed (positive integer)

Examples
""""""""


.. code-block:: LAMMPS

   pair_style dpd/ext 1.0 2.5 34387
   pair_coeff 1 1 25.0 4.5 4.5 0.5 0.5 1.2
   pair_coeff 1 2 40.0 4.5 4.5 0.5 0.5 1.2

   pair_style dpd/ext/tstat 1.0 1.0 2.5 34387
   pair_coeff 1 1 4.5 4.5 0.5 0.5 1.2
   pair_coeff 1 2 4.5 4.5 0.5 0.5 1.2

Description
"""""""""""

The style *dpd/ext* computes an extended force field for dissipative
particle dynamics (DPD) following the exposition in :ref:`(Groot)
<Groot>`, :ref:`(Junghans) <Junghans>`.

Style *dpd/ext/tstat* invokes an extended DPD thermostat on pairwise
interactions, equivalent to the non-conservative portion of the extended
DPD force field. To use *dpd/ext/tstat* as a thermostat for another pair
style, use the :doc:`pair_style hybrid/overlay <pair_hybrid>` command to
compute both the desired pair interaction and the thermostat for each
pair of particles.

For the style *dpd/ext*, the force on atom I due to atom J is given as
a sum of 3 terms

.. math::

   \mathbf{f}  = & f^C + f^D + f^R \qquad \qquad r < r_c \\
   f^C      = & A_{ij} w(r) \hat{\mathbf{r}}_{ij} \\
   f^D      = & - \gamma_{\parallel} w_{\parallel}^2(r) (\hat{\mathbf{r}}_{ij} \cdot \mathbf{v}_{ij}) \hat{\mathbf{r}}_{ij}  - \gamma_{\perp} w_{\perp}^2 (r) ( \mathbf{I} - \hat{\mathbf{r}}_{ij} \hat{\mathbf{r}}_{ij}^{\rm T} ) \mathbf{v}_{ij} \\
   f^R      = & \sigma_{\parallel} w_{\parallel}(r) \frac{\alpha}{\sqrt{\Delta t}} \hat{\mathbf{r}}_{ij}  + \sigma_{\perp} w_{\perp} (r) ( \mathbf{I} - \hat{\mathbf{r}}_{ij} \hat{\mathbf{r}}_{ij}^{\rm T} ) \frac{\mathbf{\xi}_{ij}}{\sqrt{\Delta t}}\\
   w(r)     = & 1 - r/r_c \\

where :math:`\mathbf{f}^C` is a conservative force, :math:`\mathbf{f}^D`
is a dissipative force, and :math:`\mathbf{f}^R` is a random
force. :math:`A_{ij}` is the maximum repulsion between the two atoms,
:math:`\hat{\mathbf{r}}_{ij}` is a unit vector in the direction
:math:`\mathbf{r}_i - \mathbf{r}_j`, :math:`\mathbf{v}_{ij} =
\mathbf{v}_i - \mathbf{v}_j` is the vector difference in velocities of
the two atoms, :math:`\alpha` and :math:`\mathbf{\xi}_{ij}` are Gaussian
random numbers with zero mean and unit variance, :math:`\Delta t` is the
timestep, :math:`w (r) = 1 - r / r_c` is a weight function for the
conservative interactions that varies between 0 and 1, :math:`r_c` is
the corresponding cutoff, :math:`w_{\alpha} ( r ) = ( 1 - r / \bar{r}_c
)^{s_{\alpha}}`, :math:`\alpha \equiv ( \parallel, \perp )`, are weight
functions with coefficients :math:`s_\alpha` that vary between 0 and 1,
:math:`\bar{r}_c` is the corresponding cutoff, :math:`\mathbf{I}` is the
unit matrix, :math:`\sigma_{\alpha} = \sqrt{2 k T \gamma_{\alpha}}`,
where :math:`k` is the Boltzmann constant and :math:`T` is the
temperature in the pair\_style command.

For the style *dpd/ext/tstat*, the force on atom I due to atom J is
the same as the above equation, except that the conservative
:math:`\mathbf{f}^C` term is dropped. Also, during the run, T is set
each timestep to a ramped value from Tstart to Tstop.

For the style *dpd/ext*, the pairwise energy associated with style
*dpd/ext* is only due to the conservative force term
:math:`\mathbf{f}^C`, and is shifted to be zero at the cutoff distance
:math:`r_c`. The pairwise virial is calculated using all three
terms. There is no pairwise energy for style *dpd/ext/tstat*, but the
last two terms of the formula contribute the virial.

For the style *dpd/ext/tstat*, the force on atom I due to atom J is the
same as the above equation, except that the conservative
:math:`\mathbf{f}^C` term is dropped.  Also, during the run, T is set
each timestep to a ramped value from Tstart to Tstop.

For the style *dpd/ext*, the pairwise energy associated with style
*dpd/ext* is only due to the conservative force term
:math:`\mathbf{f}^C`, and is shifted to be zero at the cutoff distance
:math:`r_c`. The pairwise virial is calculated using all three
terms. There is no pairwise energy for style *dpd/ext/tstat*, but the
last two terms of the formula contribute the virial.

For the style *dpd/ext*, the following coefficients must be defined for
each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command
as in the examples above:

* A (force units)
* :math:`\gamma_{\perp}` (force/velocity units)
* :math:`\gamma_{\parallel}` (force/velocity units)
* :math:`s_{\perp}` (unitless)
* :math:`s_{\parallel}` (unitless)
* :math:`r_c` (distance units)

The last coefficient is optional. If not specified, the global DPD
cutoff is used. Note that :math:`\sigma`'s are set equal to
:math:`\sqrt{2 k T \gamma}`, where :math:`T` is the temperature set by
the :doc:`pair_style <pair_style>` command so it does not need to be
specified.

For the style *dpd/ext/tstat*, the coefficients defined for each pair of
atoms types via the :doc:`pair_coeff <pair_coeff>` command is the same,
except that A is not included.

.. note::

   If you are modeling DPD polymer chains, you may want to use the
   :doc:`pair_style srp <pair_srp>` command in conjunction with these pair
   styles. It is a soft segmental repulsive potential (SRP) that can
   prevent DPD polymer chains from crossing each other.

.. note::

   The virial calculation for pressure when using this pair style includes
   all the components of force listed above, including the random force.

----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

The style *dpd/ext* does not support mixing. Thus, coefficients for all
I,J pairs must be specified explicitly.

The pair styles do not support the :doc:`pair_modify <pair_modify>`
shift option for the energy of the pair interaction. Note that as
discussed above, the energy due to the conservative :math:`\mathbf{f}^C`
term is already shifted to be zero at the cutoff distance :math:`r_c`.

The :doc:`pair_modify <pair_modify>` table option is not relevant for
the style *dpd/ext*.

The style *dpd/ext* does not support the :doc:`pair_modify
<pair_modify>` tail option for adding long-range tail corrections to
energy and pressure.

The pair styles can only be used via the pair keyword of the
:doc:`run_style respa <run_style>` command. They do not support the
*inner*, *middle*, and *outer*\ keywords.

The style *dpd/ext/tstat* can ramp its target temperature over multiple
runs, using the start and stop keywords of the :doc:`run <run>`
command. See the :doc:`run <run>` command for details of how to do this.

----------


Restrictions
""""""""""""

These styles are part of the DPD-BASIC package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

The default frequency for rebuilding neighbor lists is every 10 steps
(see the :doc:`neigh_modify <neigh_modify>` command). This may be too
infrequent for style *dpd/ext* simulations since particles move rapidly
and can overlap by large amounts. If this setting yields a non-zero
number of \say{dangerous} reneighborings (printed at the end of a
simulation), you should experiment with forcing reneighboring more often
and see if system energies/trajectories change.

The pair styles require to use the :doc:`comm_modify vel yes
<comm_modify>` command so that velocities are stored by ghost atoms.

The pair styles will not restart exactly when using the
:doc:`read_restart <read_restart>` command, though they should provide
statistically similar results. This is because the forces they compute
depend on atom velocities. See the :doc:`read_restart <read_restart>`
command for more details.

Related commands
""""""""""""""""

:doc:`pair_style dpd <pair_dpd>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`fix nvt <fix_nh>`, :doc:`fix langevin <fix_langevin>`,
:doc:`pair_style srp <pair_srp>`

**Default:** none

----------

.. _Groot:

**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _Junghans:

**(Junghans)** Junghans, Praprotnik and Kremer, Soft Matter 4, 156, 1119-1128 (2008).
