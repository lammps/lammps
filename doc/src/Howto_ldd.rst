Local Density Dependent Potentials via LDD
==========================================

.. contents::
   :local:

--------------------

Overview
--------

The :doc:`LDD pair style <pair_ldd>` in the BOCS package implements
potentials that are a function of the local density,
:math:`\rho_{\beta|I}`, of type :math:`\beta` particles around a central
particle, :math:`I`.  There are two types of local density dependent
potentials implemented.  There are Local Density (LD) potentials
:math:`U_{LD}` , which are direct functions of :math:`\rho_{\beta|I}`.
There are also Square Gradient (SG) potentials :math:`U_{SG}`.  SG
potentials use a direct function of the local density,
:math:`U_\nabla(\rho_{\beta|I})`, as a coefficient of the square
gradient, :math:`|\frac{\partial \rho_{\beta|I}}{\partial R_I}|^2` of a
particular local density.  Such potentials are useful for coarse grained
models and can be parameterized from trajectory data using `BOCS version
5 and higher. <https://github.com/noid-group/BOCS>`_

Local Density Definitions
-------------------------

Each site (atom) in a molecular dynamics simulation :math:`I` is
assigned a type :math:`t_I`.  These site types may be labeled
e.g. :math:`\alpha`, :math:`\beta`, ..., etc.  The local density,
:math:`\rho_{\beta|I}`, of :math:`\beta` particles around a central
particle, :math:`I`, of type :math:`t_I = \alpha` is defined:

.. math::

   \rho_{\beta|I} = \sum_{J\in S_{\beta|I}} \bar{w}_{\beta|\alpha}(r_{IJ})

where :math:`S_{\beta|I}` is the set of particles of type :math:`\beta`
that are not excluded from pairwise interactions involving site
:math:`I`.  :math:`\bar{w}_{\beta|\alpha}(r_{IJ})` is a normalized
indicator function, it determines the pair weighted contribution from
each particle :math:`J`, of type :math:`t_J = \beta` to the
:math:`\beta` local density around a central particle :math:`I` of type
:math:`t_I = \alpha`.

The normalized indicator function is defined by dividing a
continuous/differentiable non-negative non-increasing indicator function
:math:`w(r)` by its spatial integral :math:`[w]`.  Therefore,
:math:`\bar{w}(r)` is defined for each :math:`\beta|\alpha` pair as:

.. math::

   \bar{w}(r) = \frac{w(r)}{[w]} \theta(r - r_c) \\
   [w] = \int_{0}^{r_{c}} {4\pi r^2 w(r)} dr


where :math:`\theta` is the Heaviside function, :math:`r_c` is the local
density length-scale, and :math:`\bar{w}(r)` is defined such that
:math:`\bar{w}(r) = 0`, if :math:`r > r_c`.  The local density of
:math:`\beta` particles that surround a given particle :math:`I` is
defined by the sum of contributions from :math:`\beta` particles that
are within (<=) the cutoff :math:`r_c`.  The choice of indicator
function is specified by the user in the *indicator* argument of the
:doc:`ldd pair_coeff <pair_ldd>` command. Different options are listed
under different pages indexed on the :doc:`pair ldd <pair_ldd>` page.
Note that the weighting function and its length-scale can be defined
independently for each pair of site types.

If particle :math:`I` is of type :math:`t_I = \beta`, its contribution
to the local density :math:`\rho_{t_I|I}` can be included or excluded
from the sum without changing the net force on a pair of particles.
Practically, the only difference is a horizontal shift in the argument
of the LD potential function :math:`U_{\beta|\beta} (\rho)`.  Whether or
not to include this contribution is set by the user in the :doc:`ldd
pair_coeff <pair_ldd>` command via the *self yes* or *self no*
arguments.

For molecular systems, particles included in the local density around a
given particle :math:`I` are identified by the neighbor list for the
non-bonded interactions.  That is, a particle excluded by the neighbor
list for a central particle will also be excluded from the sum defining
its local density.

For example, consider a 5 bead chain that is simulated with default 1-2,
1-3, and 1-4 non-bonded exclusions:

:math:`A_{1}-B_{2}-C_{3}-D_{4}-E_{5}`

The local density of particles of type E surrounding site 1,
:math:`\rho_{E|1}` will include the contribution from the end chain
particle :math:`\bar{w}_{E|A}(r_{15})`.  Conversely, the local density
of particles of type B surrounding site 1, :math:`\rho_{B|1}` will only
include intermolecular contributions to the local density.  Whether or
not the local density of A particles surrounding site 1,
:math:`\rho_{A|1}`, will include the *self* intramolecular contribution
:math:`\bar{w}_{A|A}(0)` is determined by whether *self yes* or *self
no* is listed as an argument by the user in the :doc:`ldd pair_coeff
<pair_ldd>` command.

The choice of local density definition used for a particular interaction
between types, :math:`\beta|\alpha` is set by the user via the :doc:`ldd
pair_coeff <pair_ldd>` command with the different available options for
indicator functions.  Note, that the weighting function,
:math:`\bar{w}_{\beta|\alpha}`, and local density potential,
:math:`U_{\beta|\alpha}` are not symmetric with respect to
:math:`\alpha` and :math:`\beta`.  For instance, if :math:`\alpha = v`
corresponds to solvents and :math:`\beta = u` corresponds to solutes,
then the package allows the user to specify distinct weighting
functions, :math:`\bar{w}_{v|u} \neq \bar{w}_{u|v}` and distinct LD
potentials :math:`U_{v|u} \neq U_{u|v}`, for treating solvent around
solute and solute around solvent respectively.

--------------

Local Density/ Square Gradient Potentials
-----------------------------------------

The LDD pair style implements forces for two kinds of local density dependent potentials.

LD (Local Density) potentials are defined:

.. math::

   U_{LD}(\mathbf{R}) = \sum_{I} \sum_{\beta} U_{\beta|t_I}(\rho_{\beta|I}),

where the first sum is over all sites, :math:`I` , and the second sum is
over all site types, :math:`\beta`.  The form of
:math:`U_{\beta|\alpha}` is specified by the user in the :doc:`ldd
pair_coeff <pair_ldd>` command using the *potential* keyword.  For
instance, :math:`U_{\beta|\alpha}(\rho_{\beta|I})` can be defined as zero
(*noforce*) or as a *constant*, *linear*, or *quadratic* function of the
local density.  Alternatively, :math:`U_{\beta|\alpha}(\rho_{\beta|I})`
can be specified by a table of values using the *table/lin* or
*table/spline* keywords.  Additionally the *potential* arg *mdpd* is a
quadratic function in local density that can be combined with the dpd
pair_style to simulate pairwise mdpd interactions in the LD framework.
All of these forms are documented on the :doc:`pair_style ldd
<pair_ldd>` page.

SG (Square Gradient) potentials are defined:

.. math::

   U_{SG}(\mathbf{R}) = \sum_{I} \sum_{\beta} U_{\nabla; \beta|t_I} (\rho_{\beta|I}) \left | \nabla_{I} \rho_{\beta|I} \right|^2

where the first sum is over all sites, :math:`I`, the second sum is over
all site types, :math:`\beta` and the gradient :math:`\nabla_{I} =
\frac{\partial}{\partial\mathbf{R_I}}` is evaluated with respect to the
position of particle :math:`I`.  The form of the coefficient,
:math:`U_{\nabla; \beta|t_I}`, is specified by the user in the :doc:`ldd
pair_coeff <pair_ldd>` command using the (optional) *gradient* keyword.
Many of the available functional forms for the LD potential are also
available for the SG coefficients, see :doc:`pair_ldd <pair_ldd>` for
details.

------------------

For any given site :math:`I`, once the local densities
:math:`\rho_{\beta|I}` and the gradients of the local density
:math:`\nabla_{I} \rho_{\beta|I}` have been determined, the
corresponding forces are pairwise decomposable with pair forces
analogous to those stated in `Delyser 2021
<https://pubs.aip.org/aip/jcp/article/156/3/034106/2839866/Coarse-grained-models-for-local-density-gradients>`_
and `Delyser 2019
<https://pubs.aip.org/aip/jcp/article/151/22/224106/198468/Analysis-of-local-density-potentials>`_.

The per-atom local densities and gradients are owned by the *ldd* pair
style itself (recomputed every step, not stored in data or restart
files); no special atom style is required, and the data is exposed to the
rest of LAMMPS through :doc:`fix pair <fix_pair>` (see the :doc:`pair
ldd <pair_ldd>` page).

Interactions are specified per *species* pair.  The potential file
defines a set of species and the ``pair_coeff * * file.ldd S1 S2 ...``
command maps each atom type to a species (several types may map to the
same species).  Each of the :math:`N_{\text{species}}^2` *ordered*
species pairs is given on its own line.  In particular, the entry ``A
B`` specifies the interaction for atoms of species A surrounded by atoms
of species B, and is distinct from ``B A``, the interaction for species B
surrounded by species A.  This "surrounded by" convention (``A B``)
corresponds to the local density often written as ``B|A`` (species B
given a central particle of species A) in published works.

------------------

Input File Examples
-------------------

Example 1) A simple atomic input example using only tabulated LD potentials
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The input script uses a plain atomic atom style and maps the two atom
types to species A and B from the potential file:

.. code-block:: LAMMPS

   ## Init system/units
   units real
   dimension 3
   atom_style atomic

   newton on
   timestep 1

   read_data my_mix.data
   velocity all create 300.0 22345 dist gaussian

   ## pair_style ldd takes no cutoff (it comes from each indicator rc)
   pair_style ldd
   pair_coeff * * my_mix.ldd A B   # atom type 1 -> species A, type 2 -> species B

   ## Run / Output
   run_style verlet
   neigh_modify    delay 5 every 1 check yes
   neighbor        2.0 bin
   thermo 500

   fix 1 all nvt temp 300.0 300.0 100.0
   # expose the per-atom local densities (one column per species) + total energy
   fix lddout all pair 1 ldd local_density 0 total_energy 0
   dump 1 all custom 500 dump.txt id x y z f_lddout[*]
   dump_modify 1 sort id

   run 10000

with the potential file ``my_mix.ldd``:

.. code-block:: LAMMPS

   # all four ordered species pairs must be listed (use "ignore" for none)
   A A indicator dpd 0.0 6.5 self yes potential table/lin LD_table.1.1.dat
   A B ignore
   B A indicator dpd 0.0 5.5 self no  potential table/lin LD_table.2.1.dat
   B B indicator dpd 0.0 5.5 self yes potential table/lin LD_table.2.2.dat

Note in the above example, only species pairs A|A, A|B, and B|B interact
according to LD potentials (the A B entry, corresponding to B|A, is ignored).

Example 2) Molecular input example that layers tabulated pair, LD and SG interactions via :doc:`pair_style hybrid <pair_hybrid>`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The molecular topology uses the ordinary *full* atom style; the *ldd*
interactions are layered in with :doc:`pair_style hybrid/overlay
<pair_hybrid>`:

.. code-block:: LAMMPS

   ## Initialize system /units

   units real
   dimension 3
   atom_style full

   newton on
   timestep 1

   bond_style harmonic

   region my_box block 0 40 0 40 0 40
   create_box 2 my_box bond/types 1 extra/bond/per/atom 1

   mass 1 59.044
   mass 2 59.044

   pair_style hybrid/overlay ldd table spline 3000

   molecule DI DIMER.mol
   create_atoms 0 random 250 885401 my_box mol DI 454353

   include ff_bonded_params.lammps
   velocity all create 300.0 22345 dist gaussian

   ## Pair styles
   pair_coeff * * table lammps_nb_ALL.table nb_All 15.0
   pair_coeff * * ldd dimer.ldd A B   # type 1 -> species A, type 2 -> species B

   ## Run / Output
   run_style verlet
   neigh_modify    delay 5 every 1 check yes
   neighbor        2.0 bin
   thermo 500

   fix 1 all nvt temp 300.0 300.0 100.0
   run 10000

with the potential file ``dimer.ldd``:

.. code-block:: LAMMPS

   A A indicator lucy 0.0 6.5 self no potential table/spline LD.1.1.dat gradient table/gradspline SG.1.1.dat
   A B indicator lucy 0.0 6.5 self no potential table/spline LD.1.2.dat
   B A indicator lucy 0.0 6.5 self no potential table/spline LD.2.1.dat
   B B indicator lucy 0.0 6.5 self no potential table/spline LD.2.2.dat gradient table/gradspline SG.2.2.dat

Note in the above example LD cross interactions are both specified, but not necessarily symmetric.
Also in this example, the cross interactions opt not to use the optional *gradient* keyword.
Conversely the like-species interactions (A|A, B|B) layer SG potentials on top of the already defined LD and pair interactions.

See ``examples/PACKAGES/bocs`` for full input files for a 4-site chain example (topology like A-B-C-D).

See :doc:`pair_style ldd <pair_ldd>` for all pair_coeff *args* options that exist, all restrictions, and more examples.

-------------------------

read_data file format
------------------------------

The *ldd* pair style does not require a special atom style: the local
densities, their gradients, and the associated energies are owned by the
pair style and recomputed every step, so they are never read from or
written to a data file.  The ``read_data`` input therefore follows the
usual format for whatever atom style you choose, e.g. *atomic* for a
simple system or *full* for a molecular one.

-------------------------

Outputting the local-density data
---------------------------------

For each central particle :math:`I` and each species :math:`\beta`, there
is a local density of :math:`\beta` particles that surround :math:`I`,
:math:`\rho_{\beta|I}`, and a corresponding gradient :math:`\frac{\partial
\rho_{\beta|I}}{\partial \boldsymbol{R}_I}`.  These quantities (along with
the per-species LD/SG energies and the per-atom total energy) are made
available to the rest of LAMMPS through the :doc:`fix pair <fix_pair>`
command, which copies the chosen pair-style fields into a per-atom array.
:doc:`\pair_style ldd <pair_ldd>` describes the available options and quantities output.

For example:

.. code-block:: LAMMPS

   fix lddout all pair 1 ldd local_density 0 grad_density 0 total_energy 0
   dump 1 all custom 500 dump.txt id x y z vx vy vz fx fy fz f_lddout[*]
   dump_modify 1 sort id

writes a custom trajectory that includes the local densities and their
gradients alongside the usual positions, velocities, and forces (one
column per species, three per species for the gradient).  Such a
trajectory can be used with the `Bottom-up Open-source Coarse-graining
Software <https://github.com/noid-group/BOCS>`_, which can parameterize
LD/SG potentials from atomistic data and convert LAMMPS trajectories to
.trr files for analysis with `gromacs <https://www.gromacs.org/>`_ tools.

-------------

When to use LD/SG potentials:
-----------------------------

This package has been primarily developed in the context of
parameterizing/simulating bottom up coarse grained models.  We have
found that LD/SG potentials are particularly useful for supplementing
pair potentials when simulating CG models in the NPT ensemble and in
different interfacial environments.  In particular, such potentials
dramatically improve CG descriptions of the pressure-density equation of
state, interfacial profiles, liquid-vapor coexistence, and
transferability between bulk and interfacial systems.  SG potentials are
slower to simulate than LD potentials alone, so where applicable for
coarse graining we recommend trying to construct/simulate LD potentials
first and then [if necessary] adding SG potentials to refine further.
CG force fields with LD and SG potentials can be parameterized via the
multi-scale coarse graining method implemented in `BOCSv5 and higher
<https://github.com/noid-group/BOCS>`_.

-------------

A number of works have been published using this package (with `BOCS
<https://github.com/noid-group/BOCS>`_) in its development stage for
different CG studies and are listed below:

.. _DeLyser1:

**(DeLyser)** DeLyser, Noid. "Extending Pressure-Matching to Inhomogeneous Systems via Local-Density Potentials." The Journal of Chemical Physics, 147, no. 13: 134111 (2017)

.. _DeLyser2:

**(DeLyser)** DeLyser, Noid. "Analysis of local density potentials." The Journal of Chemical Physics, 151, no. 22:224106 (2019)

.. _DeLyser3:

**(DeLyser)** DeLyser, Noid "Bottom-up coarse-grained models for external fields and interfaces" The Journal of Chemical Physics 153, 224103 (2020)

.. _DeLyser4:

**(DeLyser)** DeLyser, Noid "Coarse-grained models for local density gradients" The Journal of Chemical Physics 156, 034106 (2021)

.. _Szukalo:

**(Szukalo)** Szukalo, Noid "A temperature-dependent length-scale for transferable local density potentials". The Journal of Chemical Physics, 159, 074104 (2023)

.. _Lesniewski1:

**(Lesniewski)** Lesniewski, Remsing, Noid "Revisiting what we lose by coarse-graining: Modeling cooperative hydrophobic phenomena with short-ranged, pair-additive forces." The Journal of Chemical Physics, 164, 204111 (2026)

.. _Dutta:

**(Dutta)** Dutta, Lesniewski, Qaisrani, Noid, Andrienko, Nikoubashman. "Accurate coarse-graining of conjugated organic molecules in melts and thin films using density-dependent potentials". The Journal of Chemical Theory and Computation, 22, 7, 3697-3708 (2026)

.. _Lesniewski2:

**(Lesniewski)** Lesniewski, DeLyser, W. G. Noid "Progress toward a better BOCS: systematic coarse-graining with local density potentials" In Prep 2026
