Gibbs ensemble Monte Carlo
==========================

.. versionadded:: TBD

This Howto explains how to compute the densities of coexisting vapor
and liquid phases with Gibbs ensemble Monte Carlo (GEMC) using the
:doc:`fix gemc <fix_gemc>` command.  As an example it uses the
Lennard-Jones (LJ) fluid truncated and shifted at a cutoff of
:math:`2.5\sigma` (LJTS), for which very accurate reference data is
available.  The complete input script is ``examples/mc/in.gemc.lj``.

Background
----------

A direct way to simulate vapor-liquid coexistence is to put a slab of
liquid in contact with its vapor in one elongated simulation box.  This
requires a large system, because the interfaces between the phases
disturb the densities near them.  GEMC avoids interfaces altogether
:ref:`(Panagiotopoulos) <Panagiotopoulos2>`.  It uses two separate,
periodic simulation boxes, one for each phase.  Three kinds of Monte
Carlo (MC) moves bring the two boxes into equilibrium with each other:

* atoms are displaced within each box (thermal equilibrium of each box),
* volume is moved from one box to the other while the total volume is
  kept constant (equal pressure in both boxes), and
* atoms are moved from one box to the other while the total number of
  atoms is kept constant (equal chemical potential in both boxes).

If the total density of the two boxes lies inside the two-phase region
at the chosen temperature, the system separates by itself: one box
ends up with the vapor density and the other with the liquid density.
The textbook by Frenkel and Smit :ref:`(Frenkel) <Frenkel4>` gives a
detailed derivation of the method and of the acceptance rules used by
fix gemc.

Running LAMMPS with two boxes
-----------------------------

In LAMMPS, each box is a separate :doc:`partition <Run_options>`.  The
number of partitions and the processors used for each are set with the
``-partition`` (or ``-p``) command-line switch.  The following commands
run the example on 2, 3, and 4 MPI processes, respectively:

.. code-block:: bash

   mpirun -np 2 lmp -in in.gemc.lj -partition 1 1
   mpirun -np 3 lmp -in in.gemc.lj -partition 1 2
   mpirun -np 4 lmp -in in.gemc.lj -partition 1 3

Since the liquid box usually contains many more atoms than the vapor
box, it is often best to give it more processors, e.g. the second box in
``-partition 1 3``.  Each partition writes its own log file
(``log.lammps.0`` and ``log.lammps.1``) and screen output
(``screen.0`` and ``screen.1``).  The screen and log file of the
whole run (``log.lammps``) receive a progress report of fix gemc.

Both partitions execute the same input script.  Everything that should
differ between the boxes can be set with a :doc:`world-style variable
<variable>`, which has one value per partition, or with the
:doc:`partition <partition>` command.

Setting up the example
----------------------

The example starts the first box at a low density and the second one
at a high density, each with a number of atoms roughly proportional to
the expected amount of each phase:

.. code-block:: LAMMPS

   variable        T index 0.9                  # temperature
   variable        rho world 0.05 0.65          # initial densities of the boxes
   variable        N world 100 400              # initial numbers of atoms
   variable        L equal (v_N/v_rho)^(1.0/3.0)

   units           lj
   atom_style      atomic
   region          box block 0 ${L} 0 ${L} 0 ${L}
   create_box      1 box
   create_atoms    1 random ${N} 4321 NULL overlap 0.9 maxtry 1000
   mass            1 1.0

   pair_style      lj/cut 2.5
   pair_modify     shift yes
   pair_coeff      1 1 1.0 1.0

   minimize        1.0e-4 1.0e-6 100 1000
   reset_timestep  0

The starting densities do not need to be accurate, but it helps to start
near the expected coexistence densities.  Random placement can create
atoms that overlap strongly, so a short :doc:`energy minimization
<minimize>` removes the overlaps before the Monte Carlo run.  Note that
the energy is shifted to zero at the cutoff with :doc:`pair_modify
shift yes <pair_modify>`.  Unlike in MD, this shift changes the results
of MC moves that change the number of pairs, so the model definition
(and the reference data used for comparison) must include it.

The Monte Carlo moves are defined by the fix gemc command:

.. code-block:: LAMMPS

   fix             mc all gemc 1 100 100 2 ${T} 0.3 0.05 29494
   compute_modify  thermo_temp dynamic/dof yes

Every timestep (*N* = 1), the fix attempts on average 100
translations, 100 atom exchanges, and 2 volume changes with a maximum
displacement of :math:`0.3\sigma` and a maximum change of the logarithm
of the volume ratio of the boxes of 0.05.  No time integration fix is
defined, so this is a pure MC simulation, and the "timestep" merely
counts MC cycles.  LAMMPS warns that atoms will not move, which can be
ignored.  Because atoms are created and deleted, :doc:`compute_modify
dynamic/dof <compute_modify>` makes the temperature output use the
current number of atoms.

Output and monitoring
---------------------

The most important quantities are the number density of each box and
the acceptance ratios of the MC moves:

.. code-block:: LAMMPS

   variable        rho_n equal atoms/vol
   thermo_style    custom step atoms vol pe v_rho_n f_mc[1] f_mc[2] &
                   f_mc[3] f_mc[4] f_mc[5] f_mc[6]
   thermo_modify   norm no
   thermo          1000
   run             20000

The columns *f_mc[1]* to *f_mc[6]* are the cumulative numbers of
attempted and successful translations, exchanges, and volume changes
of the box.  The ratio of successes to attempts is the acceptance
ratio of each kind of move:

* Translations: adjust *displace* so that 30% to 50% of the moves in the
  liquid box are accepted.  In the vapor box almost all translations are
  accepted, which is normal.
* Volume changes: adjust *maxdlogvolratio* so that 30% to 50% are
  accepted.
* Both can be adjusted automatically during equilibration with the
  *tune* keyword of :doc:`fix gemc <fix_gemc>`; its page shows how to
  continue with the adjusted values in the production run.
* Exchanges: the acceptance ratio cannot be tuned, it gets smaller for
  denser liquids.  For the LJTS fluid at :math:`T^* = 0.7` only about
  0.1% of the exchanges are accepted, while at :math:`T^* = 1.0` it is
  about 3.5%.  The exchanges are the slowest process in reaching
  equilibrium, so *X* must be large: every atom should be exchanged
  many times during the run.  With too few exchanges, the vapor density
  drifts slowly during the whole run and the averages are biased.

The number of atoms and the density of both boxes will drift at the
beginning of the run and then fluctuate around constant values.  Only
the data after this equilibration period should be analyzed.

Analyzing the results
---------------------

The coexistence densities are the average densities of the two boxes
after equilibration.  Because consecutive values are correlated, the
statistical error is best estimated by block averaging: divide the
production data into e.g. 5 to 10 blocks, average each block, and use
the standard deviation of the block averages divided by the square root
of the number of blocks.  Additional checks are that the pressures of
the two boxes agree within their (large) fluctuations and that the
average number of atoms in each box stays well above zero.

Near the critical temperature the densities of the two phases become
similar, the fluctuations become large, and the boxes may swap their
identity during a run: the box that contained the liquid then contains
the vapor and vice versa.  In that case average the density of the
denser box and of the less dense box at each output step instead of
the density of each box.  Close to the critical point a GEMC simulation
can no longer separate the phases; for LJ this happens a few percent
below the critical temperature for systems of a few hundred atoms.

Example results
---------------

The following table compares the densities from runs of 100,000 MC
cycles (last 80,000 cycles analyzed) with 500 atoms in total, using the
settings of the example input with adjusted *displace* and
*maxdlogvolratio*, to the coexistence densities from the equation of state of Thol et al.
:ref:`(Thol) <Thol1>` for the LJTS fluid.  The numbers in
parentheses are the statistical errors in units of the last digit.

.. list-table::
   :header-rows: 1
   :widths: 10 20 20 20 20

   * - :math:`T^*`
     - :math:`\rho^*_v` (GEMC)
     - :math:`\rho^*_v` (EOS)
     - :math:`\rho^*_l` (GEMC)
     - :math:`\rho^*_l` (EOS)
   * - 0.7
     - 0.00740(17)
     - 0.00746
     - 0.7883(7)
     - 0.7869
   * - 0.8
     - 0.0195(3)
     - 0.0200
     - 0.7322(7)
     - 0.7311
   * - 0.9
     - 0.0454(9)
     - 0.0453
     - 0.6649(13)
     - 0.6643
   * - 1.0
     - 0.099(3)
     - 0.0983
     - 0.573(3)
     - 0.5733

Performance
-----------

For pair styles without many-body terms, fix gemc computes the energy
change of a translation or exchange from the interactions of the moved
atom only, which is much faster than recomputing the total energy of
the box: the cost of such a move does not depend on the size of the
box.  See the :doc:`fix gemc <fix_gemc>` page for the conditions.
Volume changes always require the total energy, but they are also
needed much less often.  Additional tips:

* Use *N* = 1 and a number of moves per invocation that is comparable to
  or larger than the number of atoms.  Each invocation has an overhead
  of a neighbor list build and a force computation.
* A larger communication cutoff, set with :doc:`comm_modify cutoff
  <comm_modify>`, allows larger displacements with single-atom
  energies and makes it less frequent that the ghost atoms must be
  rebuilt during the moves.
* Give the liquid box more MPI processes than the vapor box.  Very small
  systems (a few hundred atoms) run fastest on one process per box.

Molecular fluids
----------------

With the *mol* keyword of :doc:`fix gemc <fix_gemc>`, whole molecules
are moved and exchanged.  The example ``examples/mc/in.gemc.co2``
computes the vapor-liquid coexistence of carbon dioxide with the TraPPE
model :ref:`(Potoff) <Potoff1>`: rigid linear molecules with partial
charges, computed with a long-range solver.  The molecules are created
from a molecule template, which is also given to the fix:

.. code-block:: LAMMPS

   molecule        co2mol CO2.txt
   create_atoms    0 random ${N} 7771 NULL mol co2mol 464563 overlap 2.5 maxtry 1000
   ...
   fix             mc all gemc 1 50 100 2 ${T} 0.5 0.1 7771 &
                   mol co2mol maxangle 30 tune 20 0.4 0.4

The translation moves are now translations and rotations of whole
molecules, and volume exchanges move the molecule centers, so the
molecules stay rigid in a pure MC simulation even though the bond and
angle force constants are zero.  Flexible molecules can be simulated
by combining the fix with time integration, which samples the internal
degrees of freedom; an exchanged molecule keeps its conformation and
velocities.  Since every move of a molecule requires the total energy
of the box, molecular simulations are much more expensive than atomic
ones.  Give the liquid box most of the processors.

The *tune* keyword adjusts the maximum displacement, rotation angle,
and volume change during the equilibration run; the example then
re-defines the fix with the adjusted values for the production run.

The following table compares the result of a run with 256 molecules
(3000 equilibration and 12000 production cycles) at 250 K to the
experimental saturated densities :ref:`(NIST) <NIST1>`, which the TraPPE
model reproduces closely:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20

   * -
     - vapor (g/cm\ :sup:`3`)
     - liquid (g/cm\ :sup:`3`)
   * - GEMC, TraPPE
     - 0.0461(12)
     - 1.041(7)
   * - Experiment
     - 0.0466
     - 1.046

Restrictions
------------

Fix gemc supports single atoms and, with the *mol* keyword, molecules
of a single kind.  Atoms of different types can be exchanged, so the
coexistence of simple mixtures of atoms can be computed as well.
Orthogonal and triclinic periodic boxes in 3d, charges, and long-range
solvers are supported.  Constraints (fix rigid, fix shake) are not
supported.

----------

.. _Panagiotopoulos2:

**(Panagiotopoulos)** Panagiotopoulos, Mol Phys, 61, 813-826 (1987).

.. _Frenkel4:

**(Frenkel)** Frenkel and Smit, Understanding Molecular Simulation,
3rd edition, Academic Press, London, 2023.

.. _Potoff1:

**(Potoff)** Potoff and Siepmann, AIChE J, 47, 1676-1682 (2001).

.. _NIST1:

**(NIST)** Lemmon, McLinden, and Friend, "Thermophysical Properties of
Fluid Systems" in NIST Chemistry WebBook, NIST Standard Reference
Database Number 69, https://webbook.nist.gov

.. _Thol1:

**(Thol)** Thol, Rutkai, Span, Vrabec, Lustig, Int J Thermophys, 36,
25-43 (2015).
