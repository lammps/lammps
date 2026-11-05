Benchmarks
==========

..
   FIXME:
   Current LAMMPS performance is discussed on the `Benchmarks page <https://www.lammps.org/bench.html>`_ of the `LAMMPS website <lws_>`_
   where timings and parallel efficiency are listed.  The page has
   several sections, which are briefly described below:

   * CPU performance on 5 standard problems, strong and weak scaling
   * GPU performance on same and related problems
   * Comparison of cost of interatomic potentials
   * Performance of huge, billion-atom problems

   The 5 standard problems are as follow:

   #. LJ = atomic fluid, Lennard-Jones potential with 2.5 sigma cutoff (55
      neighbors per atom), NVE integration
   #. Chain = bead-spring polymer melt of 100-mer chains, FENE bonds and LJ
      pairwise interactions with a :math:`2^{\frac{1}{6}}` sigma cutoff (5 neighbors per
      atom), NVE integration
   #. EAM = metallic solid, Cu EAM potential with 4.95 Angstrom cutoff (45
      neighbors per atom), NVE integration
   #. Chute = granular chute flow, frictional history potential with 1.1
      sigma cutoff (7 neighbors per atom), NVE integration
   #. Rhodo = rhodopsin protein in solvated lipid bilayer, CHARMM force
      field with a 10 Angstrom LJ cutoff (440 neighbors per atom),
      particle-particle particle-mesh (PPPM) for long-range Coulombics, NPT
      integration

   Input files for these 5 problems are provided in the bench directory
   of the LAMMPS distribution.  Each has 32,000 atoms and runs for 100
   timesteps.  The size of the problem (number of atoms) can be varied
   using command-line switches as described in the ``bench/README`` file.
   This is an easy way to test performance and either strong or weak
   scalability on your machine.

   The bench directory includes a few ``log.*`` files that show performance
   of these 5 problems on 1 or 4 cores of Linux desktop.  The ``bench/FERMI``
   and ``bench/KEPLER`` directories have input files and scripts and instructions
   for running the same (or similar) problems using OpenMP or GPU or Xeon
   Phi acceleration options.  See the ``README`` files in those directories and the
   :doc:`Accelerator packages <Speed_packages>` pages for instructions on how
   to build LAMMPS and run on that kind of hardware.

   The ``bench/POTENTIALS`` directory has input files which correspond to the
   table of results on the
   `Potentials <https://www.lammps.org/bench.html#potentials>`_ section of
   the Benchmarks web page.  So you can also run those test problems on
   your machine.

   The `billion-atom <https://www.lammps.org/bench.html#billion>`_ section
   of the Benchmarks web page has performance data for very large
   benchmark runs of simple Lennard-Jones (LJ) models, which use the
   ``bench/in.lj`` input script.

For all the benchmarks, a useful metric is the CPU cost per atom per
timestep.  Since performance scales roughly linearly with problem size
and timesteps for all LAMMPS models (i.e. interatomic or coarse-grained
potentials), the run time of any problem using the same model (atom
style, force field, cutoff, etc) can then be estimated.

Performance on a parallel machine can also be predicted from one-core
or one-node timings if the parallel efficiency can be estimated.  The
communication bandwidth and latency of a particular parallel machine
affects the efficiency.  On most machines LAMMPS will give a parallel
efficiency on these benchmarks above 50% so long as the number of
atoms/core is a few 100 or greater, and closer to 100% for large
numbers of atoms/core.  This is for all-MPI mode with one MPI task per
core.  For nodes with accelerator options or hardware (OpenMP, GPU,
Phi), you should first measure single node performance.  Then you can
estimate parallel performance for multi-node runs using the same logic
as for all-MPI mode, except that now you will typically need many more
atoms/node to achieve good scalability.

Floating-point operations per atom and timestep
"""""""""""""""""""""""""""""""""""""""""""""""

The performance of CPUs and High-Performance Computing (HPC) clusters
for scientific applications is often given as the number of
floating-point operations per second (FLOP/s); see for example the
`TOP500 list of the most powerful supercomputers
<https://www.top500.org/>`_.  For that list, the *sustained* FLOP/s
value of the entire supercomputer is determined from how fast a
benchmark application with a known number of floating-point operations
executes.  This is then compared to the theoretically possible (*peak*)
FLOP/s value of the machine.  The benchmark for the TOP500 list (HPL)
solves a large dense system of linear equations, for which highly
optimized code can reach a large fraction of the peak: on the June 2026
list, GPU-accelerated systems reached about 50 to 80 percent and
CPU-only systems about 60 to 90 percent of their peak.  Historic vector
computers reached similar or even higher fractions, up to 93 percent.
Some examples:

.. table::

   +--------------------------+----------------------------------------+-----------+
   | System                   | Processors                             |   Percent |
   |                          |                                        |   of peak |
   +==========================+========================================+===========+
   | *GPU-accelerated systems (June 2026 list)*                                    |
   +--------------------------+----------------------------------------+-----------+
   | El Capitan               | AMD MI300A                             |        64 |
   +--------------------------+----------------------------------------+-----------+
   | Frontier                 | AMD MI250X                             |        66 |
   +--------------------------+----------------------------------------+-----------+
   | Aurora                   | Intel GPU Max                          |        51 |
   +--------------------------+----------------------------------------+-----------+
   | Eagle                    | NVIDIA H100                            |        66 |
   +--------------------------+----------------------------------------+-----------+
   | JUPITER Booster          | NVIDIA GH200                           |        82 |
   +--------------------------+----------------------------------------+-----------+
   | *CPU-only systems (June 2026 list)*                                           |
   +--------------------------+----------------------------------------+-----------+
   | LineShine                | LX2, 304 cores                         |        80 |
   +--------------------------+----------------------------------------+-----------+
   | Fugaku                   | Fujitsu A64FX                          |        82 |
   +--------------------------+----------------------------------------+-----------+
   | Shaheen III CPU          | AMD EPYC 9654 (Zen 4)                  |        90 |
   +--------------------------+----------------------------------------+-----------+
   | Frontera                 | Intel Xeon 8280 (Cascade Lake)         |        61 |
   +--------------------------+----------------------------------------+-----------+
   | *Vector computers (historic)*                                                 |
   +--------------------------+----------------------------------------+-----------+
   | Cray Y-MP C916           | Cray vector processors                 |        90 |
   +--------------------------+----------------------------------------+-----------+
   | S-3800/480               | Hitachi vector processors              |        89 |
   +--------------------------+----------------------------------------+-----------+
   | Numerical Wind Tunnel    | Fujitsu vector processors              |        81 |
   +--------------------------+----------------------------------------+-----------+
   | Cray X1E                 | Cray vector processors                 |        85 |
   +--------------------------+----------------------------------------+-----------+
   | Earth Simulator          | NEC SX-9/E                             |        93 |
   +--------------------------+----------------------------------------+-----------+

On modern CPUs, the floating-point operations can be counted with the
hardware performance counters of the CPU, e.g. with the Linux ``perf``
tool.  For the numbers below, double precision floating-point
operations were counted, and a fused multiply-add (FMA) counts as two.
Dividing the measured FLOP/s value by the atom-steps per second from
the "Performance:" line of the LAMMPS output gives the average number of
floating-point operations (FLOP) that LAMMPS executes per atom and
timestep.  This is a measure for the computational cost of a model,
which is mostly determined by the model itself and how it is
implemented in LAMMPS: for the same input, it differed by less than 1
percent between different CPUs and between the GCC compiler
optimization settings ``-O2`` and ``-O3 -march=x86-64-v3``, while the
Intel compiler changed it by -5 to +10 percent.

The table below lists this number for most of the inputs in the
``bench/POTENTIALS`` directory, which were used to compare the cost of
interatomic potentials in :ref:`(Thompson) <Thompson2022>`.  The numbers
are averages over the 100 MD steps of each input without the setup,
measured on a single CPU core; they include the neighbor list builds
and the time integration.  They also depend on the input: the cutoff
and density (i.e. the number of neighbors per atom), the neighbor list
settings, and for some models the state of the system.  Nevertheless,
they indicate how much work a specific model requires.  The log files in
the ``bench/POTENTIALS`` directory show the timings of these inputs on
one core of three different CPUs.

.. table::

   +---------------------------------------------+--------------+
   | Model (pair style)                          | FLOP per     |
   |                                             | atom-step    |
   +=============================================+==============+
   | gran/hooke/history                          |          170 |
   +---------------------------------------------+--------------+
   | lj/cut with FENE bonds, 1.12 sigma cutoff   |          200 |
   +---------------------------------------------+--------------+
   | dpd                                         |          890 |
   +---------------------------------------------+--------------+
   | lj/cut, 2.5 sigma cutoff                    |          950 |
   +---------------------------------------------+--------------+
   | sw                                          |        1,700 |
   +---------------------------------------------+--------------+
   | eam                                         |        1,800 |
   +---------------------------------------------+--------------+
   | tersoff                                     |        3,100 |
   +---------------------------------------------+--------------+
   | rebo                                        |        3,500 |
   +---------------------------------------------+--------------+
   | lj/cut/coul/long with PPPM (SPC/E water)    |       11,000 |
   +---------------------------------------------+--------------+
   | meam                                        |       15,000 |
   +---------------------------------------------+--------------+
   | lj/charmm/coul/long with PPPM (CHARMM)      |       16,000 |
   +---------------------------------------------+--------------+
   | airebo                                      |       33,000 |
   +---------------------------------------------+--------------+
   | gayberne                                    |       56,000 |
   +---------------------------------------------+--------------+
   | reaxff with qeq/reaxff                      |      160,000 |
   +---------------------------------------------+--------------+
   | comb with qeq/comb                          |      210,000 |
   +---------------------------------------------+--------------+
   | eff/cut                                     |      260,000 |
   +---------------------------------------------+--------------+
   | snap (2J=8) with zbl                        |    1,100,000 |
   +---------------------------------------------+--------------+

Comparing the measured FLOP/s with the theoretical peak of a CPU core
(the number of double precision FLOP per clock cycle that its vector
units can execute, times the clock) shows how well the hardware is
used.  For these inputs, one core executed 0.3 to 1.3 FLOP per clock
cycle, i.e. about 2 to 8 percent of its peak on an AMD Zen 4 and an
Intel Broadwell CPU and 1 to 3 percent on an Intel Cascade Lake CPU,
whose AVX-512 units double the peak FLOP/s value.  This was measured
with LAMMPS compiled by GCC with ``-O2`` or with
``-O3 -march=x86-64-v3``.  SNAP reached up to 12 percent.

These numbers are typical for molecular dynamics and most other
scientific software: only hand-tuned code for dense linear algebra
comes close to the peak.  For comparison, the HPL benchmark reached 90
percent of the peak on the Zen 4 CPUs of Shaheen III and 61 percent on
the Cascade Lake CPUs of Frontera (see the table above).  Most
of the floating-point work in LAMMPS is scalar and depends on indirect
memory accesses through the neighbor lists.  CPUs with wider vector
units have a higher peak, but styles without explicit vectorization do
not get faster with it.  The vectorized styles of the :doc:`INTEL
package <Speed_intel>` are an exception: on the Cascade Lake CPU, with
LAMMPS compiled by the Intel oneAPI compilers, the INTEL package
versions of the lj/cut, sw, tersoff, lj/charmm/coul/long, and snap
pair styles executed 75 to 99 percent of their floating-point work
with 512-bit vector instructions and ran 1.6 to 3.3 times faster than
the corresponding plain styles.  Measured as useful work, they still
reached only about 4 to 7 percent of the peak (14 percent for SNAP),
since the vectorized code also computes masked or padded vector
elements (15 to 90 percent more floating-point operations) and the CPU
lowers its clock for AVX-512 instructions.

----------

.. _Thompson2022:

**(Thompson)** A. P. Thompson, H. M. Aktulga, R. Berger,
D. S. Bolintineanu, W. M. Brown, P. S. Crozier, P. J. in 't Veld,
A. Kohlmeyer, S. G. Moore, T. D. Nguyen, R. Shan, M. J. Stevens,
J. Tranchida, C. R. Trott, S. J. Plimpton, Comp Phys Comm, 271, 108171
(2022)

.. _lws: https://www.lammps.org
