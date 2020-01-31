Benchmarks
==========

Current LAMMPS performance is discussed on the `Benchmarks page <http://lammps.sandia.gov/bench.html>`_ of the `LAMMPS website <lws_>`_
where timings and parallel efficiency are listed.  The page has
several sections, which are briefly described below:

* CPU performance on 5 standard problems, strong and weak scaling
* GPU and Xeon Phi performance on same and related problems
* Comparison of cost of interatomic potentials
* Performance of huge, billion-atom problems

The 5 standard problems are as follow:

#. LJ = atomic fluid, Lennard-Jones potential with 2.5 sigma cutoff (55
   neighbors per atom), NVE integration
#. Chain = bead-spring polymer melt of 100-mer chains, FENE bonds and LJ
   pairwise interactions with a 2\^(1/6) sigma cutoff (5 neighbors per
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
using command-line switches as described in the bench/README file.
This is an easy way to test performance and either strong or weak
scalability on your machine.

The bench directory includes a few log.\* files that show performance
of these 5 problems on 1 or 4 cores of Linux desktop.  The bench/FERMI
and bench/KEPLER dirs have input files and scripts and instructions
for running the same (or similar) problems using OpenMP or GPU or Xeon
Phi acceleration options.  See the README files in those dirs and the
:doc:`Speed packages <Speed_packages>` doc pages for instructions on how
to build LAMMPS and run on that kind of hardware.

The bench/POTENTIALS directory has input files which correspond to the
table of results on the
`Potentials <http://lammps.sandia.gov/bench.html#potentials>`_ section of
the Benchmarks web page.  So you can also run those test problems on
your machine.

The `billion-atom <http://lammps.sandia.gov/bench.html#billion>`_ section
of the Benchmarks web page has performance data for very large
benchmark runs of simple Lennard-Jones (LJ) models, which use the
bench/in.lj input script.


----------


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
