USER-OMP package
================

The USER-OMP package was developed by Axel Kohlmeyer at Temple
University.  It provides optimized and multi-threaded versions
of many pair styles, nearly all bonded styles (bond, angle, dihedral,
improper), several Kspace styles, and a few fix styles.  It uses
the OpenMP interface for multi-threading, but can also be compiled
without OpenMP support, providing optimized serial styles in that case.

**Required hardware/software:**

To enable multi-threading, your compiler must support the OpenMP interface.
You should have one or more multi-core CPUs, as multiple threads can only be
launched by each MPI task on the local node (using shared memory).

**Building LAMMPS with the USER-OMP package:**

See the :ref:`Build extras <user-omp>` doc page for
instructions.

**Run with the USER-OMP package from the command line:**

These examples assume one or more 16-core nodes.


.. parsed-literal::

   env OMP_NUM_THREADS=16 lmp_omp -sf omp -in in.script           # 1 MPI task, 16 threads according to OMP_NUM_THREADS
   lmp_mpi -sf omp -in in.script                                  # 1 MPI task, no threads, optimized kernels
   mpirun -np 4 lmp_omp -sf omp -pk omp 4 -in in.script           # 4 MPI tasks, 4 threads/task
   mpirun -np 32 -ppn 4 lmp_omp -sf omp -pk omp 4 -in in.script   # 8 nodes, 4 MPI tasks/node, 4 threads/task

The mpirun or mpiexec command sets the total number of MPI tasks used
by LAMMPS (one or multiple per compute node) and the number of MPI
tasks used per node.  E.g. the mpirun command in MPICH does this via
its -np and -ppn switches.  Ditto for OpenMPI via -np and -npernode.

You need to choose how many OpenMP threads per MPI task will be used
by the USER-OMP package.  Note that the product of MPI tasks \*
threads/task should not exceed the physical number of cores (on a
node), otherwise performance will suffer.

As in the lines above, use the "-sf omp" :doc:`command-line switch <Run_options>`, which will automatically append "omp" to
styles that support it.  The "-sf omp" switch also issues a default
:doc:`package omp 0 <package>` command, which will set the number of
threads per MPI task via the OMP\_NUM\_THREADS environment variable.

You can also use the "-pk omp Nt" :doc:`command-line switch <Run_options>`, to explicitly set Nt = # of OpenMP threads
per MPI task to use, as well as additional options.  Its syntax is the
same as the :doc:`package omp <package>` command whose doc page gives
details, including the default values used if it is not specified.  It
also gives more details on how to set the number of threads via the
OMP\_NUM\_THREADS environment variable.

**Or run with the USER-OMP package by editing an input script:**

The discussion above for the mpirun/mpiexec command, MPI tasks/node,
and threads/MPI task is the same.

Use the :doc:`suffix omp <suffix>` command, or you can explicitly add an
"omp" suffix to individual styles in your input script, e.g.


.. parsed-literal::

   pair_style lj/cut/omp 2.5

You must also use the :doc:`package omp <package>` command to enable the
USER-OMP package.  When you do this you also specify how many threads
per MPI task to use.  The command doc page explains other options and
how to set the number of threads via the OMP\_NUM\_THREADS environment
variable.

**Speed-ups to expect:**

Depending on which styles are accelerated, you should look for a
reduction in the "Pair time", "Bond time", "KSpace time", and "Loop
time" values printed at the end of a run.

You may see a small performance advantage (5 to 20%) when running a
USER-OMP style (in serial or parallel) with a single thread per MPI
task, versus running standard LAMMPS with its standard un-accelerated
styles (in serial or all-MPI parallelization with 1 task/core).  This
is because many of the USER-OMP styles contain similar optimizations
to those used in the OPT package, described in :doc:`Section 5.3.5 <Speed_opt>`.

With multiple threads/task, the optimal choice of number of MPI
tasks/node and OpenMP threads/task can vary a lot and should always be
tested via benchmark runs for a specific simulation running on a
specific machine, paying attention to guidelines discussed in the next
sub-section.

A description of the multi-threading strategy used in the USER-OMP
package and some performance examples are `presented here <http://sites.google.com/site/akohlmey/software/lammps-icms/lammps-icms-tms2011-talk.pdf?attredirects=0&d=1>`_

**Guidelines for best performance:**

For many problems on current generation CPUs, running the USER-OMP
package with a single thread/task is faster than running with multiple
threads/task.  This is because the MPI parallelization in LAMMPS is
often more efficient than multi-threading as implemented in the
USER-OMP package.  The parallel efficiency (in a threaded sense) also
varies for different USER-OMP styles.

Using multiple threads/task can be more effective under the following
circumstances:

* Individual compute nodes have a significant number of CPU cores but
  the CPU itself has limited memory bandwidth, e.g. for Intel Xeon 53xx
  (Clovertown) and 54xx (Harpertown) quad-core processors.  Running one
  MPI task per CPU core will result in significant performance
  degradation, so that running with 4 or even only 2 MPI tasks per node
  is faster.  Running in hybrid MPI+OpenMP mode will reduce the
  inter-node communication bandwidth contention in the same way, but
  offers an additional speedup by utilizing the otherwise idle CPU
  cores.
* The interconnect used for MPI communication does not provide
  sufficient bandwidth for a large number of MPI tasks per node.  For
  example, this applies to running over gigabit ethernet or on Cray XT4
  or XT5 series supercomputers.  As in the aforementioned case, this
  effect worsens when using an increasing number of nodes.
* The system has a spatially inhomogeneous particle density which does
  not map well to the :doc:`domain decomposition scheme <processors>` or
  :doc:`load-balancing <balance>` options that LAMMPS provides.  This is
  because multi-threading achieves parallelism over the number of
  particles, not via their distribution in space.
* A machine is being used in "capability mode", i.e. near the point
  where MPI parallelism is maxed out.  For example, this can happen when
  using the :doc:`PPPM solver <kspace_style>` for long-range
  electrostatics on large numbers of nodes.  The scaling of the KSpace
  calculation (see the :doc:`kspace_style <kspace_style>` command) becomes
  the performance-limiting factor.  Using multi-threading allows less
  MPI tasks to be invoked and can speed-up the long-range solver, while
  increasing overall performance by parallelizing the pairwise and
  bonded calculations via OpenMP.  Likewise additional speedup can be
  sometimes be achieved by increasing the length of the Coulombic cutoff
  and thus reducing the work done by the long-range solver.  Using the
  :doc:`run_style verlet/split <run_style>` command, which is compatible
  with the USER-OMP package, is an alternative way to reduce the number
  of MPI tasks assigned to the KSpace calculation.


Additional performance tips are as follows:

* The best parallel efficiency from *omp* styles is typically achieved
  when there is at least one MPI task per physical CPU chip, i.e. socket
  or die.
* It is usually most efficient to restrict threading to a single
  socket, i.e. use one or more MPI task per socket.
* NOTE: By default, several current MPI implementations use a processor
  affinity setting that restricts each MPI task to a single CPU core.
  Using multi-threading in this mode will force all threads to share the
  one core and thus is likely to be counterproductive.  Instead, binding
  MPI tasks to a (multi-core) socket, should solve this issue.


Restrictions
""""""""""""


None.
