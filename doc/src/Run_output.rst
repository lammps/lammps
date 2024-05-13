Screen and logfile output
=========================

As LAMMPS reads an input script, it prints information to both the
screen and a log file about significant actions it takes to setup a
simulation.  When the simulation is ready to begin, LAMMPS performs
various initializations, and prints info about the run it is about to
perform, including the amount of memory (in MBytes per processor) that
the simulation requires.  It also prints details of the initial
thermodynamic state of the system.  During the run itself,
thermodynamic information is printed periodically, every few
timesteps.  When the run concludes, LAMMPS prints the final
thermodynamic state and a total run time for the simulation.  It also
appends statistics about the CPU time and storage requirements for the
simulation.  An example set of statistics is shown here:

.. parsed-literal::

   Loop time of 0.942801 on 4 procs for 300 steps with 2004 atoms

   Performance: 54.985 ns/day, 0.436 hours/ns, 318.201 timesteps/s, 637.674 katom-step/s
   195.2% CPU use with 2 MPI tasks x 2 OpenMP threads

   MPI task timing breakdown:
   Section \|  min time  \|  avg time  \|  max time  \|%varavg\| %total
   ---------------------------------------------------------------
   Pair    \| 0.61419    \| 0.62872    \| 0.64325    \|   1.8 \| 66.69
   Bond    \| 0.0028608  \| 0.0028899  \| 0.002919   \|   0.1 \|  0.31
   Kspace  \| 0.12652    \| 0.14048    \| 0.15444    \|   3.7 \| 14.90
   Neigh   \| 0.10242    \| 0.10242    \| 0.10242    \|   0.0 \| 10.86
   Comm    \| 0.026753   \| 0.027593   \| 0.028434   \|   0.5 \|  2.93
   Output  \| 0.00018341 \| 0.00030942 \| 0.00043542 \|   0.0 \|  0.03
   Modify  \| 0.039117   \| 0.039348   \| 0.039579   \|   0.1 \|  4.17
   Other   \|            \| 0.001041   \|            \|       \|  0.11

   Nlocal:           1002 ave        1006 max         998 min
   Histogram: 1 0 0 0 0 0 0 0 0 1
   Nghost:         8670.5 ave        8691 max        8650 min
   Histogram: 1 0 0 0 0 0 0 0 0 1
   Neighs:         354010 ave      357257 max      350763 min
   Histogram: 1 0 0 0 0 0 0 0 0 1

   Total # of neighbors = 708020
   Ave neighs/atom = 353.30339
   Ave special neighs/atom = 2.3403194
   Neighbor list builds = 26
   Dangerous builds = 0

----------

The first section provides a global loop timing summary. The *loop time*
is the total wall-clock time for the MD steps of the simulation run,
excluding the time for initialization and setup (i.e. the parts that may
be skipped with :doc:`run N pre no <run>`).  The *Performance* line is
provided for convenience to help predict how long it will take to run a
desired physical simulation and to have numbers useful for performance
comparison between different simulation settings or system sizes.  The
*CPU use* line provides the CPU utilization per MPI task; it should be
close to 100% times the number of OpenMP threads (or 1 if not using
OpenMP).  Lower numbers correspond to delays due to file I/O or
insufficient thread utilization from parts of the code that have not
been multi-threaded.

----------

The *MPI task* section gives the breakdown of the CPU run time (in
seconds) into major categories:

* *Pair* = non-bonded force computations
* *Bond* = bonded interactions: bonds, angles, dihedrals, impropers
* *Kspace* = long-range interactions: Ewald, PPPM, MSM
* *Neigh* = neighbor list construction
* *Comm* = inter-processor communication of atoms and their properties
* *Output* = output of thermodynamic info and dump files
* *Modify* = fixes and computes invoked by fixes
* *Other* = all the remaining time

For each category, there is a breakdown of the least, average and most
amount of wall time any processor spent on this category of
computation.  The "%varavg" is the percentage by which the max or min
varies from the average.  This is an indication of load imbalance.  A
percentage close to 0 is perfect load balance.  A large percentage is
imbalance.  The final "%total" column is the percentage of the total
loop time is spent in this category.

When using the :doc:`timer full <timer>` setting, an additional column
is added that also prints the CPU utilization in percent. In addition,
when using *timer full* and the :doc:`package omp <package>` command are
active, a similar timing summary of time spent in threaded regions to
monitor thread utilization and load balance is provided. A new *Thread
timings* section is also added, which lists the time spent in reducing
the per-thread data elements to the storage for non-threaded
computation. These thread timings are measured for the first MPI rank
only and thus, because the breakdown for MPI tasks can change from
MPI rank to MPI rank, this breakdown can be very different for
individual ranks. Here is an example output for this section:

.. parsed-literal::

   Thread timings breakdown (MPI rank 0):
   Total threaded time 0.6846 / 90.6%
   Section \|  min time  \|  avg time  \|  max time  \|%varavg\| %total
   ---------------------------------------------------------------
   Pair    \| 0.5127     \| 0.5147     \| 0.5167     \|   0.3 \| 75.18
   Bond    \| 0.0043139  \| 0.0046779  \| 0.0050418  \|   0.5 \|  0.68
   Kspace  \| 0.070572   \| 0.074541   \| 0.07851    \|   1.5 \| 10.89
   Neigh   \| 0.084778   \| 0.086969   \| 0.089161   \|   0.7 \| 12.70
   Reduce  \| 0.0036485  \| 0.003737   \| 0.0038254  \|   0.1 \|  0.55

----------

The third section above lists the number of owned atoms (Nlocal),
ghost atoms (Nghost), and pairwise neighbors stored per processor.
The max and min values give the spread of these values across
processors with a 10-bin histogram showing the distribution. The total
number of histogram counts is equal to the number of processors.

----------

The last section gives aggregate statistics (across all processors)
for pairwise neighbors and special neighbors that LAMMPS keeps track
of (see the :doc:`special_bonds <special_bonds>` command).  The number
of times neighbor lists were rebuilt is tallied, as is the number of
potentially *dangerous* rebuilds.  If atom movement triggered neighbor
list rebuilding (see the :doc:`neigh_modify <neigh_modify>` command),
then dangerous reneighborings are those that were triggered on the
first timestep atom movement was checked for.  If this count is
non-zero you may wish to reduce the delay factor to ensure no force
interactions are missed by atoms moving beyond the neighbor skin
distance before a rebuild takes place.

----------

If an energy minimization was performed via the
:doc:`minimize <minimize>` command, additional information is printed,
e.g.

.. parsed-literal::

   Minimization stats:
     Stopping criterion = linesearch alpha is zero
     Energy initial, next-to-last, final =
            -6372.3765206     -8328.46998942     -8328.46998942
     Force two-norm initial, final = 1059.36 5.36874
     Force max component initial, final = 58.6026 1.46872
     Final line search alpha, max atom move = 2.7842e-10 4.0892e-10
     Iterations, force evaluations = 701 1516

The first line prints the criterion that determined minimization was
converged. The next line lists the initial and final energy, as well
as the energy on the next-to-last iteration.  The next 2 lines give a
measure of the gradient of the energy (force on all atoms).  The
2-norm is the "length" of this 3N-component force vector; the largest
component (x, y, or z) of force (infinity-norm) is also given.  Then
information is provided about the line search and statistics on how
many iterations and force-evaluations the minimizer required.
Multiple force evaluations are typically done at each iteration to
perform a 1d line minimization in the search direction.  See the
:doc:`minimize <minimize>` page for more details.

----------

If a :doc:`kspace_style <kspace_style>` long-range Coulombics solver
that performs FFTs was used during the run (PPPM, Ewald), then
additional information is printed, e.g.

.. parsed-literal::

   FFT time (% of Kspce) = 0.200313 (8.34477)
   FFT Gflps 3d 1d-only = 2.31074 9.19989

The first line is the time spent doing 3d FFTs (several per timestep)
and the fraction it represents of the total KSpace time (listed
above).  Each 3d FFT requires computation (3 sets of 1d FFTs) and
communication (transposes).  The total flops performed is 5Nlog_2(N),
where N is the number of points in the 3d grid.  The FFTs are timed
with and without the communication and a Gflop rate is computed.  The
3d rate is with communication; the 1d rate is without (just the 1d
FFTs).  Thus you can estimate what fraction of your FFT time was spent
in communication, roughly 75% in the example above.
