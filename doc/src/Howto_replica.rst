Multi-replica simulations
=========================

Several commands in LAMMPS run multi-replica simulations, meaning
that multiple instances (replicas) of your simulation are run
simultaneously, with small amounts of data exchanged between replicas
periodically.

These are the relevant commands:

* :doc:`hyper <hyper>` for bond boost hyperdynamics (HD)
* :doc:`neb <neb>` for nudged elastic band calculations (NEB)
* :doc:`neb_spin <neb_spin>` for magnetic nudged elastic band calculations
* :doc:`prd <prd>` for parallel replica dynamics (PRD)
* :doc:`tad <tad>` for temperature accelerated dynamics (TAD)
* :doc:`temper <temper>` for parallel tempering with fixed volume
* :doc:`temper/npt <temper_npt>` for parallel tempering extended for NPT
* :doc:`temper/grem <temper_grem>` for parallel tempering with generalized replica exchange (gREM)
* :doc:`fix pimd <fix_pimd>` for path-integral molecular dynamics (PIMD)

NEB is a method for finding transition states and barrier potential energies.
HD, PRD, and TAD are methods for performing accelerated dynamics to find and
perform infrequent events.  Parallel tempering or replica exchange runs
different replicas at a series of temperature to facilitate rare-event
sampling.  PIMD runs different replicas whose individual particles in different
replicas are coupled together by springs to model a system of ring-polymers which
can represent the quantum nature of atom cores.

These commands can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` page for more
info.

In all these cases, you must run with one or more processors per
replica.  The processors assigned to each replica are determined at
run-time by using the :doc:`-partition command-line switch
<Run_options>` to launch LAMMPS on multiple partitions, which in this
context are the same as replicas.  E.g.  these commands:

.. code-block:: bash

   mpirun -np 16 lmp_linux -partition 8x2 -in in.temper
   mpirun -np 8 lmp_linux -partition 8x1 -in in.neb

would each run 8 replicas, on either 16 or 8 processors.  Note the use
of the :doc:`-in command-line switch <Run_options>` to specify the input
script which is required when running in multi-replica mode.

Also note that with MPI installed on a machine (e.g. your desktop), you
can run on more (virtual) processors than you have physical processors.
Thus, the above commands could be run on a single-processor (or
few-processor) desktop so that you can run a multi-replica simulation on
more replicas than you have physical processors. This is useful for
testing and debugging, since with most modern processors and MPI
libraries, the efficiency of a calculation can severely diminish when
oversubscribing processors.
