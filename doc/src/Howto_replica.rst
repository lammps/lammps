Multi-replica simulations
=========================

Several commands in LAMMPS run multi-replica simulations, meaning
that multiple instances (replicas) of your simulation are run
simultaneously, with small amounts of data exchanged between replicas
periodically.

These are the relevant commands:

* :doc:`neb <neb>` for nudged elastic band calculations
* :doc:`neb_spin <neb_spin>` for magnetic nudged elastic band calculations
* :doc:`prd <prd>` for parallel replica dynamics
* :doc:`tad <tad>` for temperature accelerated dynamics
* :doc:`temper <temper>` for parallel tempering
* :doc:`fix pimd <fix_pimd>` for path-integral molecular dynamics (PIMD)

NEB is a method for finding transition states and barrier energies.
PRD and TAD are methods for performing accelerated dynamics to find
and perform infrequent events.  Parallel tempering or replica exchange
runs different replicas at a series of temperature to facilitate
rare-event sampling.

These commands can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

PIMD runs different replicas whose individual particles are coupled
together by springs to model a system or ring-polymers.

This commands can only be used if LAMMPS was built with the USER-MISC
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

In all these cases, you must run with one or more processors per
replica.  The processors assigned to each replica are determined at
run-time by using the :doc:`-partition command-line switch <Run_options>` to launch LAMMPS on multiple partitions,
which in this context are the same as replicas.  E.g.  these commands:

.. code-block:: bash

   mpirun -np 16 lmp_linux -partition 8x2 -in in.temper
   mpirun -np 8 lmp_linux -partition 8x1 -in in.neb

would each run 8 replicas, on either 16 or 8 processors.  Note the use
of the :doc:`-in command-line switch <Run_options>` to specify the input
script which is required when running in multi-replica mode.

Also note that with MPI installed on a machine (e.g. your desktop),
you can run on more (virtual) processors than you have physical
processors.  Thus the above commands could be run on a
single-processor (or few-processor) desktop so that you can run
a multi-replica simulation on more replicas than you have
physical processors.
