Basics of running LAMMPS
========================

LAMMPS is run from the command line, reading commands from a file via
the -in command line flag, or from standard input.
Using the "-in in.file" variant is recommended:

.. code-block:: bash

   $ lmp_serial -in in.file
   $ lmp_serial < in.file
   $ /path/to/lammps/src/lmp_serial -i in.file
   $ mpirun -np 4 lmp_mpi -in in.file
   $ mpirun -np 8 /path/to/lammps/src/lmp_mpi -in in.file
   $ mpirun -np 6 /usr/local/bin/lmp -in in.file

You normally run the LAMMPS command in the directory where your input
script is located.  That is also where output files are produced by
default, unless you provide specific other paths in your input script or
on the command line.  As in some of the examples above, the LAMMPS
executable itself can be placed elsewhere.

.. note::

   The redirection operator "<" will not always work when running
   in parallel with mpirun; for those systems the -in form is required.

As LAMMPS runs it prints info to the screen and a logfile named
*log.lammps*\ .  More info about output is given on the
:doc:`Run output <Run_output>` doc page.

If LAMMPS encounters errors in the input script or while running a
simulation it will print an ERROR message and stop or a WARNING
message and continue.  See the :doc:`Errors <Errors>` doc page for a
discussion of the various kinds of errors LAMMPS can or can't detect,
a list of all ERROR and WARNING messages, and what to do about them.

----------

LAMMPS can run the same problem on any number of processors, including a
single processor.  In theory you should get identical answers on any
number of processors and on any machine.  In practice, numerical
round-off due to using floating-point math can cause slight differences
and an eventual divergence of molecular dynamics trajectories.  See the
:doc:`Errors common <Errors_common>` doc page for discussion of this.

LAMMPS can run as large a problem as will fit in the physical memory of
one or more processors.  If you run out of memory, you must run on more
processors or define a smaller problem.  The amount of memory needed and
how well it can be distributed across processors may vary based
on the models and settings and commands used.

If you run LAMMPS in parallel via mpirun, you should be aware of the
:doc:`processors <processors>` command, which controls how MPI tasks are
mapped to the simulation box, as well as mpirun options that control how
MPI tasks are assigned to physical cores of the node(s) of the machine
you are running on.  These settings can improve performance, though the
defaults are often adequate.

For example, it is often important to bind MPI tasks (processes) to
physical cores (processor affinity), so that the operating system does
not migrate them during a simulation.  If this is not the default
behavior on your machine, the mpirun option "--bind-to core" (OpenMPI)
or "-bind-to core" (MPICH) can be used.

If the LAMMPS command(s) you are using support multi-threading, you
can set the number of threads per MPI task via the environment
variable OMP_NUM_THREADS, before you launch LAMMPS:

.. code-block:: bash

   $ export OMP_NUM_THREADS=2     # bash
   $ setenv OMP_NUM_THREADS 2     # csh or tcsh

This can also be done via the :doc:`package <package>` command or via
the :doc:`-pk command-line switch <Run_options>` which invokes the
package command.  See the :doc:`package <package>` command or
:doc:`Speed <Speed>` doc pages for more details about which accelerator
packages and which commands support multi-threading.

----------

You can experiment with running LAMMPS using any of the input scripts
provided in the examples or bench directory.  Input scripts are named
in.\* and sample outputs are named log.\*.P where P is the number of
processors it was run on.

Some of the examples or benchmarks require LAMMPS to be built with
optional packages.
