Fenix online recovery
=====================

The :ref:`Fenix package <PKG-FENIX>` enables the use of new and emerging User
Level Failure Mitigation (ULFM) MPI functionality described at the `ULFM
research hub <https://fault-tolerance.org>`_, in the `Open MPI ULFM
documentation <https://docs.open-mpi.org/en/v5.0.x/features/ulfm.html>`_, and in
:ref:`(Bland2013) <Bland2013>`. The `Fenix library
<https://github.com/sandialabs/fenix>`_ uses the low-level ULFM functionality to
enable high-level online recovery of MPI applications. The benefits of Fenix
are described in :ref:`(Whitlock2022) <Whitlock2022>` and :ref:`(Whitlock2024)
<Whitlock2024>`.

As of writing this text, Fenix requires Open MPI (v5+) for the fault tolerance
features it uses. Newer versions of Open MPI may have significant fault
tolerance performance improvements for some systems. To run Open MPI with fault
tolerance enabled, users must include the :code:`--with-ft=ulfm` flag to the
:code:`mpirun` command.

At a high level, Fenix provides LAMMPS with a resilient communicator that is
repaired after a process is lost. Processes may be lost due to hardware failure,
application bugs, memory starvation, network issues, or any number of other
root causes, all of which are handled in the same manner using Fenix. Fenix
repairs the MPI communicator without requiring a full relaunch of MPI, which
allows applications to avoid the expensive startup costs of MPI (often minutes
per restart on large scale runs). By default, Fenix performs a "shrinking" repair
- meaning the repaired communicator simply excludes the lost process(es) and is
therefore smaller than the pre-failure communicator.

.. code-block:: LAMMPS

   # Example of using shrinking repair with Fenix
   # Print some status information at the start of the file

   if $(is_defined(variable,fenix_restarted)) then &
     "print 'Restarting with $(extract_setting(world_size)) ranks'" &
   else &
     "print 'Starting with $(extract_setting(world_size)) ranks'"

   fenix

   <... simulation setup and running commands ...>

   # Sample output, assuming one rank fails:
   #   Starting with 100 ranks
   #   <... simulation output ...>
   #   Restarting with 99 ranks
   #   <... simulation output ...>

Fenix can avoid this shrinking if the user offers some number of "spare" ranks.
Spare ranks are MPI processes in the Fenix input communicator that are held
aside by Fenix to replace lost processes and avoid shrinking the communicator
during recovery. These spare ranks do nothing and are not seen by the simulator
until they are needed to replace a failed rank. These spare ranks have been
configured to sleep for short periods between checks to see if they are needed
for repair, so they will not utilize resources heavily even if you allocate
multiple spare processes per node.

.. code-block:: LAMMPS

   # Example using spare ranks
   print 'Starting with $(extract_setting(world_size)) ranks'

   fenix spares 10
   print 'Running with $(extract_setting(world_size)) ranks'

   <... simulation setup and running commands ...>

   # Sample output, assuming 1 rank fails and then 10 ranks fail
   #   Starting with 110 ranks
   #   Running with 100 ranks
   #   <... simulation output ...>
   #   Starting with 100 ranks
   #   Running with 100 ranks
   #   <... simulation output ...>
   #   Starting with 99 ranks
   #   Running with 99 ranks

When running LAMMPS with only a single partition, Fenix will leverage the highest N
ranks as spares when N spares are requested. When running with multiple
partitions, Fenix will default to operating on the world comm for each partition
as if it were the only communicator. This will cause problems if your different
partitions communicate with each other over the universe communicator. In that
case, you should use the :code:`universal` argument to initialize Fenix in
universal mode. In universal mode, shrinking recovery is not supported and Fenix
will claim the final partition as spares.

.. code-block:: LAMMPS

   # Example using universal mode
   print 'Partition $(part) starting with $(extract_setting(world_size)) ranks'

   fenix universal
   print 'Partition $(part) running with $(extract_setting(world_size)) ranks'

   <... simulation setup and running commands ...>

   # Sample output assuming 1 rank fails, then another
   #   Partition 0 starting with 10 ranks
   #   Partition 1 starting with 10 ranks
   #   Partition 2 starting with 1 ranks
   #   Partition 0 running with 10 ranks
   #   Partition 1 running with 10 ranks
   #   <... simulation output ...>
   #   Partition 0 starting with 10 ranks
   #   Partition 1 starting with 10 ranks
   #   Partition 0 running with 10 ranks
   #   Partition 1 running with 10 ranks
   #   <... simulation output ...>
   #   ERROR: not enough spare ranks to maintain initial communicator size (src/FENIX/fenix.cpp:217)


Once a failure has been recovered from, the LAMMPS state is entirely wiped.
With the exception of any established Fenix configuration remaining, the
application resumes as if it had just been started for the first time (correctly
accounting for any command-line arguments). By default, this means the LAMMPS
runtime will attempt to rewind the initial input file or STDIN and repeat the
run exactly as before. To establish a different behavior when repairing, users
can configure a recovery jump file and/or label. This is useful to jump past
simulation configuration steps and straight to reading a restart file. It is
valid behavior to invoke the Fenix command again during a restart, but only the
configuration for the recovery jump file and/or label will be updated.


--------

.. _Bland2013:

**(Bland2013)**: Bland, Bouteiller, Herault, Bosilca, and Dongarra (2013). IJHPCA. https://doi.org/10.1177/109434201348823

.. _Whitlock2022:

**(Whitlock2022)** Whitlock, Morales, Bosilca, Bouteiller, Nicolae, Teranishi, Giem, Sarkar (2022). CLUSTER. https://doi.org/10.1109/CLUSTER51413.2022.00052

.. _Whitlock2024:

**(Whitlock2024)** Whitlock (2024). Georgia Institute of Technology. https://doi.org/1853/77831
