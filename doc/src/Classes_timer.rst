LAMMPS Timer Class
******************

The Timer class provides performance timing for various stages of a LAMMPS
simulation.  It tracks both CPU time and wall-clock time for different
components of the simulation loop.

Key features include:

- Timing of pair, bond, kspace, neighbor, and communication operations
- Multiple levels of timing detail (off, loop, normal, full)
- Optional MPI synchronization before timing for accurate measurements
- Timeout support for limiting total simulation wall time
- Separate timing categories for replica-exchange and NEB operations

Timer detail levels control the amount of timing overhead:

- **OFF**: No timing performed
- **LOOP**: Only total loop time tracked
- **NORMAL**: Major simulation components timed
- **FULL**: All detailed timing including MPI sync points

The timing information is reported at the end of each run and can be used
to identify performance bottlenecks.  Timer settings can be modified using
the :doc:`timer <timer>` command.

.. doxygenclass:: LAMMPS_NS::Timer
   :project: progguide
   :members:
