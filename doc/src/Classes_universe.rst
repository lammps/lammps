Universe Class
**************

The Universe class manages the parallel environment for multi-partition LAMMPS
simulations.  When LAMMPS is run with the ``-partition`` command-line flag,
the total MPI communicator is divided into separate worlds, each running an
independent simulation.  These partitions can communicate with each other
through the universe communicator for coordinated multi-replica methods such
as replica exchange, parallel tempering, or nudged elastic band calculations.

For single-partition runs (the common case), the universe communicator is
equivalent to the world communicator, and there is only one world.

.. doxygenclass:: LAMMPS_NS::Universe
   :project: progguide
   :members:
