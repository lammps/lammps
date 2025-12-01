Comm Class
**********

The Comm class is the abstract base class for managing inter-processor
communication in parallel LAMMPS simulations.  It handles the decomposition
of the simulation box across MPI processes and the exchange of atom data
between processors.

Key responsibilities include:

- Domain decomposition across MPI processes (brick or tiled layouts)
- Forward communication of atom coordinates to ghost atoms
- Reverse communication of forces from ghost atoms to owning processors
- Exchange of atoms that migrate between processor domains
- Setting up communication patterns for ghost atom acquisition
- Communication operations for Pair, Bond, Fix, Compute, and Dump classes

This is an abstract base class with two concrete implementations:

- **CommBrick**: Traditional 6-way stencil communication for regular grids
- **CommTiled**: Communication for irregular domain decompositions (from RCB)

The communication style and domain layout can be selected via the
:doc:`comm_style <comm_style>` and :doc:`comm_modify <comm_modify>` commands.

.. doxygenclass:: LAMMPS_NS::Comm
   :project: progguide
   :members:
