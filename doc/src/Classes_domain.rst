LAMMPS Domain Class
*******************

The Domain class manages the simulation box geometry and properties in LAMMPS.
It handles both orthogonal and triclinic (tilted) simulation boxes, and provides
the geometric framework for domain decomposition across MPI processes.

Key responsibilities include:

- Defining global simulation box dimensions and boundaries
- Managing periodic, fixed, and shrink-wrap boundary conditions
- Handling triclinic boxes (both restricted and general forms)
- Coordinate transformations between Cartesian and lamda (fractional) coordinates
- Minimum image convention calculations for periodic systems
- Image flag manipulation for atoms crossing periodic boundaries
- Managing simulation regions (geometric shapes for atom selection)
- Managing the lattice for atom creation

The class supports both 2D and 3D simulations and provides methods for
coordinate remapping, unwrapping, and periodic image calculations.  The
simulation box is defined by the :doc:`boundary <boundary>` command and
created by the :doc:`create_box <create_box>` or :doc:`read_data <read_data>`
commands.

.. doxygenclass:: LAMMPS_NS::Domain
   :project: progguide
   :members:
