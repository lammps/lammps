LAMMPS Group Class
******************

The Group class manages named groups of atoms in LAMMPS.  Groups provide a way
to select subsets of atoms for various operations such as applying fixes,
computing properties, or generating output.  Each atom can belong to multiple
groups simultaneously.

Key features include:

- Definition of groups by atom type, region, molecule ID, or other criteria
- Dynamic groups that automatically update membership during simulation
- Efficient bitmask-based membership testing
- Calculation of group properties (count, mass, charge, center of mass, etc.)
- Support for up to 32 groups (limited by 32-bit bitmask)

The special group "all" (index 0) always contains all atoms and cannot be
deleted.  Group membership is stored as a bitmask in each atom's mask property,
allowing efficient O(1) membership testing.

Groups are created and modified using the :doc:`group <group>` command.

.. doxygenclass:: LAMMPS_NS::Group
   :project: progguide
   :members:
