Force Class
***********

The Force class manages all force field components in a LAMMPS simulation,
including pair potentials, bonded interactions (bond, angle, dihedral,
improper), and long-range electrostatics (KSpace).  It serves as a container
and factory for these interaction style instances.

Key responsibilities include:

- Creating and managing pair, bond, angle, dihedral, improper, and kspace styles
- Storing physical constants for unit conversions
- Managing Newton's third law settings for pairwise and bonded interactions
- Storing special bonds settings for 1-2, 1-3, 1-4 neighbor exclusions and weights
- Providing factory maps for dynamic creation of force field styles

The class maintains pointers to the currently active force field styles and
provides methods to create, match, and query these styles based on user input
from the :doc:`pair_style <pair_style>`, :doc:`bond_style <bond_style>`,
:doc:`angle_style <angle_style>`, :doc:`dihedral_style <dihedral_style>`,
:doc:`improper_style <improper_style>`, and :doc:`kspace_style <kspace_style>`
commands.

.. doxygenclass:: LAMMPS_NS::Force
   :project: progguide
   :members:
