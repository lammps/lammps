Scatter/gather operations
=========================

This section has functions which gather per-atom data from one or more
processors into a contiguous global list ordered by atom ID.  The same
list is returned to all calling processors.  It also contains
functions which scatter per-atom data from a contiguous global list
across the processors that own those atom IDs.  It also has a
create_atoms() function which can create new atoms by scattering them
appropriately to owning processors in the LAMMPS spatial
decomposition.

It documents the following functions:

- :cpp:func:`lammps_gather_atoms`
- :cpp:func:`lammps_gather_atoms_concat`
- :cpp:func:`lammps_gather_atoms_subset`
- :cpp:func:`lammps_scatter_atoms`
- :cpp:func:`lammps_scatter_atoms_subset`
- :cpp:func:`lammps_gather_bonds`
- :cpp:func:`lammps_gather`
- :cpp:func:`lammps_gather_concat`
- :cpp:func:`lammps_gather_subset`
- :cpp:func:`lammps_scatter`
- :cpp:func:`lammps_scatter_subset`
- :cpp:func:`lammps_create_atoms`

-----------------------

.. doxygenfunction:: lammps_gather_atoms
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_gather_atoms_concat
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_gather_atoms_subset
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_scatter_atoms
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_scatter_atoms_subset
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_gather_bonds
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_gather
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_gather_concat
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_gather_subset
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_scatter
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_scatter_subset
   :project: progguide

-----------------------

.. doxygenfunction:: lammps_create_atoms(void *handle, int n, const int *id, const int *type, const double *x, const double *v, const int *image, int bexpand)
   :project: progguide


