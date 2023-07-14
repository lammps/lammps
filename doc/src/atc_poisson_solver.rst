.. index:: fix_modify AtC poisson_solver

fix_modify AtC poisson_solver command
=====================================

Syntax
""""""

.. code-block:: LAMMPS

   fix_modify <AtC fixID> poisson_solver mesh create <nx> <ny> <nz> <region-ID> <f|p> <f|p> <f|p>

* AtC fixID = ID of :doc:`fix atc <fix_atc>` instance
* poisson_solver = name of the AtC sub-command
* *nx* *ny* *nz* = number of elements in x, y, and z
* region-id = id of region to be meshed
* *f* or *p* = periodicity flags for x, y, and z


Examples
""""""""

.. code-block:: LAMMPS

   fix_modify AtC poisson_solver mesh create 10 1 1 feRegion p p p

Description
"""""""""""

Creates a uniform mesh in a rectangular region.

Restrictions
""""""""""""

Creates only uniform rectangular grids in rectangular regions.

Related AtC commands
""""""""""""""""""""

- :ref:`fix_modify AtC command overview <atc_fix_modify>`

Default
"""""""

None.
