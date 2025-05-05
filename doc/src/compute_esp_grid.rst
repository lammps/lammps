.. index:: compute esp/grid
.. index:: compute esp/grid/kk

compute esp/grid/local command
===================================

Accelerator Variants: *esp/grid/kk*

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID esp/grid spacing s

* ID, group-ID are documented in :doc:`compute <compute>` command
* esp/grid = style name of this compute command
* *spacing* = grid spacing (positive real in distance units)


Examples
""""""""

.. code-block:: LAMMPS

    compute esp all esp/grid
    compute esp all esp/grid spacing 0.5

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates electrostatic potential (ESP) on a grid with spacing *s* (default 0.3A). For each grid cell, the ESP is calculated at the center point. Grid points within less than 1.4A of an atom are excluded.

----------

.. include:: accel_styles.rst

----------

Output info
"""""""""""

Compute *esp/grid* evaluates a local array.
The array contains one row for each of the
local grid points, looping over the global index *ix* fastest,
then *iy*, and *iz* slowest.  The array contains math :math:`ntypes+6` columns,
where *ntypes* is the number of LAMMPS types. The first three columns are
the global indexes *ix*, *iy*, and *iz*, followed by the *x*, *y*,
and *z* coordinates of the grid point, followed by the *ntypes* columns
containing the values of the Gaussians for each type.

This compute calculates a per-grid vector of the electrostatic potential in electric field :doc:`units <units>`. The length of the vector is Nx * Ny * Nz where

.. code-block::

    nx=std::max(1,int(std::ceil(domain->prd[0]/spacing)));
    ny=std::max(1,int(std::ceil(domain->prd[1]/spacing)));
    nz=std::max(1,int(std::ceil(domain->prd[2]/spacing)));

For access by other commands, the name of the single grid produced by this command is "grid". The name of its per-grid data is "esp".

Restrictions
""""""""""""

This compute is part of the EXTRA-COMPUTE package.  They are only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute property/grid <compute_property_grid>`, :doc:`dump grid <dump>`

Default
"""""""

The default value for the *spacing* keyword is 0.3A.



