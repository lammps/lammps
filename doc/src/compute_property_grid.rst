.. index:: compute property/grid

compute property/grid command
=============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID property/grid Nx Ny Nz input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* property/grid = style name of this compute command
* Nx, Ny, Nz = grid size in each dimension
* input1,etc = one or more attributes

  .. parsed-literal::

       attributes = id, ix, iy, iz, x, y, z, xs, ys, zs, xc, yc, zc, xsc, ysc, zsc
         id = ID of grid cell, x fastest, y next, z slowest
         ix,iy,iz = grid indices in each dimension (1 to N inclusive)
         x,y,z = coords of lower left corner of grid cell
         xs,ys,zs = scaled coords of lower left corner of grid cell (0.0 to 1.0)
         xc,yc,zc = coords of center point of grid cell
         xsc,ysc,zsc = scaled coords of center point of grid cell (0.0 to 1.0)

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all property/grid id ix iy iz
   compute 1 all property/grid id xc yc zc

Description
"""""""""""

Define a computation that stores the specified attributes of a
distributed grid.  In LAMMPS, distributed grids are regular 2d or 3d
grids which overlay a 2d or 3d simulation domain.  Each processor owns
the grid cells whose center points lie within its subdomain.  See the
:doc:`Howto grid <Howto_grid>` doc page for details of how distributed
grids can be defined by various commands and referenced.

This compute stores the specified attributes of grids as per-grid data
so they can be accessed by other :doc:`output commands <Howto_output>`
such as :doc:`dump grid <dump>`.

*Nx*, *Ny*, and *Nz* define the size of the grid.  For a 2d simulation
*Nz* must be 1.  When this compute is used by :doc:`dump grid <dump>`,
to output per-grid values from other computes of fixes, the grid size
specified for this command must be consistent with the grid sizes
used by the other commands.

The *id* attribute stores the grid ID for each grid cell.  For a
global grid of size Nx by Ny by Nz (in 3d simulations) the grid IDs
range from 1 to Nx*Ny*Nz.  They are ordered with the X index of the 3d
grid varying fastest, then Y, then Z slowest.  For 2d grids (in 2d
simulations), the grid IDs range from 1 to Nx*Ny, with X varying
fastest and Y slowest.

The *ix*, *iy*, *iz* attributes are the indices of a grid cell in
each dimension.  They range from 1 to Nx inclusive in the X dimension,
and similar for Y and Z.

The *x*, *y*, *z* attributes are the coordinates of the lower left
corner point of each grid cell.

The *xs*, *ys*, *zs* attributes are also coordinates of the lower left
corner point of each grid cell, except in scaled coordinates, where
the lower-left corner of the entire simulation box is (0,0,0) and the
upper right corner is (1,1,1).

The *xc*, *yc*, *zc* attributes are the coordinates of the center
point of each grid cell.

The *xsc*, *ysc*, *zsc* attributes are also coordinates of the center
point each grid cell, except in scaled coordinates, where the
lower-left corner of the entire simulation box is (0,0,0) and the upper
right corner is (1,1,1).

For :doc:`triclinic simulation boxes <Howto_triclinic>`, the grid
point coordinates for (x,y,z) and (xc,yc,zc) will reflect the
triclinic geometry.  For (xs,yz,zs) and (xsc,ysc,zsc), the coordinates
are the same for orthogonal versus triclinic boxes.

Output info
"""""""""""

This compute calculates a per-grid vector or array depending on the
number of input values.  The length of the vector or number of array
rows (distributed across all processors) is Nx * Ny * Nz.  For access
by other commands, the name of the single grid produced by this
command is "grid".  The name of its per-grid data is "data".

The (x,y,z) and (xc,yc,zc) coordinates are in distance :doc:`units
<units>`.

Restrictions
""""""""""""

For 2d simulations, the attributes which refer to
the Z dimension cannot be used.

Related commands
""""""""""""""""

:doc:`dump grid <dump>`

Default
"""""""

none
