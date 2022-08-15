.. index:: compute property/grid

compute property/grid command
=============================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID property/grid input1 input2 ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* property/grid = style name of this compute command
* input1,etc = one or more attributes

  .. parsed-literal::

       attributes = id, ix, iy, iz, x, y, z, xs, ys, zs, xc, yc, zc, xsc, ysc, zsc
         id = ID of grid point, x fastest, y next, z slowest
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
the grid points within its subdomain.

See the :doc:`Howto grid <Howto_grid>` doc page for details of how
distributed grids can be defined by various commands and examples of
how they can be used to measure properties of a system.

This compute stores the specified attributes of grids as per-grid data
so they can be accessed by other :doc:`output commands <Howto_output>`
such as :doc:`dump grid <dump>`.

The *id* attribute stores the grid ID for each grid point.  For a
global grid of size Nx by Ny by Nz (in 3d simulations) the grid IDs
range from 1 to Nx*Ny*Nz.  They are ordered with the X index of the 3d
grid varying fastest, then Y, then Z slowest.  For 2d grids (in 2d
simulations), the grid IDs range from 1 to Nx*Ny, with X varying
fastest and Y slowest.

The *ix*, *iy*, *iz* attributes are the indices of a grid point in
each dimension.  They range from 1 to Nx inclusive in the X dimension,
and similar for Y and Z.

The *x*, *y*, *z* attributes are the coordinates of the lower left
corner point of each grid cell.

The *xs*, *ys*, *zs* attributes are also coordinates of the lower left
corner point of each grid cell, except in scaled coordinates, where
the lower-left corner of the entire simulation box is (0,0,0) and he
upper right corner is (1,1,1).


Only in triclinic.

For 2d simulations, none of the attributes which refer to
the Z dimension can be used.



Output info
"""""""""""

This compute calculates a per-grid vector or array depending on the
number of input values.  The length of the vector or number of rows
for each processor is the the number of grid points it owns.

This compute calculates a global vector or global array where the
number of rows = the number of chunks *Nchunk* as calculated by the
specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  If a
single input is specified, a global vector is produced.  If two or
more inputs are specified, a global array is produced where the number
of columns = the number of inputs.  The vector or array can be
accessed by any command that uses global values from a compute as
input.  See the :doc:`Howto output <Howto_output>` page for an
overview of LAMMPS output options.

The vector or array values are "intensive".  The values will be
unitless or in the units discussed above.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix ave/grid <fix_ave_grid>`

Default
"""""""

none
