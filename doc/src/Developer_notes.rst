Notes for Developers and Code Maintainers
-----------------------------------------

This section documents how a few large sections of code with LAMMPS
work at a conceptual level.  Comments on code in source files
typically document what a variable stores, what a small section of
code does, or what a function does or its input/outputs.  The topics
on this page are intended to document code at a higher level.

KSpace PPPM FFT grids
^^^^^^^^^^^^^^^^^^^^^

The various :doc:`KSpace PPPM <kspace_style>` styles in LAMMPS use
FFTs to solve Poisson's equation.  This subsection describes:

* how FFT grids are defined
* how they are decomposed across processors
* how they are indexed by each processor
* how particle charge and electric field values are mapped to/from
  the grid

An FFT grid cell is a 3d volume; grid points are corners of a grid
cell and the code stores values assigned to grid points in vectors or
3d arrays.  A global 3d FFT grid has points indexed 0 to N-1 inclusive
in each dimension.

Each processor owns two subsets of the grid, each subset is
brick-shaped.  Depending on how it is used, these subsets are
allocated as a 1d vector or 3d array.  Either way, the ordering of
values within contiguous memory x fastest, then y, z slowest.

For the ``3d decomposition`` of the grid, the global grid is
partitioned into bricks that correspond to the subdomains of the
simulation box that each processor owns.  Often, this is a regular 3d
array (Px by Py by Pz) of bricks, where P = number of processors =
Px * Py * Pz.  More generally it can be a tiled decomposition, where
each processor owns a brick and the union of all the bricks is the
global grid.  Tiled decompositions are produced by load balancing with
the RCB algorithm; see the :doc:`balance rcb <balance>` command.

For the ``FFT decompostion`` of the grid, each processor owns a brick
that spans the entire x dimension of the grid while the y and z
dimensions are partitioned as a regular 2d array (P1 by P2), where P =
P1 * P2.

The following indices store the inclusive bounds of the brick a
processor owns, within the global grid:

* nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in = 3d decomposition brick
* nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft = FFT decomposition brick
* nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out = 3d
  decomposition brick + ghost cells

The ``in`` and ``fft`` indices are from 0 to N-1 inclusive in each
dimension, where N is the grid size.

The ``out`` indices index an array which stores the ``in`` subset of
the grid plus ghost cells that surround it.  These indices can thus be
< 0 or >= N.

The number of ghost cells a processor owns in each of the 6 directions
is a function of:

* neighbor skin distance (since atoms can move outside a proc subdomain)
* qdist = offset or charge from atom due to TIP4P fictitious charge
* order = mapping stencil size
* shift = factor used when order is an even number (see below)

Here is an explanation of how the PPPM variables order, nlower/nupper,
shift, and OFFSET work.

These are the relevant variables that determine how atom charge is
mapped to grid points and how field values are mapped from grid points
to atoms:

* order = # of nearby grid points in each dim that atom charge/field are mapped to/from
* nlower,nupper = extent of stencil around the grid point an atom is assigned to
* OFFSET = large integer added/subtracted when mapping to avoid
  int(-0.75) = 0 when -1 is the desired result
  
The particle_map() method assigns each atom to a grid point.

If order is even, say 4:

* atom is assigned to grid point to its left (in each dim)
* shift = OFFSET
* nlower = -1, nupper = 2, which are offsets from assigned grid point
* window of mapping grid pts is thus 2 grid points to left of atom, 2 to right
  
If order is odd, say 5:

* atom is assigned to left/right grid pt it is closest to (in each dim)
* shift = OFFSET + 0.5
* nlower = 2, nupper = 2
* if point is in left half of cell, then
  window of affected grid pts is 3 grid points to left of atom, 2 to right
* if point is in right half of cell, then
  window of affected grid pts is 2 grid points to left of atom, 3 to right

These settings apply to each dimension, so that if order = 5, an
atom's charge is mapped to 125 grid points that surround the atom.
