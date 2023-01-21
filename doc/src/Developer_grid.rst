Use of distributed grids within style classes
---------------------------------------------

.. versionadded:: 22Dec2022

The LAMMPS source code includes two classes which facilitate the
creation and use of distributed grids.  These are the Grid2d and
Grid3d classes in the src/grid2d.cpp.h and src/grid3d.cpp.h files
respectively.  As the names imply, they are used for 2d or 3d
simulations, as defined by the :doc:`dimension <dimension>` command.

The :doc:`Howto_grid <Howto_grid>` page gives an overview of how
distributed grids are defined from a user perspective, lists LAMMPS
commands which use them, and explains how grid cell data is referenced
from an input script.  Please read that page first as it motivates the
coding details discussed here.

This doc page is for users who wish to write new styles (input script
commands) which use distributed grids.  There are a variety of
material models and analysis methods which use atoms (or
coarse-grained particles) and grids in tandem.

A *distributed* grid means each processor owns a subset of the grid
cells.  In LAMMPS, the subset for each processor will be a sub-block
of grid cells with low and high index bounds in each dimension of the
grid.  The union of the sub-blocks across all processors is the global
grid.

More specifically, a grid point is defined for each cell (by default
the center point), and a processor owns a grid cell if its point is
within the processor's spatial sub-domain.  The union of processor
sub-domains is the global simulation box.  If a grid point is on the
boundary of two sub-domains, the lower processor owns the grid cell.  A
processor may also store copies of ghost cells which surround its
owned cells.

----------

Style commands
^^^^^^^^^^^^^^

Style commands which can define and use distributed grids include the
:doc:`compute <compute>`, :doc:`fix <fix>`, :doc:`pair <pair_style>`,
and :doc:`kspace <kspace_style>` styles.  If you wish grid cell data
to persist across timesteps, then use a fix.  If you wish grid cell
data to be accessible by other commands, then use a fix or compute.
Currently in LAMMPS, the :doc:`pair_style amoeba <pair_amoeba>`,
:doc:`kspace_style pppm <kspace_style>`, and :doc:`kspace_style msm
<kspace_style>` commands use distributed grids but do not require
either of these capabilities; they thus create and use distributed
grids internally.  Note that a pair style which needs grid cell data
to persist could be coded to work in tandem with a fix style which
provides that capability.

The *size* of a grid is specified by the number of grid cells in each
dimension of the simulation domain.  In any dimension the size can be
any value >= 1.  Thus a 10x10x1 grid for a 3d simulation is
effectively a 2d grid, where each grid cell spans the entire
z-dimension.  A 1x100x1 grid for a 3d simulation is effectively a 1d
grid, where grid cells are a series of thin xz slabs in the
y-dimension.  It is even possible to define a 1x1x1 3d grid, though it
may be inefficient to use it in a computational sense.

Note that the choice of grid size is independent of the number of
processors or their layout in a grid of processor sub-domains which
overlays the simulations domain.  Depending on the distributed grid
size, a single processor may own many 1000s or no grid cells.

A command can define multiple grids, each of a different size.  Each
grid is an instantiation of the Grid2d or Grid3d class.

The command also defines what data it will store for each grid it
creates and it allocates the multi-dimensional array(s) needed to
store the data.  No grid cell data is stored within the Grid2d or
Grid3d classes.

If a single value per grid cell is needed, the data array will have
the same dimension as the grid, i.e. a 2d array for a 2d grid,
likewise for 3d.  If multiple values per grid cell are needed, the
data array will have one more dimension than the grid, i.e. a 3d array
for a 2d grid, or 4d array for a 3d grid.  A command can choose to
define multiple data arrays for each grid it defines.

----------

Grid data allocation and access
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest way for a command to allocate and access grid cell data
is to use the *create_offset()* methods provided by the Memory class.
Arguments for these methods can be values returned by the
*setup_grid()* method (described below), which define the extent of
the grid cells (owned+ghost) the processor owns.  These 4 methods
allocate memory for 2d (first two) and 3d (second two) grid data.  The
two methods that end in "_one" allocate an array which stores a single
value per grid cell.  The two that end in "_multi" allocate an array
which stores *Nvalues* per grid cell.

.. code-block:: c++

   // single value per cell for a 2d grid = 2d array
   memory->create2d_offset(data2d_one, nylo_out, nyhi_out,
                           nxlo_out, nxhi_out, "data2d_one");

   // nvalues per cell for a 2d grid = 3d array
   memory->create3d_offset_last(data2d_multi, nylo_out, nyhi_out,
                                nxlo_out, nxhi_out, nvalues, "data2d_multi");

   // single value per cell for a 3d grid = 3d array
   memory->create3d_offset(data3d_one, nzlo_out, nzhi_out, nylo_out,
                           nyhi_out, nxlo_out, nxhi_out, "data3d_one");

   // nvalues per cell for a 3d grid = 4d array
   memory->create4d_offset_last(data3d_multi, nzlo_out, nzhi_out, nylo_out,
                                nyhi_out, nxlo_out, nxhi_out, nvalues,
                                "data3d_multi");

Note that these multi-dimensional arrays are allocated as contiguous
chunks of memory where the x-index of the grid varies fastest, then y,
and the z-index slowest.  For multiple values per grid cell, the
Nvalues are contiguous, so their index varies even faster than the
x-index.

The key point is that the "offset" methods create arrays which are
indexed by the range of indices which are the bounds of the sub-block
of the global grid owned by this processor.  This means loops like
these can be written in the caller code to loop over owned grid cells,
where the "i" loop bounds are the range of owned grid cells for the
processor.  These are the bounds returned by the *setup_grid()*
method:

.. code-block:: c++

    for (int iy = iylo; iy <= iyhi; iy++)
      for (int ix = ixlo; ix <= ixhi; ix++)
        data2d_one[iy][ix] = 0.0;

    for (int iy = iylo; iy <= iyhi; iy++)
      for (int ix = ixlo; ix <= ixhi; ix++)
        for (int m = 0; m < nvalues; m++)
          data2d_multi[iy][ix][m] = 0.0;

    for (int iz = izlo; iz <= izhi; iz++)
      for (int iy = iylo; iy <= iyhi; iy++)
        for (int ix = ixlo; ix <= ixhi; ix++)
          data3d_one[iz][iy][ix] = 0.0;

    for (int iz = izlo; iz <= izhi; iz++)
      for (int iy = iylo; iy <= iyhi; iy++)
        for (int ix = ixlo; ix <= ixhi; ix++)
           for (int m = 0; m < nvalues; m++)
              data3d_multi[iz][iy][ix][m] = 0.0;

Simply replacing the "i" bounds with "o" bounds, also returned by the
*setup_grid()* method, would alter this code to loop over owned+ghost
cells (the entire allocated grid).

----------

Grid class constructors
^^^^^^^^^^^^^^^^^^^^^^^

The following sub-sections describe the public methods of the Grid3d
class which a style command can invoke.  The Grid2d methods are
similar; simply remove arguments which refer to the z-dimension.

There are 2 constructors which can be used.  They differ in the extra
i/o xyz lo/hi arguments:

.. code-block:: c++

   Grid3d(class LAMMPS *lmp, MPI_Comm gcomm, int gnx, int gny, int gnz)
   Grid3d(class LAMMPS *lmp, MPI_Comm gcomm, int gnx, int gny, int gnz,
          int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
          int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi)

Both constructors take the LAMMPS instance pointer and a communicator
over which the grid will be distributed.  Typically this is the
*world* communicator the LAMMPS instance is using.  The
:doc:`kspace_style msm <kspace_style>` command creates a series of
grids, each of different size, which are partitioned across different
sub-communicators of processors.  Both constructors are also passed
the global grid size: *gnx* by *gny* by *gnz*.

The first constructor is used when the caller wants the Grid class to
partition the global grid across processors; the Grid class defines
which grid cells each processor owns and also which it stores as ghost
cells.  A subsequent call to *setup_grid()*, discussed below, returns
this info to the caller.

The second constructor allows the caller to define the extent of owned
and ghost cells, and pass them to the Grid class.  The 6 arguments
which start with "i" are the inclusive lower and upper index bounds of
the owned (inner) grid cells this processor owns in each of the 3
dimensions within the global grid.  Owned grid cells are indexed from
0 to N-1 in each dimension.

The 6 arguments which start with "o" are the inclusive bounds of the
owned+ghost (outer) grid cells it stores.  If the ghost cells are on
the other side of a periodic boundary, then these indices may be < 0
or >= N in any dimension, so that oxlo <= ixlo and ixhi >= ixhi is
always the case.

For example, if Nx = 100, then a processor might pass ixlo=50,
ixhi=60, oxlo=48, oxhi=62 to the Grid class.  Or ixlo=0, ixhi=10,
oxlo=-2, oxhi=13.  If a processor owns no grid cells in a dimension,
then the ihi value should be specified as one less than the ilo value.

Note that the only reason to use the second constructor is if the
logic for assigning ghost cells is too complex for the Grid class to
compute, using the various set() methods described next.  Currently
only the kspace_style pppm/electrode and kspace_style msm commands use
the second constructor.

----------

Grid class set methods
^^^^^^^^^^^^^^^^^^^^^^

The following methods affect how the Grid class computes which owned
and ghost cells are assigned to each processor.  *Set_shift_grid()* is
the only method which influences owned cell assignment; all the rest
influence ghost cell assignment.  These methods are only used with the
first constructor; they are ignored if the second constructor is used.
These methods must be called before the *setup_grid()* method is
invoked, because they influence its operation.

.. code-block:: c++

   void set_shift_grid(double shift);
   void set_distance(double distance);
   void set_stencil_atom(int lo, int hi);
   void set_shift_atom(double shift_lo, double shift_hi);
   void set_stencil_grid(int lo, int hi);
   void set_zfactor(double factor);

Processors own a grid cell if a point within the grid cell is inside
the processor's sub-domain.  By default this is the center point of the
grid cell.  The *set_shift_grid()* method can change this.  The *shift*
argument is a value from 0.0 to 1.0 (inclusive) which is the offset of
the point within the grid cell in each dimension.  The default is 0.5
for the center of the cell.  A value of 0.0 is the lower left corner
point; a value of 1.0 is the upper right corner point.  There is
typically no need to change the default as it is optimal for
minimizing the number of ghost cells needed.

If a processor maps its particles to grid cells, it needs to allow for
its particles being outside its sub-domain between reneighboring.  The
*distance* argument of the *set_distance()* method sets the furthest
distance outside a processor's sub-domain which a particle can move.
Typically this is half the neighbor skin distance, assuming
reneighboring is done appropriately.  This distance is used in
determining how many ghost cells a processor needs to store to enable
its particles to be mapped to grid cells.  The default value is 0.0.

Some commands, like the :doc:`kspace_style pppm <kspace_style>`
command, map values (charge in the case of PPPM) to a stencil of grid
cells beyond the grid cell the particle is in.  The stencil extent may
be different in the low and high directions.  The *set_stencil_atom()*
method defines the maximum values of those 2 extents, assumed to be
the same in each of the 3 dimensions.  Both the lo and hi values are
specified as positive integers.  The default values are both 0.

Some commands, like the :doc:`kspace_style pppm <kspace_style>`
command, shift the position of an atom when mapping it to a grid cell,
based on the size of the stencil used to map values to the grid
(charge in the case of PPPM).  The lo and hi arguments of the
*set_shift_atom()* method are the minimum shift in the low direction
and the maximum shift in the high direction, assumed to be the same in
each of the 3 dimensions.  The shifts should be fractions of a grid
cell size with values between 0.0 and 1.0 inclusive.  The default
values are both 0.0.  See the src/pppm.cpp file for examples of these
lo/hi values for regular and staggered grids.

Some methods like the :doc:`fix ttm/grid <fix_ttm>` command, perform
finite difference kinds of operations on the grid, to diffuse electron
heat in the case of the two-temperature model (TTM).  This operation
uses ghost grid values beyond the owned grid values the processor
updates.  The *set_stencil_grid()* method defines the extent of this
stencil in both directions, assumed to be the same in each of the 3
dimensions.  Both the lo and hi values are specified as positive
integers.  The default values are both 0.

The kspace_style pppm commands allow a grid to be defined which
overlays a volume which extends beyond the simulation box in the z
dimension.  This is for the purpose of modeling a 2d-periodic slab
(non-periodic in z) as if it were a larger 3d periodic system,
extended (with empty space) in the z dimension.  The
:doc:`kspace_modify slab <kspace_modify>` command is used to specify
the ratio of the larger volume to the simulation volume; a volume
ratio of ~3 is typical.  For this kind of model, the PPPM caller sets
the global grid size *gnz* ~3x larger than it would be otherwise.
This same ratio is passed by the PPPM caller as the *factor* argument
to the Grid class via the *set_zfactor()* method (*set_yfactor()* for
2d grids).  The Grid class will then assign ownership of the 1/3 of
grid cells that overlay the simulation box to the processors which
also overlay the simulation box.  The remaining 2/3 of the grid cells
are assigned to processors whose sub-domains are adjacent to the upper
z boundary of the simulation box.

----------

Grid class setup_grid method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The *setup_grid()* method is called after the first constructor
(above) to partition the grid across processors, which determines
which grid cells each processor owns.  It also calculates how many
ghost grid cells in each dimension and each direction each processor
needs to store.

Note that this method is NOT called if the second constructor above is
used.  In that case, the caller assigns owned and ghost cells to each
processor.

Also note that this method must be invoked after any *set_*()* methods have
been used, since they can influence the assignment of owned and ghost
cells.

.. code-block:: c++

   void setup_grid(int &ixlo, int &ixhi, int &iylo, int &iyhi, int &izlo, int &izhi,
                   int &oxlo, int &oxhi, int &oylo, int &oyhi, int &ozlo, int &ozhi)

The 6 return arguments which start with "i" are the inclusive lower
and upper index bounds of the owned (inner) grid cells this processor
owns in each of the 3 dimensions within the global grid.  Owned grid
cells are indexed from 0 to N-1 in each dimension.

The 6 return arguments which start with "o" are the inclusive bounds of
the owned+ghost cells it owns.  If the ghost cells are on the other
side of a periodic boundary, then these indices may be < 0 or >= N in
any dimension, so that oxlo <= ixlo and ixhi >= ixhi is always the
case.

----------

More grid class set methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following 2 methods can be used to override settings made by the
constructors above.  If used, they must be called called before the
*setup_comm()* method is invoked, since it uses the settings that
these methods override.  In LAMMPS these methods are called by by the
:doc:`kspace_style msm <kspace_style>` command for the grids it
instantiates using the 2nd constructor above.

.. code-block:: c++

   void set_proc_neighs(int pxlo, int pxhi, int pylo, int pyhi, int pzlo, int pzhi)
   void set_caller_grid(int fxlo, int fxhi, int fylo, int fyhi, int fzlo, int fzhi)

The *set_proc_neighs()* method sets the processor IDs of the 6
neighboring processors for each processor.  Normally these would match
the processor grid neighbors which LAMMPS creates to overlay the
simulation box (the default).  However, MSM excludes non-participating
processors from coarse grid communication when less processors are
used.  This method allows MSM to override the default values.

The *set_caller_grid()* method species the size of the data arrays the
caller allocates.  Normally these would match the extent of the ghost
grid cells (the default).  However the MSM caller allocates a larger
data array (more ghost cells) for its finest-level grid, for use in
other operations besides owned/ghost cell communication.  This method
allows MSM to override the default values.


----------

Grid class get methods
^^^^^^^^^^^^^^^^^^^^^^

The following methods allow the caller to query the settings for a
specific grid, whether it created the grid or another command created
it.

.. code-block:: c++

   void get_size(int &nxgrid, int &nygrid, int &nzgrid);
   void get_bounds_owned(int &xlo, int &xhi, int &ylo, int &yhi, int &zlo, int &zhi)
   void get_bounds_ghost(int &xlo, int &xhi, int &ylo, int &yhi, int &zlo, int &zhi)

The *get_size()* method returns the size of the global grid in each dimension.

The *get_bounds_owned()* method return the inclusive index bounds of
the grid cells this processor owns.  The values range from 0 to N-1 in
each dimension.  These values are the same as the "i" values returned
by *setup_grid()*.

The *get_bounds_ghost()* method return the inclusive index bounds of
the owned+ghost grid cells this processor stores.  The owned cell
indices range from 0 to N-1, so these indices may be less than 0 or
greater than or equal to N in each dimension.  These values are the
same as the "o" values returned by *setup_grid()*.

----------

Grid class owned/ghost communication
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If needed by the command, the following methods setup and perform
communication of grid data to/from neighboring processors.  The
*forward_comm()* method sends owned grid cell data to the
corresponding ghost grid cells on other processors.  The
*reverse_comm()* method sends ghost grid cell data to the
corresponding owned grid cells on another processor.  The caller can
choose to sum ghost grid cell data to the owned grid cell or simply
copy it.

.. code-block:: c++

   void setup_comm(int &nbuf1, int &nbuf2)
   void forward_comm(int caller, void *ptr, int which, int nper, int nbyte,
                     void *buf1, void *buf2, MPI_Datatype datatype);
   void reverse_comm(int caller, void *ptr, int which, int nper, int nbyte,
                     void *buf1, void *buf2, MPI_Datatype datatype)
   int ghost_adjacent();

The *setup_comm()* method must be called one time before performing
*forward* or *reverse* communication (multiple times if needed).  It
returns two integers, which should be used to allocate two buffers.
The *nbuf1* and *nbuf2* values are the number of grid cells whose data
will be stored in two buffers by the Grid class when *forward* or
*reverse* communication is performed.  The caller should thus allocate
them to a size large enough to hold all the data used in any single
forward or reverse communication operation it performs.  Note that the
caller may allocate and communicate multiple data arrays for a grid it
instantiates.  This size includes the bytes needed for the data type
of the grid data it stores, e.g. double precision values.

The *forward_comm()* and *reverse_comm()* methods send grid cell data
from owned to ghost cells, or ghost to owned cells, respectively, as
described above.  The *caller* argument should be one of these values
-- Grid3d::COMPUTE, Grid3d::FIX, Grid3d::KSPACE, Grid3d::PAIR --
depending on the style of the caller class.  The *ptr* argument is the
"this" pointer to the caller class.  These 2 arguments are used to
call back to pack()/unpack() functions in the caller class, as
explained below.

The *which* argument is a flag the caller can set which is passed to
the caller's pack()/unpack() methods.  This allows a single callback
method to pack/unpack data for several different flavors of
forward/reverse communication, e.g. operating on different grids or
grid data.

The *nper* argument is the number of values per grid cell to be
communicated.  The *nbyte* argument is the number of bytes per value,
e.g. 8 for double-precision values.  The *buf1* and *buf2* arguments
are the two allocated buffers described above.  So long as they are
allocated for the maximum size communication, they can be re-used for
any *forward_comm()/reverse_comm()* call.  The *datatype* argument is
the MPI_Datatype setting, which should match the buffer allocation and
the *nbyte* argument.  E.g. MPI_DOUBLE for buffers storing double
precision values.

To use the *forward_grid()* method, the caller must provide two
callback functions; likewise for use of the *reverse_grid()* methods.
These are the 4 functions, their arguments are all the same.

.. code-block:: c++

   void pack_forward_grid(int which, void *vbuf, int nlist, int *list);
   void unpack_forward_grid(int which, void *vbuf, int nlist, int *list);
   void pack_reverse_grid(int which, void *vbuf, int nlist, int *list);
   void unpack_reverse_grid(int which, void *vbuf, int nlist, int *list);

The *which* argument is set to the *which* value of the
*forward_comm()* or *reverse_comm()* calls.  It allows the pack/unpack
function to select what data values to pack/unpack.  *Vbuf* is the
buffer to pack/unpack the data to/from.  It is a void pointer so that
the caller can cast it to whatever data type it chooses, e.g. double
precision values.  *Nlist* is the number of grid cells to pack/unpack
and *list* is a vector (nlist in length) of offsets to where the data
for each grid cell resides in the caller's data arrays, which is best
illustrated with an example from the src/EXTRA-FIX/fix_ttm_grid.cpp
class which stores the scalar electron temperature for 3d system in a
3d grid (one value per grid cell):

.. code-block:: c++

   void FixTTMGrid::pack_forward_grid(int /*which*/, void *vbuf, int nlist, int *list)
   {
     auto buf = (double *) vbuf;
     double *src = &T_electron[nzlo_out][nylo_out][nxlo_out];
     for (int i = 0; i < nlist; i++) buf[i] = src[list[i]];
   }

In this case, the *which* argument is not used, *vbuf* points to a
buffer of doubles, and the electron temperature is stored by the
FixTTMGrid class in a 3d array of owned+ghost cells called T_electron.
That array is allocated by the *memory->create_3d_offset()* method
described above so that the first grid cell it stores is indexed as
T_electron[nzlo_out][nylo_out][nxlo_out].  The *nlist* values in
*list* are integer offsets from that first grid cell.  Setting *src*
to the address of the first cell allows those offsets to be used to
access the temperatures to pack into the buffer.

Here is a similar portion of code from the src/fix_ave_grid.cpp class
which can store two kinds of data, a scalar count of atoms in a grid
cell, and one or more grid-cell-averaged atom properties.  The code
from its *unpack_reverse_grid()* function for 2d grids and multiple
per-atom properties per grid cell (*nvalues*) is shown here:

.. code-block:: c++

   void FixAveGrid::unpack_reverse_grid(int /*which*/, void *vbuf, int nlist, int *list)
   {
     auto buf = (double *) vbuf;
     double *count,*data,*values;
     count = &count2d[nylo_out][nxlo_out];
     data = &array2d[nylo_out][nxlo_out][0];
     m = 0;
     for (i = 0; i < nlist; i++) {
       count[list[i]] += buf[m++];
       values = &data[nvalues*list[i]];
       for (j = 0; j < nvalues; j++)
        values[j] += buf[m++];
     }
   }

Both the count and the multiple values per grid cell are communicated
in *vbuf*.  Note that *data* is now a pointer to the first value in
the first grid cell.  And *values* points to where the first value in
*data* is for an offset of grid cells, calculated by multiplying
*nvalues* by *list[i]*.  Finally, because this is reverse
communication, the communicated buffer values are summed to the caller
values.

The *ghost_adjacent()* method returns a 1 if every processor can
perform the necessary owned/ghost communication with only its nearest
neighbor processors (4 in 2d, 6 in 3d).  It returns a 0 if any
processor's ghost cells extend further than nearest neighbor
processors.

This can be checked by callers who have the option to change the
global grid size to insure more efficient nearest-neighbor-only
communication if they wish.  In this case, they instantiate a grid of
a given size (resolution), then invoke *setup_comm()* followed by
*ghost_adjacent()*.  If the ghost cells are not adjacent, they destroy
the grid instance and start over with a higher-resolution grid.
Several of the :doc:`kspace_style pppm <kspace_style>` command
variants have this option.

----------

Grid class remap methods for load balancing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following methods are used when a load-balancing operation,
triggered by the :doc:`balance <balance>` or :doc:`fix balance
<fix_balance>` commands, changes the partitioning of the simulation
domain into processor sub-domains.

In order to work with load-balancing, any style command (compute, fix,
pair, or kspace style) which allocates a grid and stores per-grid data
should define a *reset_grid()* method; it takes no arguments.  It will
be called by the two balance commands after they have reset processor
sub-domains and migrated atoms (particles) to new owning processors.
The *reset_grid()* method will typically perform some or all of the
following operations.  See the src/fix_ave_grid.cpp and
src/EXTRA_FIX/fix_ttm_grid.cpp files for examples of *reset_grid()*
methods, as well as the *pack_remap_grid()* and *unpack_remap_grid()*
functions.

First, the *reset_grid()* method can instantiate new grid(s) of the
same global size, then call *setup_grid()* to partition them via the
new processor sub-domains.  At this point, it can invoke the
*identical()* method which compares the owned and ghost grid cell
index bounds between two grids, the old grid passed as a pointer
argument, and the new grid whose *identical()* method is being called.
It returns 1 if the indices match on all processors, otherwise 0.  If
they all match, then the new grids can be deleted; the command can
continue to use the old grids.

If not, then the command should allocate new grid data array(s) which
depend on the new partitioning.  If the command does not need to
persist its grid data from the old partitioning to the new one, then
the command can simply delete the old data array(s) and grid
instance(s).  It can then return.

If the grid data does need to persist, then the data for each grid
needs to be "remapped" from the old grid partitioning to the new grid
partitioning.  The *setup_remap()* and *remap()* methods are used for
that purpose.

.. code-block:: c++

   int identical(Grid3d *old);
   void setup_remap(Grid3d *old, int &nremap_buf1, int &nremap_buf2)
   void remap(int caller, void *ptr, int which, int nper, int nbyte,
              void *buf1, void *buf2, MPI_Datatype datatype)

The arguments to these methods are identical to those for
the *setup_comm()* and *forward_comm()* or *reverse_comm()* methods.
However the returned *nremap_buf1* and *nremap2_buf* values will be
different than the *nbuf1* and *nbuf2* values.  They should be used to
allocate two different remap buffers, separate from the owned/ghost
communication buffers.

To use the *remap()* method, the caller must provide two
callback functions:

.. code-block:: c++

   void pack_remap_grid(int which, void *vbuf, int nlist, int *list);
   void unpack_remap_grid(int which, void *vbuf, int list, int *list);

Their arguments are identical to those for the *pack_forward_grid()*
and *unpack_forward_grid()* callback functions (or the reverse
variants) discussed above.  Normally, both these methods pack/unpack
all the data arrays for a given grid.  The *which* argument of the
*remap()* method sets the *which* value for the pack/unpack functions.
If the command instantiates multiple grids (of different sizes), it
can be used within the pack/unpack methods to select which grid's data
is being remapped.

Note that the *pack_remap_grid()* function must copy values from the
OLD grid data arrays into the *vbuf* buffer. The *unpack_remap_grid()*
function must copy values from the *vbuf* buffer into the NEW grid
data arrays.

After the remap operation for grid cell data has been performed, the
*reset_grid()* method can deallocate the two remap buffers it created,
and can then exit.

----------

Grid class I/O methods
^^^^^^^^^^^^^^^^^^^^^^

There are two I/O methods in the Grid classes which can be used to
read and write grid cell data to files.  The caller can decide on the
precise format of each file, e.g. whether header lines are prepended
or comment lines are allowed.  Fundamentally, the file should contain
one line per grid cell for the entire global grid.  Each line should
contain identifying info as to which grid cell it is, e.g. a unique
grid cell ID or the ix,iy,iz indices of the cell within a 3d grid.
The line should also contain one or more data values which are stored
within the grid data arrays created by the command

For grid cell IDs, the LAMMPS convention is that the IDs run from 1 to
N, where N = Nx * Ny for 2d grids and N = Nx * Ny * Nz for 3d grids.
The x-index of the grid cell varies fastest, then y, and the z-index
varies slowest.  So for a 10x10x10 grid the cell IDs from 901-1000
would be in the top xy layer of the z dimension.

The *read_file()* method does something simple.  It reads a chunk of
consecutive lines from the file and passes them back to the caller to
process.  The caller provides a *unpack_read_grid()* function for this
purpose.  The function checks the grid cell ID or indices and only
stores grid cell data for the grid cells it owns.

The *write_file()* method does something slightly more complex.  Each
processor packs the data for its owned grid cells into a buffer.  The
caller provides a *pack_write_grid()* function for this purpose.  The
*write_file()* method then loops over all processors and each sends
its buffer one at a time to processor 0, along with the 3d (or 2d)
index bounds of its grid cell data within the global grid.  Processor
0 calls back to the *unpack_write_grid()* function provided by the
caller with the buffer.  The function writes one line per grid cell to
the file.

See the src/EXTRA_FIX/fix_ttm_grid.cpp file for examples of now both
these methods are used to read/write electron temperature values
from/to a file, as well as for implementations of the the pack/unpack
functions described below.

Here are the details of the two I/O methods and the 3 callback
functions.  See the src/fix_ave_grid.cpp file for examples of all of
them.

.. code-block:: c++

   void read_file(int caller, void *ptr, FILE *fp, int nchunk, int maxline)
   void write_file(int caller, void *ptr, int which,
                   int nper, int nbyte, MPI_Datatype datatype

The *caller* argument in both methods should be one of these values --
Grid3d::COMPUTE, Grid3d::FIX, Grid3d::KSPACE, Grid3d::PAIR --
depending on the style of the caller class.  The *ptr* argument in
both methods is the "this" pointer to the caller class.  These 2
arguments are used to call back to pack()/unpack() functions in the
caller class, as explained below.

For the *read_file()* method, the *fp* argument is a file pointer to
the file to be read from, opened on processor 0 by the caller.
*Nchunk* is the number of lines to read per chunk, and *maxline* is
the maximum number of characters per line.  The Grid class will
allocate a buffer for storing chunks of lines based on these values.

For the *write_file()* method, the *which* argument is a flag the
caller can set which is passed back to the caller's pack()/unpack()
methods.  If the command instantiates multiple grids (of different
sizes), this flag can be used within the pack/unpack methods to select
which grid's data is being written out (presumably to different
files).  the *nper* argument is the number of values per grid cell to
be written out.  The *nbyte* argument is the number of bytes per
value, e.g. 8 for double-precision values.  The *datatype* argument is
the MPI_Datatype setting, which should match the *nbyte* argument.
E.g. MPI_DOUBLE for double precision values.

To use the *read_grid()* method, the caller must provide one callback
function.  To use the *write_grid()* method, it provides two callback
functions:

.. code-block:: c++

   int unpack_read_grid(int nlines, char *buffer)
   void pack_write_grid(int which, void *vbuf)
   void unpack_write_grid(int which, void *vbuf, int *bounds)

For *unpack_read_grid()* the *nlines* argument is the number of lines
of character data read from the file and contained in *buffer*.  The
lines each include a newline character at the end.  When the function
processes the lines, it may choose to skip some of them (header or
comment lines).  It returns an integer count of the number of grid
cell lines it processed.  This enables the Grid class *read_file()*
method to know when it has read the correct number of lines.

For *pack_write_grid()* and *unpack_write_grid()*, the *vbuf* argument
is the buffer to pack/unpack data to/from.  It is a void pointer so
that the caller can cast it to whatever data type it chooses,
e.g. double precision values.  the *which* argument is set to the
*which* value of the *write_file()* method.  It allows the caller to
choose which grid data to operate on.

For *unpack_write_grid()*, the *bounds* argument is a vector of 4 or 6
integer grid indices (4 for 2d, 6 for 3d).  They are the
xlo,xhi,ylo,yhi,zlo,zhi index bounds of the portion of the global grid
which the *vbuf* holds owned grid cell data values for.  The caller
should loop over the values in *vbuf* with a double loop (2d) or
triple loop (3d), similar to the code snippets listed above.  The
x-index varies fastest, then y, and the z-index slowest.  If there are
multiple values per grid cell, the index for those values varies
fastest of all.  The caller can add the x,y,z indices of the grid cell
(or the corresponding grid cell ID) to the data value(s) written as
one line to the output file.

----------

Style class grid access methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A style command can enable its grid cell data to be accessible from
other commands.  For example :doc:`fix ave/grid <fix_ave_grid>` or
:doc:`dump grid <dump>` or :doc:`dump grid/vtk <dump>`.  Those
commands access the grid cell data by using a *grid reference* in
their input script syntax, as described on the :doc:`Howto_grid
<Howto_grid>` doc page.  They look like this:

* c_ID:gname:dname
* c_ID:gname:dname[I]
* f_ID:gname:dname
* f_ID:gname:dname[I]

Each grid a command instantiates has a unique *gname*, defined by the
command.  Likewise each grid cell data structure (scalar or vector)
associated with a grid has a unique *dname*, also defined by the
command.

To provide access to its grid cell data, a style command needs to
implement the following 4 methods:

.. code-block:: c++

   int get_grid_by_name(const std::string &name, int &dim);
   void *get_grid_by_index(int index);
   int get_griddata_by_name(int igrid, const std::string &name, int &ncol);
   void *get_griddata_by_index(int index);

Currently only computes and fixes can implement these methods.  If it
does so, the compute of fix should also set the variable
*pergrid_flag* to 1.  See any of the compute or fix commands which set
"pergrid_flag = 1" for examples of how these 4 functions can be
implemented.

The *get_grid_by_name()* method takes a grid name as input and returns
two values.  The *dim* argument is returned as 2 or 3 for the
dimensionality of the grid.  The function return is a grid index from
0 to G-1 where *G* is the number of grids the command instantiates.  A
value of -1 is returned if the grid name is not recognized.

The *get_grid_by_index()* method is called after the
*get_grid_by_name()* method, using the grid index it returned as its
argument.  This method will return a pointer to the Grid2d or Grid3d
class.  The caller can use this to query grid attributes, such as the
global size of the grid.  The :doc:`dump grid <dump>` to insure each
its grid reference arguments are for grids of the same size.

The *get_griddata_by_name()* method takes a grid index *igrid* and a
data name as input.  It returns two values.  The *ncol* argument is
returned as a 0 if the grid data is a single value (scalar) per grid
cell, or an integer M > 0 if there are M values (vector) per grid
cell.  Note that even if M = 1, it is still a 1-length vector, not a
scalar.  The function return is a data index from 0 to D-1 where *D*
is the number of data sets associated with that grid by the command.
A value of -1 is returned if the data name is not recognized.

The *get_griddata_by_index()* method is called after the
*get_griddata_by_name()* method, using the data index it returned as
its argument.  This method will return a pointer to the
multi-dimensional array which stores the requested data.

As in the discussion above of the Memory class *create_offset()*
methods, the dimensionality of the array associated with the returned
pointer depends on whether it is a 2d or 3d grid and whether there is
a single or multiple values stored for each grid cell:

* single value per cell for a 2d grid = 2d array pointer
* multiple values per cell for a 2d grid = 3d array pointer
* single value per cell for a 3d grid = 3d array pointer
* multiple values per cell for a 3d grid = 4d array pointer

The caller will typically access the data by casting the void pointer
to the corresponding array pointer and using nested loops in x,y,z
between owned or ghost index bounds returned by the
*get_bounds_owned()* or *get_bounds_ghost()* methods to index into the
array.  Example code snippets with this logic were listed above,

----------

Final notes
^^^^^^^^^^^

Finally, here are some additional issues to pay attention to for
writing any style command which uses distributed grids via the Grid2d
or Grid3d class.

The command destructor should delete all instances of the Grid class,
any buffers it allocated for forward/reverse or remap communication,
and any data arrays it allocated to store grid cell data.

If a command is intended to work for either 2d or 3d simulations, then
it should have logic to instantiate either 2d or 3d grids and their
associated data arrays, depending on the dimension of the simulation
box.  The :doc:`fix ave/grid <fix_ave_grid>` command is an example of
such a command.

When a command maps its particles to the grid and updates grid cell
values, it should check that it is not updating or accessing a grid
cell value outside the range of its owned+ghost cells, and generate an
error message if that is the case.  This could happen, for example, if
a particle has moved further than half the neighbor skin distance,
because the neighbor list update criterion are not adequate to prevent
it from happening.  See the src/KSPACE/pppm.cpp file and its
*particle_map()* method for an example of this kind of error check.
