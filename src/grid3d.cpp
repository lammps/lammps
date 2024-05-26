// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "grid3d.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "irregular.h"
#include "pair.h"
#include "kspace.h"
#include "fix.h"
#include "math_extra.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int DELTA = 16;

static constexpr int OFFSET = 16384;

/* ----------------------------------------------------------------------
   NOTES:
   if o indices for ghosts are < 0 or hi indices are >= N,
     then grid is treated as periodic in that dimension,
     comm is done across the periodic boundaries
   tiled implementations only work for RCB, not general tilings
     b/c RCB tree is used to find neighboring tiles
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   constructor to create a 3d distributed grid
   Grid3d assigns owned/ghost cells to each proc via setup_grid()
     it MUST be called after constructor
   gcomm = caller's communicator
   gnx,gny,gnz = global grid size
------------------------------------------------------------------------- */

Grid3d::Grid3d(LAMMPS *lmp, MPI_Comm gcomm, int gnx, int gny, int gnz) :
  Pointers(lmp), swap(nullptr), requests(nullptr), srequest(nullptr), rrequest(nullptr),
    sresponse(nullptr), rresponse(nullptr), send(nullptr), recv(nullptr), copy(nullptr),
    send_remap(nullptr), recv_remap(nullptr), overlap_procs(nullptr), xsplit(nullptr),
    ysplit(nullptr), zsplit(nullptr), grid2proc(nullptr), rcbinfo(nullptr), overlap_list(nullptr)

{
  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;
  nz = gnz;

  // default settings, can be overridden by set() methods
  // these affect assignment of owned and ghost cells

  maxdist = 0.0;
  stencil_grid_lo = stencil_grid_hi = 0;
  stencil_atom_lo = stencil_atom_hi = 0;
  shift_grid = 0.5;
  shift_atom_lo = shift_atom_hi = 0.0;
  zextra = 0;
  zfactor = 1.0;
}

/* ----------------------------------------------------------------------
   alternate constructor to create a 3d distributed grid
   caller assigns owned/ghost cells to each proc
     setup_grid() must NOT be called
   used by MSM and PPPM/Electrode b/c their definition of ghost cells is complex
   gcomm = caller's communicator
   gnx,gny,gnz = global grid size
   i xyz lo/hi = extent of owned grid cells on this proc
   o xyz lo/hi = extent of owned+ghost grid cells on this proc
   owned and ghost indices are inclusive
     owned indices range from 0 to N-1
     ghost indices can extend < 0 or >= N
------------------------------------------------------------------------- */

Grid3d::Grid3d(LAMMPS *lmp, MPI_Comm gcomm, int gnx, int gny, int gnz,
               int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
               int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi) :
  Pointers(lmp), swap(nullptr), requests(nullptr), srequest(nullptr), rrequest(nullptr),
    sresponse(nullptr), rresponse(nullptr), send(nullptr), recv(nullptr), copy(nullptr),
    send_remap(nullptr), recv_remap(nullptr), overlap_procs(nullptr), xsplit(nullptr),
    ysplit(nullptr), zsplit(nullptr), grid2proc(nullptr), rcbinfo(nullptr), overlap_list(nullptr)
{
  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;
  nz = gnz;

  // store owned/ghost indices provided by caller

  inxlo = ixlo;
  inxhi = ixhi;
  inylo = iylo;
  inyhi = iyhi;
  inzlo = izlo;
  inzhi = izhi;

  outxlo = oxlo;
  outxhi = oxhi;
  outylo = oylo;
  outyhi = oyhi;
  outzlo = ozlo;
  outzhi = ozhi;

  // additional intialization
  // other constructor invokes this from setup_grid()

  Grid3d::initialize();
}

/* ---------------------------------------------------------------------- */

Grid3d::~Grid3d()
{
  // brick comm data structs

  for (int i = 0; i < nswap; i++) {
    memory->destroy(swap[i].packlist);
    memory->destroy(swap[i].unpacklist);
  }
  memory->sfree(swap);

  delete [] xsplit;
  delete [] ysplit;
  delete [] zsplit;
  memory->destroy(grid2proc);

  // tiled comm data structs

  for (int i = 0; i < nsend; i++)
    memory->destroy(send[i].packlist);
  memory->sfree(send);

  for (int i = 0; i < nrecv; i++)
    memory->destroy(recv[i].unpacklist);
  memory->sfree(recv);

  for (int i = 0; i < ncopy; i++) {
    memory->destroy(copy[i].packlist);
    memory->destroy(copy[i].unpacklist);
  }
  memory->sfree(copy);

  delete [] requests;
  delete [] requests_remap;

  memory->sfree(rcbinfo);

  // remap data structs

  deallocate_remap();
}

// ----------------------------------------------------------------------
// set Grid parameters
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   maxdist = max distance outside proc subdomain a particle can be
   used to determine extent of ghost cells
------------------------------------------------------------------------- */

void Grid3d::set_distance(double distance)
{
  maxdist = distance;
}

/* ----------------------------------------------------------------------
   # of grid cells beyond an owned grid cell that caller accesses
   used by FixTTMGrid for a finite different stencil
   can be different in lo vs hi direction
------------------------------------------------------------------------- */

void Grid3d::set_stencil_grid(int lo, int hi)
{
  stencil_grid_lo = lo;
  stencil_grid_hi = hi;
}

/* ----------------------------------------------------------------------
   # of grid cells beyond a particle's grid cell that caller accesses
   used by PPPM for smearing a point charge to the grid
   can be different in lo vs hi direction
------------------------------------------------------------------------- */

void Grid3d::set_stencil_atom(int lo, int hi)
{
  stencil_atom_lo = lo;
  stencil_atom_hi = hi;
}

/* ----------------------------------------------------------------------
   shift_grid = offset within grid cell of position of grid point
   0.5 = cell center (default), 0.0 = lower-left corner of cell
   used to decide which proc owns a grid cell (grid pt within subdomain)
------------------------------------------------------------------------- */

void Grid3d::set_shift_grid(double shift)
{
  shift_grid = shift;
}

/* ----------------------------------------------------------------------
   shift_atom = offset added to atoms when caller maps them to grid cells
   0.5 = half a grid cell, 0.0 = no offset
   used to compute maximum possible ghost extents
   use of lo/hi allows max ghost extent on each side to be different
   PPPM uses 0.5 when stencil order is odd, 0.0 when order is even
   PPPM/stagger applies different shift values for 2 stagger iterations
------------------------------------------------------------------------- */

void Grid3d::set_shift_atom(double shift_lo, double shift_hi)
{
  shift_atom_lo = shift_lo;
  shift_atom_hi = shift_hi;
}

/* ----------------------------------------------------------------------
   enable extra grid cells in Z
   factor = muliplication factor on box size Z and thus grid size
   factor > 1.0 when grid extends beyond Z box size (3.0 = tripled in size)
   used by PPPM for 2d periodic slab geometries
   only enable zextra if factor > 1.0
   default zextra = 0, factor = 1.0 (no extra grid cells in Z)
------------------------------------------------------------------------- */

void Grid3d::set_zfactor(double factor)
{
  if (factor == 1.0) zextra = 0;
  else zextra = 1;
  zfactor = factor;
}

/* ----------------------------------------------------------------------
   set IDs of proc neighbors used in uniform local owned/ghost comm
   must be called BEFORE setup_comm() to override
     the processor neighbors stored by extract_comm_info()
   used by MSM to exclude non-participating procs for coarse grid comm
------------------------------------------------------------------------- */

void Grid3d::set_proc_neighs(int pxlo, int pxhi, int pylo, int pyhi,
                             int pzlo, int pzhi)
{
  procxlo = pxlo;
  procxhi = pxhi;
  procylo = pylo;
  procyhi = pyhi;
  proczlo = pzlo;
  proczhi = pzhi;
}

/* ----------------------------------------------------------------------
   set allocation dimensions of caller grid used by indices() to setup pack/unpack
   must be called BEFORE setup_comm() to override
     the caller grid size used in indices()
   used by MSM to allow a larger level 0 grid to be allocated
     with more ghost cells for other operations
------------------------------------------------------------------------- */

void Grid3d::set_caller_grid(int fxlo, int fxhi, int fylo, int fyhi,
                             int fzlo, int fzhi)
{
  fullxlo = fxlo;
  fullxhi = fxhi;
  fullylo = fylo;
  fullyhi = fyhi;
  fullzlo = fzlo;
  fullzhi = fzhi;
}

// ----------------------------------------------------------------------
// retrieve Grid parameters
// ----------------------------------------------------------------------

int Grid3d::identical(Grid3d *grid2)
{
  int inxlo2,inxhi2,inylo2,inyhi2,inzlo2,inzhi2;
  int outxlo2,outxhi2,outylo2,outyhi2,outzlo2,outzhi2;

  grid2->get_bounds_owned(inxlo2,inxhi2,inylo2,inyhi2,inzlo2,inzhi2);
  grid2->get_bounds_ghost(outxlo2,outxhi2,outylo2,outyhi2,outzlo2,outzhi2);

  int flag = 0;
  if (inxlo != inxlo2 || inxhi != inxhi2 ||
      inylo != inylo2 || inyhi != inyhi2 ||
      inzlo != inzlo2 || inzhi != inzhi2) flag = 1;
  if (outxlo != outxlo2 || outxhi != outxhi2 ||
      outylo != outylo2 || outyhi != outyhi2 ||
      outzlo != outzlo2 || outzhi != outzhi2) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,gridcomm);

  if (flagall) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

void Grid3d::get_size(int &nxgrid, int &nygrid, int &nzgrid)
{
  nxgrid = nx;
  nygrid = ny;
  nzgrid = nz;
}

/* ---------------------------------------------------------------------- */

void Grid3d::get_bounds_owned(int &xlo, int &xhi, int &ylo, int &yhi,
                              int &zlo, int &zhi)
{
  xlo = inxlo;
  xhi = inxhi;
  ylo = inylo;
  yhi = inyhi;
  zlo = inzlo;
  zhi = inzhi;
}

/* ---------------------------------------------------------------------- */

void Grid3d::get_bounds_ghost(int &xlo, int &xhi, int &ylo, int &yhi,
                              int &zlo, int &zhi)
{
  xlo = outxlo;
  xhi = outxhi;
  ylo = outylo;
  yhi = outyhi;
  zlo = outzlo;
  zhi = outzhi;
}

// ----------------------------------------------------------------------
// define owned and ghost grid cells
// also store comm and grid partitioning info
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   setup grid partition for each proc = owned + ghost cells
   return:
     i xyz lohi = portion of global grid this proc owns, 0 <= index < N
     o xyz lohi = owned + ghost grid cells in all directions
   for periodic dims, o indices can be < 0 or >= N
   for non-periodic dims, o indices will be >= 0 and < N
     since no grid comm is done across non-periodic boundaries
------------------------------------------------------------------------- */

void Grid3d::setup_grid(int &ixlo, int &ixhi, int &iylo, int &iyhi,
                        int &izlo, int &izhi,
                        int &oxlo, int &oxhi, int &oylo, int &oyhi,
                        int &ozlo, int &ozhi)
{
  // owned grid cells = those whose grid point is within proc subdomain
  // shift_grid = 0.5 for grid point at cell center, 0.0 for lower-left corner

  double fraclo,frachi;

  if (comm->layout != Comm::LAYOUT_TILED) {
    fraclo = comm->xsplit[comm->myloc[0]];
    frachi = comm->xsplit[comm->myloc[0]+1];
    partition_grid(nx,fraclo,frachi,shift_grid,0,inxlo,inxhi);
    fraclo = comm->ysplit[comm->myloc[1]];
    frachi = comm->ysplit[comm->myloc[1]+1];
    partition_grid(ny,fraclo,frachi,shift_grid,0,inylo,inyhi);
    fraclo = comm->zsplit[comm->myloc[2]];
    frachi = comm->zsplit[comm->myloc[2]+1];
    partition_grid(nz,fraclo,frachi,shift_grid,zextra,inzlo,inzhi);
  } else {
    fraclo = comm->mysplit[0][0];
    frachi = comm->mysplit[0][1];
    partition_grid(nx,fraclo,frachi,shift_grid,0,inxlo,inxhi);
    fraclo = comm->mysplit[1][0];
    frachi = comm->mysplit[1][1];
    partition_grid(ny,fraclo,frachi,shift_grid,0,inylo,inyhi);
    fraclo = comm->mysplit[2][0];
    frachi = comm->mysplit[2][1];
    partition_grid(nz,fraclo,frachi,shift_grid,zextra,inzlo,inzhi);
  }

  // extend owned grid bounds with ghost grid cells in each direction

  ghost_grid();

  // additional intialization
  // other constructor invokes this directly

  initialize();

  // return values

  ixlo = inxlo;
  ixhi = inxhi;
  iylo = inylo;
  iyhi = inyhi;
  izlo = inzlo;
  izhi = inzhi;

  oxlo = outxlo;
  oxhi = outxhi;
  oylo = outylo;
  oyhi = outyhi;
  ozlo = outzlo;
  ozhi = outzhi;
}

/* ----------------------------------------------------------------------
   additional one-time setup common to both constructors
 ---------------------------------------------------------------------- */

void Grid3d::initialize()
{
  // error check on size of grid stored by this proc

  bigint total = (bigint)
    (outxhi - outxlo + 1) * (outyhi - outylo + 1) * (outzhi - outzlo + 1);
  if (total > MAXSMALLINT)
    error->one(FLERR, "Too many owned+ghost grid3d points");

  // default = caller grid is allocated to ghost grid
  // used when computing pack/unpack lists in indices()
  // these values can be overridden using set_caller_grid()

  fullxlo = outxlo;
  fullxhi = outxhi;
  fullylo = outylo;
  fullyhi = outyhi;
  fullzlo = outzlo;
  fullzhi = outzhi;

  // initialize data structs

  nswap = maxswap = 0;
  swap = nullptr;

  nsend = nrecv = ncopy = 0;
  send = nullptr;
  recv = nullptr;
  copy = nullptr;
  requests = nullptr;
  requests_remap = nullptr;

  xsplit = ysplit = zsplit = nullptr;
  grid2proc = nullptr;
  rcbinfo = nullptr;

  nsend_remap = nrecv_remap = self_remap = 0;
  send_remap = nullptr;
  recv_remap = nullptr;

  // store info about Comm decomposition needed for remap operation
  // two Grid instances will exist for duration of remap
  // each must know Comm decomp at time Grid instance was created

  extract_comm_info();
}

/* ----------------------------------------------------------------------
   partition a global regular grid into one brick-shaped sub-grid per proc
   if grid point is inside my sub-domain I own it,
     this includes sub-domain lo boundary but excludes hi boundary
   ngrid = extent of global grid in a dimension
     indices into the global grid range from 0 to Ngrid-1 in that dim
   shift determines position of grid pt within grid cell
     shift = 0.5 for cell center, 0.0 for lower-left corner
   extra = 0 if grid exactly covers the simulation box
   extra = 1 if grid extends beyond the +z boundary by zfactor (PPPM slab)
     effectively maps proc partitions to the box-size subset of the grid
   lo/hi = inclusive lo/hi bounds for brick of global grid cells I own
     lo grid index = first grid pt >= fraclo*bound
     hi grid index = last grid pt < frachi*bound
     if proc owns no grid cells in a dim, then inlo > inhi
   special case: 2 procs share boundary which a grid point is exactly on
     2 if test equalties ensure a consistent decision as to which proc owns it
------------------------------------------------------------------------- */

void Grid3d::partition_grid(int ngrid, double fraclo, double frachi,
                            double shift, int extra, int &lo, int &hi)
{
  if (extra == 0) {
    lo = static_cast<int> (fraclo * ngrid);
    while (lo+shift < fraclo*ngrid) lo++;
    hi = static_cast<int> (frachi * ngrid);
    while (hi+shift >= frachi*ngrid) hi--;
  } else {
    lo = static_cast<int> (fraclo * ngrid/zfactor);
    while (lo+shift < fraclo*ngrid/zfactor) lo++;
    hi = static_cast<int> (frachi * ngrid/zfactor);
    while (hi+shift >= frachi*ngrid/zfactor) hi--;
  }
}

/* ----------------------------------------------------------------------
   extend ghost grid cells in each direction beyond owned grid
   indices into the global grid range from 0 to N-1 in each dim
   ghost cell indices for periodic systems can be < 0 or >= N
------------------------------------------------------------------------- */

void Grid3d::ghost_grid()
{
  double *prd,*boxlo,*sublo,*subhi;
  int triclinic = domain->triclinic;

  if (triclinic == 0) {
    prd = domain->prd;
    boxlo = domain->boxlo;
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    prd = domain->prd_lamda;
    boxlo = domain->boxlo_lamda;
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // for triclinic, maxdist = different value in each orthogonal direction

  double dist[3] = {0.0,0.0,0.0};
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = maxdist;
  else MathExtra::tribbox(domain->h,maxdist,&dist[0]);

  // lo/hi = min/max index of global grid cells my owned atoms can be mapped to
  //   includes effects of maxdist and shift_atom settings
  // lo/hi can be further extended by stencil_atom and stencil_grid settings
  // all those settings are set by caller
  // ghost cell layers needed in each dim/dir = max of two extension effects
  // OFFSET allows generation of negative indices with static_cast
  // out xyz lo/hi = index range of owned + ghost cells
  // if zextra, nz and effective prd[2] are both larger, so dzinv is the same

  double dxinv = nx / prd[0];
  double dyinv = ny / prd[1];
  double dzinv = nz / prd[2];
  if (zextra) dzinv = nz / (prd[2] * zfactor);

  int lo, hi;

  lo = static_cast<int>((sublo[0]-dist[0]-boxlo[0]) * dxinv + shift_atom_lo + OFFSET) - OFFSET;
  hi = static_cast<int>((subhi[0]+dist[0]-boxlo[0]) * dxinv + shift_atom_hi + OFFSET) - OFFSET;
  outxlo = MIN(lo - stencil_atom_lo, inxlo - stencil_grid_lo);
  outxhi = MAX(hi + stencil_atom_hi, inxhi + stencil_grid_hi);

  lo = static_cast<int>((sublo[1]-dist[1]-boxlo[1]) * dyinv + shift_atom_lo + OFFSET) - OFFSET;
  hi = static_cast<int>((subhi[1]+dist[1]-boxlo[1]) * dyinv + shift_atom_hi + OFFSET) - OFFSET;
  outylo = MIN(lo - stencil_atom_lo, inylo - stencil_grid_lo);
  outyhi = MAX(hi + stencil_atom_hi, inyhi + stencil_grid_hi);

  lo = static_cast<int>((sublo[2]-dist[2]-boxlo[2]) * dzinv + shift_atom_lo + OFFSET) - OFFSET;
  hi = static_cast<int>((subhi[2]+dist[2]-boxlo[2]) * dzinv + shift_atom_hi + OFFSET) - OFFSET;
  outzlo = MIN(lo - stencil_atom_lo, inzlo - stencil_grid_lo);
  outzhi = MAX(hi + stencil_atom_hi, inzhi + stencil_grid_hi);

  // if zextra = 1:
  //   adjust grid boundaries for processors at +z end,
  //   to include added empty grid cells between periodically repeating slabs
  // in this case:
  //   want grid data forward communicated from +z proc to -z proc, but not vice versa
  //   want grid data reverse communicated from -z proc to +z proc, but not vice versa
  // this is accomplished by inzhi = outzhi on +z end (no ghost cells)
  // also ensure no other procs use ghost cells beyond +z limit

  if (zextra) {
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[2] == comm->procgrid[2]-1) inzhi = outzhi = nz - 1;
    } else {
      if (comm->mysplit[2][1] == 1.0) inzhi = outzhi = nz - 1;
    }
    outzhi = MIN(outzhi,nz-1);
  }

  // limit out xyz lo/hi indices to global grid for non-periodic dims
  // if zextra = 1 (e.g. PPPM), treat z-dimension as if periodic

  int *periodicity = domain->periodicity;

  if (!periodicity[0]) {
    outxlo = MAX(0,outxlo);
    outxhi = MIN(nx-1,outxhi);
  }

  if (!periodicity[1]) {
    outylo = MAX(0,outylo);
    outyhi = MIN(ny-1,outyhi);
  }

  if (!periodicity[2] && !zextra) {
    outzlo = MAX(0,outzlo);
    outzhi = MIN(nz-1,outzhi);
  }
}

/* ----------------------------------------------------------------------
   store copy of info from Comm class about processor partitioning
   used when a remap is performed between two Grid instances, old and new
------------------------------------------------------------------------- */

void Grid3d::extract_comm_info()
{
  // for non TILED layout:
  // proc xyz lohi = my 6 neighbor procs in this MPI_Comm
  //   these proc IDs can be overridden by caller using set_proc_neighs()
  // xyz split = copy of 1d vectors in Comm
  // grid2proc = copy of 3d array in Comm

  if (comm->layout != Comm::LAYOUT_TILED) {
    procxlo = comm->procneigh[0][0];
    procxhi = comm->procneigh[0][1];
    procylo = comm->procneigh[1][0];
    procyhi = comm->procneigh[1][1];
    proczlo = comm->procneigh[2][0];
    proczhi = comm->procneigh[2][1];

    xsplit = new double[comm->procgrid[0]+1];
    ysplit = new double[comm->procgrid[1]+1];
    zsplit = new double[comm->procgrid[2]+1];
    memcpy(xsplit, comm->xsplit, sizeof(double) * (comm->procgrid[0]+1));
    memcpy(ysplit, comm->ysplit, sizeof(double) * (comm->procgrid[1]+1));
    memcpy(zsplit, comm->zsplit, sizeof(double) * (comm->procgrid[2]+1));

    memory->create(grid2proc,comm->procgrid[0],comm->procgrid[1],comm->procgrid[2],
                   "grid3d:grid2proc");
    memcpy(&grid2proc[0][0][0],&comm->grid2proc[0][0][0],
           sizeof(int) * comm->procgrid[0] * comm->procgrid[1] * comm->procgrid[2]);
  }

  // for TILED layout:
  // create RCB tree of grid partitioning info for grid decomp
  // Comm provides dim info for this proc, stored as RCBinfo.dim
  //   dim is -1 for proc 0, but never accessed
  // RCBinfo.cut = this proc's inlo in that dim
  // Allgather creates the tree of dims and cuts

  if (comm->layout == Comm::LAYOUT_TILED) {
    rcbinfo = (RCBinfo *)
      memory->smalloc(nprocs*sizeof(RCBinfo),"grid3d:rcbinfo");
    RCBinfo rcbone;
    rcbone.dim = comm->rcbcutdim;
    if (rcbone.dim <= 0) rcbone.cut = inxlo;
    else if (rcbone.dim == 1) rcbone.cut = inylo;
    else if (rcbone.dim == 2) rcbone.cut = inzlo;
    MPI_Allgather(&rcbone,sizeof(RCBinfo),MPI_CHAR,
                  rcbinfo,sizeof(RCBinfo),MPI_CHAR,gridcomm);
  }
}

// ----------------------------------------------------------------------
// setup of local owned/ghost grid comm
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   setup commmunication of owned/ghost grid cells
     either for brick decomp or tiled decomp
   return sizes of two buffers needed for communication
     nbuf1 = largest pack or unpack in any Send or Recv or Copy
     nbuf2 = larget of sum of all packs or unpacks in Send or Recv
   for brick comm, nbuf1 = nbuf2
   for tiling comm, nbuf2 >= nbuf2
   nbuf1,nbuf2 are counts of grid cells
     caller converts them to message sizes for grid data it stores
------------------------------------------------------------------------- */

void Grid3d::setup_comm(int &nbuf1, int &nbuf2)
{
  if (comm->layout != Comm::LAYOUT_TILED) setup_comm_brick(nbuf1,nbuf2);
  else setup_comm_tiled(nbuf1,nbuf2);
}

/* ----------------------------------------------------------------------
   setup owned/ghost comm for brick comm
   each proc has 6 neighbors
   comm pattern = series of swaps with one of those 6 procs
   can be multiple swaps with same proc if ghost extent is large
   swap may not be symmetric if both procs do not need same layers of ghosts
   all procs perform same # of swaps in a direction, even if some don't need it
------------------------------------------------------------------------- */

void Grid3d::setup_comm_brick(int &nbuf1, int &nbuf2)
{
  int nsent,sendfirst,sendlast,recvfirst,recvlast;
  int sendplanes,recvplanes;
  int notdoneme,notdone;

  // notify 6 neighbor procs how many ghost grid planes I need from them
  // ghost xyz lo = # of my lower grid planes that proc xyz lo needs as its ghosts
  // ghost xyz hi = # of my upper grid planes that proc xyz hi needs as its ghosts
  // if this proc is its own neighbor across periodic bounary, value is from self

  int nplanes = inxlo - outxlo;
  if (procxlo != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,procxlo,0,
                   &ghostxhi,1,MPI_INT,procxhi,0,gridcomm,MPI_STATUS_IGNORE);
  else ghostxhi = nplanes;

  nplanes = outxhi - inxhi;
  if (procxhi != me)
      MPI_Sendrecv(&nplanes,1,MPI_INT,procxhi,0,
                   &ghostxlo,1,MPI_INT,procxlo,0,gridcomm,MPI_STATUS_IGNORE);
  else ghostxlo = nplanes;

  nplanes = inylo - outylo;
  if (procylo != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,procylo,0,
                 &ghostyhi,1,MPI_INT,procyhi,0,gridcomm,MPI_STATUS_IGNORE);
  else ghostyhi = nplanes;

  nplanes = outyhi - inyhi;
  if (procyhi != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,procyhi,0,
                 &ghostylo,1,MPI_INT,procylo,0,gridcomm,MPI_STATUS_IGNORE);
  else ghostylo = nplanes;

  nplanes = inzlo - outzlo;
  if (proczlo != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,proczlo,0,
                 &ghostzhi,1,MPI_INT,proczhi,0,gridcomm,MPI_STATUS_IGNORE);
  else ghostzhi = nplanes;

  nplanes = outzhi - inzhi;
  if (proczhi != me)
    MPI_Sendrecv(&nplanes,1,MPI_INT,proczhi,0,
                 &ghostzlo,1,MPI_INT,proczlo,0,gridcomm,MPI_STATUS_IGNORE);
  else ghostzlo = nplanes;

  // setup swaps = exchange of grid data with one of 6 neighobr procs
  // can be more than one in a direction if ghost region extends beyond neigh proc
  // all procs have same swap count, but swapsize npack/nunpack can be empty

  nswap = 0;

  // send own grid pts to -x processor, recv ghost grid pts from +x processor

  nsent = 0;
  sendfirst = inxlo;
  sendlast = inxhi;
  recvfirst = inxhi+1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) grow_swap();

    swap[nswap].sendproc = procxlo;
    swap[nswap].recvproc = procxhi;
    sendplanes = MIN(sendlast-sendfirst+1,ghostxlo-nsent);
    swap[nswap].npack =
      indices(swap[nswap].packlist,
              sendfirst,sendfirst+sendplanes-1,inylo,inyhi,inzlo,inzhi);

    if (procxlo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procxlo,0,
                   &recvplanes,1,MPI_INT,procxhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              recvfirst,recvfirst+recvplanes-1,inylo,inyhi,inzlo,inzhi);

    nsent += sendplanes;
    sendfirst += sendplanes;
    sendlast += recvplanes;
    recvfirst += recvplanes;
    nswap++;

    if (nsent < ghostxlo) notdoneme = 1;
    else notdoneme = 0;
    MPI_Allreduce(&notdoneme,&notdone,1,MPI_INT,MPI_SUM,gridcomm);
  }

  // send own grid pts to +x processor, recv ghost grid pts from -x processor

  nsent = 0;
  sendfirst = inxlo;
  sendlast = inxhi;
  recvlast = inxlo-1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) grow_swap();

    swap[nswap].sendproc = procxhi;
    swap[nswap].recvproc = procxlo;
    sendplanes = MIN(sendlast-sendfirst+1,ghostxhi-nsent);
    swap[nswap].npack =
      indices(swap[nswap].packlist,
              sendlast-sendplanes+1,sendlast,inylo,inyhi,inzlo,inzhi);

    if (procxhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procxhi,0,
                   &recvplanes,1,MPI_INT,procxlo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              recvlast-recvplanes+1,recvlast,inylo,inyhi,inzlo,inzhi);

    nsent += sendplanes;
    sendfirst -= recvplanes;
    sendlast -= sendplanes;
    recvlast -= recvplanes;
    nswap++;

    if (nsent < ghostxhi) notdoneme = 1;
    else notdoneme = 0;
    MPI_Allreduce(&notdoneme,&notdone,1,MPI_INT,MPI_SUM,gridcomm);
  }

  // send own grid pts to -y processor, recv ghost grid pts from +y processor

  nsent = 0;
  sendfirst = inylo;
  sendlast = inyhi;
  recvfirst = inyhi+1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) grow_swap();

    swap[nswap].sendproc = procylo;
    swap[nswap].recvproc = procyhi;
    sendplanes = MIN(sendlast-sendfirst+1,ghostylo-nsent);
    swap[nswap].npack =
      indices(swap[nswap].packlist,
              outxlo,outxhi,sendfirst,sendfirst+sendplanes-1,inzlo,inzhi);

    if (procylo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procylo,0,
                   &recvplanes,1,MPI_INT,procyhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              outxlo,outxhi,recvfirst,recvfirst+recvplanes-1,inzlo,inzhi);

    nsent += sendplanes;
    sendfirst += sendplanes;
    sendlast += recvplanes;
    recvfirst += recvplanes;
    nswap++;

    if (nsent < ghostylo) notdoneme = 1;
    else notdoneme = 0;
    MPI_Allreduce(&notdoneme,&notdone,1,MPI_INT,MPI_SUM,gridcomm);
  }

  // send own grid pts to +y processor, recv ghost grid pts from -y processor

  nsent = 0;
  sendfirst = inylo;
  sendlast = inyhi;
  recvlast = inylo-1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) grow_swap();

    swap[nswap].sendproc = procyhi;
    swap[nswap].recvproc = procylo;
    sendplanes = MIN(sendlast-sendfirst+1,ghostyhi-nsent);
    swap[nswap].npack =
      indices(swap[nswap].packlist,
              outxlo,outxhi,sendlast-sendplanes+1,sendlast,inzlo,inzhi);

    if (procyhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procyhi,0,
                   &recvplanes,1,MPI_INT,procylo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              outxlo,outxhi,recvlast-recvplanes+1,recvlast,inzlo,inzhi);

    nsent += sendplanes;
    sendfirst -= recvplanes;
    sendlast -= sendplanes;
    recvlast -= recvplanes;
    nswap++;

    if (nsent < ghostyhi) notdoneme = 1;
    else notdoneme = 0;
    MPI_Allreduce(&notdoneme,&notdone,1,MPI_INT,MPI_SUM,gridcomm);
  }

  // send own grid pts to -z processor, recv ghost grid pts from +z processor

  nsent = 0;
  sendfirst = inzlo;
  sendlast = inzhi;
  recvfirst = inzhi+1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) grow_swap();

    swap[nswap].sendproc = proczlo;
    swap[nswap].recvproc = proczhi;
    sendplanes = MIN(sendlast-sendfirst+1,ghostzlo-nsent);
    swap[nswap].npack =
      indices(swap[nswap].packlist,
              outxlo,outxhi,outylo,outyhi,sendfirst,sendfirst+sendplanes-1);

    if (proczlo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,proczlo,0,
                   &recvplanes,1,MPI_INT,proczhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              outxlo,outxhi,outylo,outyhi,recvfirst,recvfirst+recvplanes-1);

    nsent += sendplanes;
    sendfirst += sendplanes;
    sendlast += recvplanes;
    recvfirst += recvplanes;
    nswap++;

    if (nsent < ghostzlo) notdoneme = 1;
    else notdoneme = 0;
    MPI_Allreduce(&notdoneme,&notdone,1,MPI_INT,MPI_SUM,gridcomm);
  }

  // send own grid pts to +z processor, recv ghost grid pts from -z processor

  nsent = 0;
  sendfirst = inzlo;
  sendlast = inzhi;
  recvlast = inzlo-1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) grow_swap();

    swap[nswap].sendproc = proczhi;
    swap[nswap].recvproc = proczlo;
    sendplanes = MIN(sendlast-sendfirst+1,ghostzhi-nsent);
    swap[nswap].npack =
      indices(swap[nswap].packlist,
              outxlo,outxhi,outylo,outyhi,sendlast-sendplanes+1,sendlast);

    if (proczhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,proczhi,0,
                   &recvplanes,1,MPI_INT,proczlo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              outxlo,outxhi,outylo,outyhi,recvlast-recvplanes+1,recvlast);

    nsent += sendplanes;
    sendfirst -= recvplanes;
    sendlast -= sendplanes;
    recvlast -= recvplanes;
    nswap++;

    if (nsent < ghostzhi) notdoneme = 1;
    else notdoneme = 0;
    MPI_Allreduce(&notdoneme,&notdone,1,MPI_INT,MPI_SUM,gridcomm);
  }

  // ngrid = max of any forward/reverse pack/unpack grid points

  int ngrid = 0;
  for (int i = 0; i < nswap; i++) {
    ngrid = MAX(ngrid,swap[i].npack);
    ngrid = MAX(ngrid,swap[i].nunpack);
  }

  nbuf1 = nbuf2 = ngrid;
}

/* ----------------------------------------------------------------------
   setup owned/ghost comm for tiled comm
   each proc has arbitrary # of neighbors that overlap its ghost extent
   identify which procs will send me ghost cells, and vice versa
   may not be symmetric if both procs do not need same layers of ghosts
   comm pattern = post recvs for all my ghosts, send my owned, wait on recvs
     no exchanges by dimension, unlike CommTiled forward/reverse comm of particles
------------------------------------------------------------------------- */

void Grid3d::setup_comm_tiled(int &nbuf1, int &nbuf2)
{
  int i,m;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int ghostbox[6],pbc[3];

  // find overlaps of my extended ghost box with all owned boxes
  // accounts for ghost box overlapping periodic boundaries
  // noverlap = # of overlaps, including self
  // overlap = vector of overlap info using Overlap data struct

  ghostbox[0] = outxlo;
  ghostbox[1] = outxhi;
  ghostbox[2] = outylo;
  ghostbox[3] = outyhi;
  ghostbox[4] = outzlo;
  ghostbox[5] = outzhi;

  pbc[0] = pbc[1] = pbc[2] = 0;

  Overlap *overlap;
  int noverlap = compute_overlap(1,ghostbox,pbc,overlap);

  // send each proc an overlap message
  // content: me, index of my overlap, box that overlaps with its owned cells
  // ncopy = # of overlaps with myself across a periodic boundary
  //         skip copy to self when non-PBC

  int *proclist;
  memory->create(proclist,noverlap,"grid3d:proclist");
  srequest = (Request *)
    memory->smalloc(noverlap*sizeof(Request),"grid3d:srequest");

  int nsend_request = 0;
  ncopy = 0;

  for (m = 0; m < noverlap; m++) {
    if (overlap[m].proc == me) {
      if (overlap[m].pbc[0] == 0 && overlap[m].pbc[1] == 0 &&
          overlap[m].pbc[2] == 0) continue;
      ncopy++;
    } else {
      proclist[nsend_request] = overlap[m].proc;
      srequest[nsend_request].sender = me;
      srequest[nsend_request].index = m;
      for (i = 0; i < 6; i++)
        srequest[nsend_request].box[i] = overlap[m].box[i];
      nsend_request++;
    }
  }

  auto irregular = new Irregular(lmp);
  int nrecv_request = irregular->create_data(nsend_request,proclist,1);
  auto rrequest = (Request *) memory->smalloc(nrecv_request*sizeof(Request),"grid3d:rrequest");
  irregular->exchange_data((char *) srequest,sizeof(Request),(char *) rrequest);
  irregular->destroy_data();

  // compute overlaps between received ghost boxes and my owned box
  // overlap box used to setup my Send data struct and respond to requests

  send = (Send *) memory->smalloc(nrecv_request*sizeof(Send),"grid3d:send");
  sresponse = (Response *) memory->smalloc(nrecv_request*sizeof(Response),"grid3d:sresponse");
  memory->destroy(proclist);
  memory->create(proclist,nrecv_request,"grid3d:proclist");

  for (m = 0; m < nrecv_request; m++) {
    send[m].proc = rrequest[m].sender;
    xlo = MAX(rrequest[m].box[0],inxlo);
    xhi = MIN(rrequest[m].box[1],inxhi);
    ylo = MAX(rrequest[m].box[2],inylo);
    yhi = MIN(rrequest[m].box[3],inyhi);
    zlo = MAX(rrequest[m].box[4],inzlo);
    zhi = MIN(rrequest[m].box[5],inzhi);
    send[m].npack = indices(send[m].packlist,xlo,xhi,ylo,yhi,zlo,zhi);

    proclist[m] = rrequest[m].sender;
    sresponse[m].index = rrequest[m].index;
    sresponse[m].box[0] = xlo;
    sresponse[m].box[1] = xhi;
    sresponse[m].box[2] = ylo;
    sresponse[m].box[3] = yhi;
    sresponse[m].box[4] = zlo;
    sresponse[m].box[5] = zhi;
  }

  nsend = nrecv_request;

  // reply to each Request message with a Response message
  // content: index for the overlap on requestor, overlap box on my owned grid

  int nsend_response = nrecv_request;
  int nrecv_response = irregular->create_data(nsend_response,proclist,1);
  auto rresponse = (Response *) memory->smalloc(nrecv_response*sizeof(Response),"grid3d:rresponse");
  irregular->exchange_data((char *) sresponse,sizeof(Response),(char *) rresponse);
  irregular->destroy_data();
  delete irregular;

  // process received responses
  // box used to setup my Recv data struct after unwrapping via PBC
  // adjacent = 0 if any box of ghost cells does not adjoin my owned cells

  recv = (Recv *) memory->smalloc(nrecv_response*sizeof(Recv),"grid3d:recv");
  adjacent = 1;

  for (i = 0; i < nrecv_response; i++) {
    m = rresponse[i].index;
    recv[i].proc = overlap[m].proc;
    xlo = rresponse[i].box[0] + overlap[m].pbc[0] * nx;
    xhi = rresponse[i].box[1] + overlap[m].pbc[0] * nx;
    ylo = rresponse[i].box[2] + overlap[m].pbc[1] * ny;
    yhi = rresponse[i].box[3] + overlap[m].pbc[1] * ny;
    zlo = rresponse[i].box[4] + overlap[m].pbc[2] * nz;
    zhi = rresponse[i].box[5] + overlap[m].pbc[2] * nz;
    recv[i].nunpack = indices(recv[i].unpacklist,xlo,xhi,ylo,yhi,zlo,zhi);

    if (xlo != inxhi+1 && xhi != inxlo-1 &&
        ylo != inyhi+1 && yhi != inylo-1 &&
        zlo != inzhi+1 && zhi != inzlo-1) adjacent = 0;
  }

  nrecv = nrecv_response;

  // create Copy data struct from overlaps with self
  // skip copy to self when non-PBC

  copy = (Copy *) memory->smalloc(ncopy*sizeof(Copy),"grid3d:copy");

  ncopy = 0;
  for (m = 0; m < noverlap; m++) {
    if (overlap[m].proc != me) continue;
    if (overlap[m].pbc[0] == 0 && overlap[m].pbc[1] == 0 &&
        overlap[m].pbc[2] == 0) continue;
    xlo = overlap[m].box[0];
    xhi = overlap[m].box[1];
    ylo = overlap[m].box[2];
    yhi = overlap[m].box[3];
    zlo = overlap[m].box[4];
    zhi = overlap[m].box[5];
    copy[ncopy].npack = indices(copy[ncopy].packlist,xlo,xhi,ylo,yhi,zlo,zhi);
    xlo = overlap[m].box[0] + overlap[m].pbc[0] * nx;
    xhi = overlap[m].box[1] + overlap[m].pbc[0] * nx;
    ylo = overlap[m].box[2] + overlap[m].pbc[1] * ny;
    yhi = overlap[m].box[3] + overlap[m].pbc[1] * ny;
    zlo = overlap[m].box[4] + overlap[m].pbc[2] * nz;
    zhi = overlap[m].box[5] + overlap[m].pbc[2] * nz;
    copy[ncopy].nunpack = indices(copy[ncopy].unpacklist,xlo,xhi,ylo,yhi,zlo,zhi);
    ncopy++;
  }

  // set offsets for received data

  int offset = 0;
  for (m = 0; m < nsend; m++) {
    send[m].offset = offset;
    offset += send[m].npack;
  }

  offset = 0;
  for (m = 0; m < nrecv; m++) {
    recv[m].offset = offset;
    offset += recv[m].nunpack;
  }

  // length of MPI requests vector is max of nsend, nrecv

  int nrequest = MAX(nsend,nrecv);
  requests = new MPI_Request[nrequest];

  // clean-up

  clean_overlap();
  memory->destroy(proclist);
  memory->sfree(srequest);
  memory->sfree(rrequest);
  memory->sfree(sresponse);
  memory->sfree(rresponse);

  // nbuf1 = largest pack or unpack in any Send or Recv or Copy
  // nbuf2 = larget of sum of all packs or unpacks in Send or Recv

  nbuf1 = 0;

  for (m = 0; m < ncopy; m++) {
    nbuf1 = MAX(nbuf1,copy[m].npack);
    nbuf1 = MAX(nbuf1,copy[m].nunpack);
  }

  int nbufs = 0;
  for (m = 0; m < nsend; m++) {
    nbuf1 = MAX(nbuf1,send[m].npack);
    nbufs += send[m].npack;
  }

  int nbufr = 0;
  for (m = 0; m < nrecv; m++) {
    nbuf1 = MAX(nbuf1,recv[m].nunpack);
    nbufr += recv[m].nunpack;
  }

  nbuf2 = MAX(nbufs,nbufr);
}

// ----------------------------------------------------------------------
// query locality of forward/reverse grid comm
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   check if all procs only need ghost info from adjacent procs
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int Grid3d::ghost_adjacent()
{
  if (comm->layout != Comm::LAYOUT_TILED) return ghost_adjacent_brick();
  return ghost_adjacent_tiled();
}

/* ----------------------------------------------------------------------
   adjacent = 0 if a proc's ghost xyz lohi values exceed its subdomain size
   return 0 if adjacent=0 for any proc, else 1
------------------------------------------------------------------------- */

int Grid3d::ghost_adjacent_brick()
{
  adjacent = 1;
  if (ghostxlo > inxhi-inxlo+1) adjacent = 0;
  if (ghostxhi > inxhi-inxlo+1) adjacent = 0;
  if (ghostylo > inyhi-inylo+1) adjacent = 0;
  if (ghostyhi > inyhi-inylo+1) adjacent = 0;
  if (ghostzlo > inzhi-inzlo+1) adjacent = 0;
  if (ghostzhi > inzhi-inzlo+1) adjacent = 0;

  int adjacent_all;
  MPI_Allreduce(&adjacent,&adjacent_all,1,MPI_INT,MPI_MIN,gridcomm);
  return adjacent_all;
}

/* ----------------------------------------------------------------------
   adjacent = 0 if a proc's received ghosts were flagged
     as non-adjacent in setup_tiled()
   return 0 if adjacent=0 for any proc, else 1
------------------------------------------------------------------------- */

int Grid3d::ghost_adjacent_tiled()
{
  int adjacent_all;
  MPI_Allreduce(&adjacent,&adjacent_all,1,MPI_INT,MPI_MIN,gridcomm);
  return adjacent_all;
}

// ----------------------------------------------------------------------
// forward/reverse comm of owned/ghost grid data via callbacks
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   forward comm of my owned cells to other's ghost cells
------------------------------------------------------------------------- */

void Grid3d::forward_comm(int caller, void *ptr, int which, int nper, int nbyte,
                            void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (comm->layout != Comm::LAYOUT_TILED) {
    if (caller == KSPACE)
      forward_comm_brick<KSpace>((KSpace *) ptr,which,nper,nbyte,
                                 buf1,buf2,datatype);
    else if (caller == PAIR)
      forward_comm_brick<Pair>((Pair *) ptr,which,nper,nbyte,
                               buf1,buf2,datatype);
    else if (caller == FIX)
      forward_comm_brick<Fix>((Fix *) ptr,which,nper,nbyte,
                              buf1,buf2,datatype);
  } else {
    if (caller == KSPACE)
      forward_comm_tiled<KSpace>((KSpace *) ptr,which,nper,nbyte,
                                 buf1,buf2,datatype);
    else if (caller == PAIR)
      forward_comm_tiled<Pair>((Pair *) ptr,which,nper,nbyte,
                               buf1,buf2,datatype);
    else if (caller == FIX)
      forward_comm_tiled<Fix>((Fix *) ptr,which,nper,nbyte,
                              buf1,buf2,datatype);
  }
}

/* ----------------------------------------------------------------------
   forward comm for brick decomp via list of swaps with 6 neighbor procs
------------------------------------------------------------------------- */

template < class T >
void Grid3d::
forward_comm_brick(T *ptr, int which, int nper, int /*nbyte*/,
                   void *buf1, void *buf2, MPI_Datatype datatype)
{
  int m;
  MPI_Request request;

  for (m = 0; m < nswap; m++) {
    if (swap[m].sendproc == me)
      ptr->pack_forward_grid(which,buf2,swap[m].npack,swap[m].packlist);
    else
      ptr->pack_forward_grid(which,buf1,swap[m].npack,swap[m].packlist);

    if (swap[m].sendproc != me) {
      if (swap[m].nunpack) MPI_Irecv(buf2,nper*swap[m].nunpack,datatype,
                                     swap[m].recvproc,0,gridcomm,&request);
      if (swap[m].npack) MPI_Send(buf1,nper*swap[m].npack,datatype,
                                  swap[m].sendproc,0,gridcomm);
      if (swap[m].nunpack) MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    ptr->unpack_forward_grid(which,buf2,swap[m].nunpack,swap[m].unpacklist);
  }
}

/* ----------------------------------------------------------------------
   forward comm for tiled decomp via Send/Recv lists of each neighbor proc
------------------------------------------------------------------------- */

template < class T >
void Grid3d::
forward_comm_tiled(T *ptr, int which, int nper, int nbyte,
                   void *buf1, void *vbuf2, MPI_Datatype datatype)
{
  int i,m,offset;

  auto buf2 = (char *) vbuf2;

  // post all receives

  for (m = 0; m < nrecv; m++) {
    offset = nper * recv[m].offset * nbyte;
    MPI_Irecv((void *) &buf2[offset],nper*recv[m].nunpack,datatype,
              recv[m].proc,0,gridcomm,&requests[m]);
  }

  // perform all sends to other procs

  for (m = 0; m < nsend; m++) {
    ptr->pack_forward_grid(which,buf1,send[m].npack,send[m].packlist);
    MPI_Send(buf1,nper*send[m].npack,datatype,send[m].proc,0,gridcomm);
  }

  // perform all copies to self

  for (m = 0; m < ncopy; m++) {
    ptr->pack_forward_grid(which,buf1,copy[m].npack,copy[m].packlist);
    ptr->unpack_forward_grid(which,buf1,copy[m].nunpack,copy[m].unpacklist);
  }

  // unpack all received data

  for (i = 0; i < nrecv; i++) {
    MPI_Waitany(nrecv,requests,&m,MPI_STATUS_IGNORE);
    offset = nper * recv[m].offset * nbyte;
    ptr->unpack_forward_grid(which,(void *) &buf2[offset],
                             recv[m].nunpack,recv[m].unpacklist);
  }
}

/* ----------------------------------------------------------------------
   reverse comm of my ghost cells to sum to owner cells
------------------------------------------------------------------------- */

void Grid3d::reverse_comm(int caller, void *ptr, int which, int nper, int nbyte,
                            void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (comm->layout != Comm::LAYOUT_TILED) {
    if (caller == KSPACE)
      reverse_comm_brick<KSpace>((KSpace *) ptr,which,nper,nbyte,
                                buf1,buf2,datatype);
    else if (caller == PAIR)
      reverse_comm_brick<Pair>((Pair *) ptr,which,nper,nbyte,
                               buf1,buf2,datatype);
    else if (caller == FIX)
      reverse_comm_brick<Fix>((Fix *) ptr,which,nper,nbyte,
                              buf1,buf2,datatype);
  } else {
    if (caller == KSPACE)
      reverse_comm_tiled<KSpace>((KSpace *) ptr,which,nper,nbyte,
                                 buf1,buf2,datatype);
    else if (caller == PAIR)
      reverse_comm_tiled<Pair>((Pair *) ptr,which,nper,nbyte,
                               buf1,buf2,datatype);
    else if (caller == FIX)
      reverse_comm_tiled<Fix>((Fix *) ptr,which,nper,nbyte,
                              buf1,buf2,datatype);
  }
}

/* ----------------------------------------------------------------------
   reverse comm for brick decomp via list of swaps with 6 neighbor procs
------------------------------------------------------------------------- */

template < class T >
void Grid3d::
reverse_comm_brick(T *ptr, int which, int nper, int /*nbyte*/,
                   void *buf1, void *buf2, MPI_Datatype datatype)
{
  int m;
  MPI_Request request;

  for (m = nswap-1; m >= 0; m--) {
    if (swap[m].recvproc == me)
      ptr->pack_reverse_grid(which,buf2,swap[m].nunpack,swap[m].unpacklist);
    else
      ptr->pack_reverse_grid(which,buf1,swap[m].nunpack,swap[m].unpacklist);

    if (swap[m].recvproc != me) {
      if (swap[m].npack) MPI_Irecv(buf2,nper*swap[m].npack,datatype,
                                   swap[m].sendproc,0,gridcomm,&request);
      if (swap[m].nunpack) MPI_Send(buf1,nper*swap[m].nunpack,datatype,
                                     swap[m].recvproc,0,gridcomm);
      if (swap[m].npack) MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    ptr->unpack_reverse_grid(which,buf2,swap[m].npack,swap[m].packlist);
  }
}

/* ----------------------------------------------------------------------
   reverse comm for tiled decomp via Send/Recv lists of each neighbor proc
------------------------------------------------------------------------- */

template < class T >
void Grid3d::
reverse_comm_tiled(T *ptr, int which, int nper, int nbyte,
                   void *buf1, void *vbuf2, MPI_Datatype datatype)
{
  int i,m,offset;

  auto buf2 = (char *) vbuf2;

  // post all receives

  for (m = 0; m < nsend; m++) {
    offset = nper * send[m].offset * nbyte;
    MPI_Irecv((void *) &buf2[offset],nper*send[m].npack,datatype,
              send[m].proc,0,gridcomm,&requests[m]);
  }

  // perform all sends to other procs

  for (m = 0; m < nrecv; m++) {
    ptr->pack_reverse_grid(which,buf1,recv[m].nunpack,recv[m].unpacklist);
    MPI_Send(buf1,nper*recv[m].nunpack,datatype,recv[m].proc,0,gridcomm);
  }

  // perform all copies to self

  for (m = 0; m < ncopy; m++) {
    ptr->pack_reverse_grid(which,buf1,copy[m].nunpack,copy[m].unpacklist);
    ptr->unpack_reverse_grid(which,buf1,copy[m].npack,copy[m].packlist);
  }

  // unpack all received data

  for (i = 0; i < nsend; i++) {
    MPI_Waitany(nsend,requests,&m,MPI_STATUS_IGNORE);
    offset = nper * send[m].offset * nbyte;
    ptr->unpack_reverse_grid(which,(void *) &buf2[offset],
                             send[m].npack,send[m].packlist);
  }
}

// ----------------------------------------------------------------------
// remap comm between 2 old/new grid decomposition of owned grid data
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   setup remap from old grid decomposition to this grid decomposition
   return sizes of two buffers needed for communication
   either for brick decomp or tiling decomp
   nbuf1 = largest pack or unpack in any Send or Recv or Copy
   nbuf2 = larget of sum of all packs or unpacks in Send or Recv
   for brick comm, nbuf1 = nbuf2
   for tiled comm, nbuf2 >= nbuf2
   nbuf1,nbuf2 are just count of grid points
     caller converts them to message size for grid data it stores
------------------------------------------------------------------------- */

void Grid3d::setup_remap(Grid3d *old, int &nremap_buf1, int &nremap_buf2)
{
  int m;
  int pbc[3];
  int *box;

  // deallocated existing remap data structs

  deallocate_remap();

  // overlaps of my old decomp owned box with all owned boxes in new decomp
  // noverlap_old = # of overlaps, including self
  // overlap_old = vector of overlap info in Overlap data struct

  int oldbox[6];
  old->get_bounds_owned(oldbox[0],oldbox[1],oldbox[2],oldbox[3],
                        oldbox[4],oldbox[5]);
  pbc[0] = pbc[1] = pbc[2] = 0;

  Overlap *overlap_old;
  int noverlap_old = compute_overlap(0,oldbox,pbc,overlap_old);

  // use overlap_old to construct send and copy lists
  // skip overlaps that contain no grid cells

  self_remap = 0;

  nsend_remap = 0;
  for (m = 0; m < noverlap_old; m++) {
    box = overlap_old[m].box;
    if (box[0] > box[1] || box[2] > box[3] || box[4] > box[5]) continue;
    if (overlap_old[m].proc == me) self_remap = 1;
    else nsend_remap++;
  }

  send_remap = new Send[nsend_remap];

  nsend_remap = 0;
  for (m = 0; m < noverlap_old; m++) {
    box = overlap_old[m].box;
    if (box[0] > box[1] || box[2] > box[3] || box[4] > box[5]) continue;
    if (overlap_old[m].proc == me) {
      copy_remap.npack =
        old->indices(copy_remap.packlist,
                     box[0],box[1],box[2],box[3],box[4],box[5]);
    } else {
      send_remap[nsend_remap].proc = overlap_old[m].proc;
      send_remap[nsend_remap].npack =
        old->indices(send_remap[nsend_remap].packlist,
                     box[0],box[1],box[2],box[3],box[4],box[5]);
      nsend_remap++;
    }
  }

  // overlaps of my new decomp owned box with all owned boxes in old decomp
  // noverlap_new = # of overlaps, including self
  // overlap_new = vector of overlap info in Overlap data struct

  int newbox[6];
  get_bounds_owned(newbox[0],newbox[1],newbox[2],newbox[3],newbox[4],newbox[5]);
  pbc[0] = pbc[1] = pbc[2] = 0;

  Overlap *overlap_new;
  int noverlap_new = old->compute_overlap(0,newbox,pbc,overlap_new);

  // use overlap_new to construct recv and copy lists
  // skip overlaps that contain no grid cells
  // set offsets for Recv data

  nrecv_remap = 0;
  for (m = 0; m < noverlap_new; m++) {
    box = overlap_new[m].box;
    if (box[0] > box[1] || box[2] > box[3] || box[4] > box[5]) continue;
    if (overlap_new[m].proc != me) nrecv_remap++;
  }

  recv_remap = new Recv[nrecv_remap];

  nrecv_remap = 0;
  for (m = 0; m < noverlap_new; m++) {
    box = overlap_new[m].box;
    if (box[0] > box[1] || box[2] > box[3] || box[4] > box[5]) continue;
    if (overlap_new[m].proc == me) {
      copy_remap.nunpack =
        indices(copy_remap.unpacklist,
                box[0],box[1],box[2],box[3],box[4],box[5]);
    } else {
      recv_remap[nrecv_remap].proc = overlap_new[m].proc;
      recv_remap[nrecv_remap].nunpack =
        indices(recv_remap[nrecv_remap].unpacklist,
                box[0],box[1],box[2],box[3],box[4],box[5]);
      nrecv_remap++;
    }
  }

  // set offsets for received data

  int offset = 0;
  for (m = 0; m < nrecv_remap; m++) {
    recv_remap[m].offset = offset;
    offset += recv_remap[m].nunpack;
  }

  // length of MPI requests vector = nrecv_remap

  delete [] requests_remap;
  requests_remap = new MPI_Request[nrecv_remap];

  // clean-up

  clean_overlap();
  old->clean_overlap();

  // nremap_buf1 = largest pack or unpack in any Send or Recv or Copy
  // nremap_buf2 = sum of all unpacks in Recv

  nremap_buf1 = 0;

  if (self_remap) {
    nremap_buf1 = MAX(nremap_buf1,copy_remap.npack);
    nremap_buf1 = MAX(nremap_buf1,copy_remap.nunpack);
  }

  for (m = 0; m < nsend_remap; m++)
    nremap_buf1 = MAX(nremap_buf1,send_remap[m].npack);

  nremap_buf2 = 0;
  for (m = 0; m < nrecv_remap; m++) {
    nremap_buf1 = MAX(nremap_buf1,recv_remap[m].nunpack);
    nremap_buf2 += recv_remap[m].nunpack;
  }
}

/* ----------------------------------------------------------------------
   perform remap from old grid decomposition to this grid decomposition
   pack/unpack operations are performed by caller via callbacks
------------------------------------------------------------------------- */

void Grid3d::remap(int caller, void *ptr, int which, int nper, int nbyte,
                   void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (caller == FIX) remap_style<Fix>((Fix *) ptr,which,nper,nbyte,buf1,buf2,datatype);
}

/* ------------------------------------------------------------------------- */

template < class T >
void Grid3d::remap_style(T *ptr, int which, int nper, int nbyte,
                         void *buf1, void *vbuf2, MPI_Datatype datatype)
{
  int i,m,offset;

  auto buf2 = (char *) vbuf2;

  // post all receives

  for (m = 0; m < nrecv_remap; m++) {
    offset = nper * recv_remap[m].offset * nbyte;
    MPI_Irecv((void *) &buf2[offset],nper*recv_remap[m].nunpack,datatype,
              recv_remap[m].proc,0,gridcomm,&requests_remap[m]);
  }

  // perform all sends to other procs

  for (m = 0; m < nsend_remap; m++) {
    ptr->pack_remap_grid(which,buf1,send_remap[m].npack,send_remap[m].packlist);
    MPI_Send(buf1,nper*send_remap[m].npack,datatype,send_remap[m].proc,0,gridcomm);
  }

  // perform remap to self if defined

  if (self_remap) {
    ptr->pack_remap_grid(which,buf1,copy_remap.npack,copy_remap.packlist);
    ptr->unpack_remap_grid(which,buf1,copy_remap.nunpack,copy_remap.unpacklist);
  }

  // unpack all received data

  for (i = 0; i < nrecv_remap; i++) {
    MPI_Waitany(nrecv_remap,requests_remap,&m,MPI_STATUS_IGNORE);
    offset = nper * recv_remap[m].offset * nbyte;
    ptr->unpack_remap_grid(which,(void *) &buf2[offset],
                           recv_remap[m].nunpack,recv_remap[m].unpacklist);
  }
}

// ----------------------------------------------------------------------
// grid I/O methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   read grid values from a file
------------------------------------------------------------------------- */

void Grid3d::read_file(int caller, void *ptr, FILE *fp, int nchunk, int maxline)
{
  if (caller == FIX)
    read_file_style<Fix>((Fix *) ptr,fp,nchunk,maxline);
}

/* ----------------------------------------------------------------------
   proc 0 reads one chunk of lines at a time from file
   broadcast chunk buffer to other procs
   call back to caller so it can process the chunk of lines
   caller returns count of grid-value lines in chunk
------------------------------------------------------------------------- */

template < class T >
void Grid3d::read_file_style(T *ptr, FILE *fp, int nchunk, int maxline)
{
  auto buffer = new char[nchunk * maxline];
  bigint ntotal = (bigint) nx * ny * nz;
  bigint nread = 0;

  while (nread < ntotal) {
    int mychunk = MIN(ntotal - nread, nchunk);
    int eof = utils::read_lines_from_file(fp, mychunk, maxline, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of grid data file");

    nread += ptr->unpack_read_grid(nchunk,buffer);
  }

  delete [] buffer;
}

/* ----------------------------------------------------------------------
   write grid values to a file
------------------------------------------------------------------------- */

void Grid3d::write_file(int caller, void *ptr, int which,
                        int nper, int nbyte, MPI_Datatype datatype)
{
  if (caller == FIX)
    write_file_style<Fix>((Fix *) ptr, which, nper, nbyte, datatype);
}

/* ----------------------------------------------------------------------
   proc 0 reads one chunk of lines at a time from file
   broadcast chunk buffer to other procs
   call back to caller so it can process the chunk of lines
   caller returns count of grid-value lines in chunk
------------------------------------------------------------------------- */

template < class T >
void Grid3d::write_file_style(T *ptr, int which,
                              int nper, int nbyte, MPI_Datatype datatype)
{
  // maxsize = max size of grid data owned by any proc

  int mysize = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1);
  mysize *= nper;
  int maxsize;
  MPI_Allreduce(&mysize,&maxsize,1,MPI_INT,MPI_MAX,world);

  // pack my grid data via callback to caller

  char *onebuf;
  if (me == 0) memory->create(onebuf,maxsize*nbyte,"grid3d:onebuf");
  else memory->create(onebuf,mysize*nbyte,"grid3d:nebuf");
  ptr->pack_write_grid(which,onebuf);

  // ping each proc for its grid data
  // call back to caller with each proc's grid data

  int tmp;
  int bounds[6];

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(onebuf,maxsize,datatype,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Recv(bounds,6,MPI_INT,iproc,0,world,&status);
      } else {
        bounds[0] = inxlo;
        bounds[1] = inxhi;
        bounds[2] = inylo;
        bounds[3] = inyhi;
        bounds[4] = inzlo;
        bounds[5] = inzhi;
      }

      ptr->unpack_write_grid(which,onebuf,bounds);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(onebuf,mysize,datatype,0,0,world);
    bounds[0] = inxlo;
    bounds[1] = inxhi;
    bounds[2] = inylo;
    bounds[3] = inyhi;
    bounds[4] = inzlo;
    bounds[5] = inzhi;
    MPI_Send(bounds,6,MPI_INT,0,0,world);
  }

  // clean up

  memory->destroy(onebuf);
}

// ----------------------------------------------------------------------
// overlap methods for brick and tiled RCB decompositions
// overlap = overlap of owned or owned+ghost box with all boxes of a decomposition
// for owned/ghost grid comm, called only by tiled decomposition
//   brick decomp uses one or more comm passes with neigh procs
//   like forward/reverse comm for atoms
// for remap, called by both brick and tiled decompositions
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   compute list of overlaps between box and the owned grid boxes of all procs
   ghostflag = 1 if box includes ghost grid pts, called by setup_tiled()
   ghostflag = 0 if box has no ghost grid pts, called by setup_remap()
   layout != LAYOUT_TILED is only invoked by setup_remap()
   for brick decomp of Grid, done using xyz split + grid2proc copied from Comm
   for tiled decomp of Grid, done via recursive box drop on RCB tree
   box = 6 integers = (xlo,xhi,ylo,yhi,zlo,zhi)
     box can be owned cells or owned + ghost cells
   pbc = flags for grid periodicity in each dim
     if box includes ghost cells, it can overlap PBCs (only for setup_tiled)
     each lo/hi value may extend beyond 0 to N-1 into another periodic image
   return # of overlaps including with self, caller handles self overlaps as needed
   return list of overlaps
     for setup_tiled() this is what box_drop() computes
       entire box for each overlap
       caller will determine extent of overlap using PBC info
     for setup_remap(), return extent of overlap (no PBC info involved)
       use proc_box_uniform() or tiled() and MAX/MIN to determine this
------------------------------------------------------------------------- */

int Grid3d::compute_overlap(int ghostflag, int *box, int *pbc, Overlap *&overlap)
{
  int obox[6];

  memory->create(overlap_procs,nprocs,"grid3d:overlap_procs");
  noverlap_list = maxoverlap_list = 0;
  overlap_list = nullptr;

    // skip overlap check if box contains no grid cells

  if (box[0] > box[1] || box[2] > box[3] || box[4] > box[5]) {
    overlap = overlap_list;
    return noverlap_list;
  }

  if (comm->layout != Comm::LAYOUT_TILED) {

    // find comm->procgrid indices in each dim for box bounds

    int iproclo = proc_index_uniform(box[0],nx,shift_grid,0,xsplit);
    int iprochi = proc_index_uniform(box[1],nx,shift_grid,0,xsplit);
    int jproclo = proc_index_uniform(box[2],ny,shift_grid,1,ysplit);
    int jprochi = proc_index_uniform(box[3],ny,shift_grid,1,ysplit);
    int kproclo = proc_index_uniform(box[4],nz,shift_grid,2,zsplit);
    int kprochi = proc_index_uniform(box[5],nz,shift_grid,2,zsplit);

    // compute extent of overlap of box with with each proc's obox

    for (int k = kproclo; k <= kprochi; k++)
      for (int j = jproclo; j <= jprochi; j++)
        for (int i = iproclo; i <= iprochi; i++) {
          partition_grid(nx,xsplit[i],xsplit[i+1],shift_grid,0,obox[0],obox[1]);
          partition_grid(ny,ysplit[j],ysplit[j+1],shift_grid,0,obox[2],obox[3]);
          partition_grid(nz,zsplit[k],zsplit[k+1],shift_grid,zextra,obox[4],obox[5]);

          if (noverlap_list == maxoverlap_list) grow_overlap();
          overlap_list[noverlap_list].proc = grid2proc[i][j][k];
          overlap_list[noverlap_list].box[0] = MAX(box[0],obox[0]);
          overlap_list[noverlap_list].box[1] = MIN(box[1],obox[1]);
          overlap_list[noverlap_list].box[2] = MAX(box[2],obox[2]);
          overlap_list[noverlap_list].box[3] = MIN(box[3],obox[3]);
          overlap_list[noverlap_list].box[4] = MAX(box[4],obox[4]);
          overlap_list[noverlap_list].box[5] = MIN(box[5],obox[5]);
          noverlap_list++;
        }

  } else {
    box_drop(box,pbc);

    // compute extent of overlap of box with with each proc's obox

    if (ghostflag == 0) {
      for (int m = 0; m < noverlap_list; m++) {
        obox[0] = 0;
        obox[1] = nx-1;
        obox[2] = 0;
        obox[3] = ny-1;
        obox[4] = 0;
        obox[5] = nz-1;

        partition_tiled(overlap_list[m].proc,0,nprocs-1,obox);

        overlap_list[m].box[0] = MAX(box[0],obox[0]);
        overlap_list[m].box[1] = MIN(box[1],obox[1]);
        overlap_list[m].box[2] = MAX(box[2],obox[2]);
        overlap_list[m].box[3] = MIN(box[3],obox[3]);
        overlap_list[m].box[4] = MAX(box[4],obox[4]);
        overlap_list[m].box[5] = MIN(box[5],obox[5]);
      }
    }
  }

  overlap = overlap_list;
  return noverlap_list;
}

/* ----------------------------------------------------------------------
   deallocate data created by recursive overlap computation
------------------------------------------------------------------------- */

void Grid3d::clean_overlap()
{
  memory->destroy(overlap_procs);
  memory->sfree(overlap_list);
}

/* ----------------------------------------------------------------------
   recursively split a box until it doesn't overlap any periodic boundaries
   box = 6 integers = (xlo,xhi,ylo,yhi,zlo,zhi)
     each lo/hi value may extend beyonw 0 to N-1 into another periodic image
   pbc = flags in each dim of which periodic image the caller box was in
   when a box straddles a periodic bounadry, split it in two
   when a box does not straddle, drop it down RCB tree
     add all the procs it overlaps with to Overlap list
------------------------------------------------------------------------- */

void Grid3d::box_drop(int *box, int *pbc)
{
  int i,m;

  // newbox12 and newpbc are initially copies of caller box and pbc

  int newbox1[6],newbox2[6],newpbc[3];

  for (i = 0; i < 6; i++) newbox1[i] = newbox2[i] = box[i];
  for (i = 0; i < 3; i++) newpbc[i] = pbc[i];

  // 6 if tests to see if box needs to be split across a periodic boundary
  // newbox1 and 2 = new split boxes, newpbc increments current pbc
  // final else is no split

  int splitflag = 1;

  if (box[0] < 0) {
    newbox1[0] = 0;
    newbox2[0] = box[0] + nx;
    newbox2[1] = nx - 1;
    newpbc[0]--;
  } else if (box[1] >= nx) {
    newbox1[1] = nx - 1;
    newbox2[0] = 0;
    newbox2[1] = box[1] - nx;
    newpbc[0]++;
  } else if (box[2] < 0) {
    newbox1[2] = 0;
    newbox2[2] = box[2] + ny;
    newbox2[3] = ny - 1;
    newpbc[1]--;
  } else if (box[3] >= ny) {
    newbox1[3] = ny - 1;
    newbox2[2] = 0;
    newbox2[3] = box[3] - ny;
    newpbc[1]++;
  } else if (box[4] < 0) {
    newbox1[4] = 0;
    newbox2[4] = box[4] + nz;
    newbox2[5] = nz - 1;
    newpbc[2]--;
  } else if (box[5] >= nz) {
    newbox1[5] = nz - 1;
    newbox2[4] = 0;
    newbox2[5] = box[5] - nz;
    newpbc[2]++;

  // box is not split, drop on RCB tree
  // returns np = # of procs it overlaps, including self
  // returns proc_overlap = list of proc IDs it overlaps
  // add each overlap to overlap list

  } else {
    splitflag = 0;
    int np = 0;
    box_drop_grid(box,0,nprocs-1,np,overlap_procs);
    for (m = 0; m < np; m++) {
      if (noverlap_list == maxoverlap_list) grow_overlap();
      overlap_list[noverlap_list].proc = overlap_procs[m];
      for (i = 0; i < 6; i++) overlap_list[noverlap_list].box[i] = box[i];
      for (i = 0; i < 3; i++) overlap_list[noverlap_list].pbc[i] = pbc[i];
      noverlap_list++;
    }
  }

  // recurse with 2 split boxes

  if (splitflag) {
    box_drop(newbox1,pbc);
    box_drop(newbox2,newpbc);
  }
}

/* ----------------------------------------------------------------------
   recursively drop a box down the RCB tree to find all procs it overlaps with
   box = 6 integers = (xlo,xhi,ylo,yhi,zlo,zhi)
     each lo/hi value ranges from 0 to N-1 in a dim, N = grid size in that dim
     box is guaranteed to be wholly within the global domain
   return Np = # of procs, plist = proc IDs
------------------------------------------------------------------------- */

void Grid3d::box_drop_grid(int *box, int proclower, int procupper,
                           int &np, int *plist)
{
  // end recursion when partition is a single proc
  // add proclower to plist

  if (proclower == procupper) {
    plist[np++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use < and >= criteria so does not include a box it only touches
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // cut = index of first grid cell in upper partition
  // dim = 0,1,2 dimension of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int dim = rcbinfo[procmid].dim;
  int cut = rcbinfo[procmid].cut;

  if (box[2*dim] < cut) box_drop_grid(box,proclower,procmid-1,np,plist);
  if (box[2*dim+1] >= cut) box_drop_grid(box,procmid,procupper,np,plist);
}

// ----------------------------------------------------------------------
// miscellaneous methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   grow list of swaps by DELTA
------------------------------------------------------------------------- */

void Grid3d::grow_swap()
{
  maxswap += DELTA;
  swap = (Swap *) memory->srealloc(swap,maxswap*sizeof(Swap),"grid3d:swap");
}

/* ----------------------------------------------------------------------
   grow list of overlaps by DELTA
------------------------------------------------------------------------- */

void Grid3d::grow_overlap()
{
  maxoverlap_list += DELTA;
  overlap_list = (Overlap *)
    memory->srealloc(overlap_list,maxoverlap_list*sizeof(Overlap),"grid3d:overlap_list");
}

/* ----------------------------------------------------------------------
   deallocate remap data structs
------------------------------------------------------------------------- */

void Grid3d::deallocate_remap()
{
  for (int i = 0; i < nsend_remap; i++)
    memory->destroy(send_remap[i].packlist);
  delete [] send_remap;

  for (int i = 0; i < nrecv_remap; i++)
    memory->destroy(recv_remap[i].unpacklist);
  delete [] recv_remap;

  if (self_remap) {
    memory->destroy(copy_remap.packlist);
    memory->destroy(copy_remap.unpacklist);
  }
}

/* ----------------------------------------------------------------------
   create 1d list of offsets into 3d array section (xlo:xhi,ylo:yhi,zlo:zhi)
   assume caller's 3d array is allocated as
     (fullxlo:fullxhi,fullylo:fullyhi,fullzlo:fullzhi)
------------------------------------------------------------------------- */

int Grid3d::indices(int *&list,
                    int xlo, int xhi, int ylo, int yhi, int zlo, int zhi)
{
  int nmax = (xhi-xlo+1) * (yhi-ylo+1) * (zhi-zlo+1);
  memory->create(list,nmax,"grid3d:indices");
  if (nmax == 0) return 0;

  int nx = (fullxhi-fullxlo+1);
  int ny = (fullyhi-fullylo+1);

  int n = 0;
  int ix,iy,iz;
  for (iz = zlo; iz <= zhi; iz++)
    for (iy = ylo; iy <= yhi; iy++)
      for (ix = xlo; ix <= xhi; ix++)
        list[n++] = (iz-fullzlo)*ny*nx + (iy-fullylo)*nx + (ix-fullxlo);

  return nmax;
}

/* ----------------------------------------------------------------------
   find the comm->procgrid index = which proc owns the igrid index
   igrid = grid index (0 to N-1) in dim
   n = # of grid points in dim
   shift determines position of grid pt within grid cell
     shift = 0.5 for cell center, 0.0 for lower-left corner
   dim = which dimension (0,1,2)
   split = comm->x/y/z split for fractional bounds of each proc domain
------------------------------------------------------------------------- */

int Grid3d::proc_index_uniform(int igrid, int n, double shift, int dim, double *split)
{
  int lo,hi;
  double fraclo,frachi;

  // loop over # of procs in this dime
  // compute the grid bounds for that proc
  // if igrid falls within those bounds, return m = proc index
  // same logic as in partition_grid()

  int m;
  for (m = 0; m < comm->procgrid[dim]; m++) {
    fraclo = split[m];
    frachi = split[m+1];

    lo = static_cast<int> (fraclo * n);
    while (lo+shift < fraclo*n) lo++;
    hi = static_cast<int> (frachi * n);
    if (hi+shift >= frachi*n) hi--;

    if (igrid >= lo && igrid <= hi) break;
  }

  return m;
}

/* ----------------------------------------------------------------------
   compute the grid box owned by proc within tiled decomposition
   performed recursively until proclower = procupper = proc
   return box = lo/hi bounds of proc's box in 3 dims
------------------------------------------------------------------------- */

void Grid3d::partition_tiled(int proc, int proclower, int procupper, int *box)
{
  // end recursion when partition is a single proc

  if (proclower == procupper) return;

  // split processor partition
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // cut = index of first grid cell in upper partition
  // dim = 0,1,2 dimension of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int dim = rcbinfo[procmid].dim;
  int cut = rcbinfo[procmid].cut;

  // adjust box to reflect which half of partition the proc is in

  if (proc < procmid) {
    box[2*dim+1] = cut-1;
    partition_tiled(proc,proclower,procmid-1,box);
  } else {
    box[2*dim] = cut;
    partition_tiled(proc,procmid,procupper,box);
  }
}
