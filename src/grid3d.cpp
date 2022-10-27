// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

using namespace LAMMPS_NS;

#define DELTA 16

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
   constructor called by all classes except PPPM and MSM
   gcomm = world communicator
   gnx, gny, gnz = size of global grid
   maxdist = max distance outside of proc domain a particle will be
   extra = additional ghost grid pts needed in each dim, e.g. for stencil
   shift = 0.0 for grid pt in lower-left corner of grid cell, 0.5 for center
   return:
     i xyz lohi = portion of global grid this proc owns, 0 <= index < N
     o xyz lohi = owned + ghost grid cells needed in all directions
   for non-periodic dims, o indices will not be < 0 or >= N,
     since no grid comm is done across non-periodic boundaries
------------------------------------------------------------------------- */

Grid3d::Grid3d(LAMMPS *lmp, MPI_Comm gcomm,
               int gnx, int gny, int gnz,
               double maxdist, int extra, double shift,
               int &ixlo, int &ixhi, int &iylo, int &iyhi, int &izlo, int &izhi,
               int &oxlo, int &oxhi, int &oylo, int &oyhi, int &ozlo, int &ozhi)
  : Pointers(lmp)
{
  // store commnicator and global grid size
  // set layout mode

  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;
  nz = gnz;

  ngrid[0] = nx; ngrid[1] = ny; ngrid[2] = nz;
  layout = comm->layout;

  // partition global grid across procs
  // i xyz lo/hi = lower/upper bounds of global grid this proc owns
  // indices range from 0 to N-1 inclusive in each dim

  comm->partition_grid(nx, ny, nz, 0.0, ixlo, ixhi, iylo, iyhi, izlo, izhi);

  // nlo,nhi = min/max index of global grid pt my owned atoms can be mapped to
  // finite difference stencil requires extra grid pt around my owned grid pts
  // max of these 2 quantities is the ghost cells needed in each dim
  // o xyz lo/hi = owned + ghost cells

  int triclinic = domain->triclinic;

  double *prd,*boxlo,*sublo,*subhi;

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

  double dist[3] = {0.0,0.0,0.0};
  if (triclinic == 0) dist[0] = dist[1] = dist[2] = maxdist;
  else MathExtra::tribbox(domain->h,maxdist,&dist[0]);

  double dxinv = nx / prd[0];
  double dyinv = ny / prd[1];
  double dzinv = nz / prd[2];
  double SHIFT = OFFSET + shift;
  int nlo, nhi;

  nlo = static_cast<int>((sublo[0]-dist[0]-boxlo[0]) * dxinv + SHIFT) - OFFSET;
  nhi = static_cast<int>((subhi[0]+dist[0]-boxlo[0]) * dxinv + SHIFT) - OFFSET;
  oxlo = MIN(nlo, ixlo - extra);
  oxhi = MAX(nhi, ixhi + extra);

  nlo = static_cast<int>((sublo[1]-dist[1]-boxlo[1]) * dyinv + SHIFT) - OFFSET;
  nhi = static_cast<int>((subhi[1]+dist[1]-boxlo[1]) * dyinv + SHIFT) - OFFSET;
  oylo = MIN(nlo, iylo - extra);
  oyhi = MAX(nhi, iyhi + extra);

  nlo = static_cast<int>((sublo[2]-dist[2]-boxlo[2]) * dzinv + SHIFT) - OFFSET;
  nhi = static_cast<int>((subhi[2]+dist[2]-boxlo[2]) * dzinv + SHIFT) - OFFSET;
  ozlo = MIN(nlo, izlo - extra);
  ozhi = MAX(nhi, izhi + extra);

  // limit o xyz lo/hi indices for non-periodic dimensions

  int *periodicity = domain->periodicity;

  if (!periodicity[0]) {
    oxlo = MAX(0,oxlo);
    oxhi = MIN(nx-1,oxhi);
  }

  if (!periodicity[1]) {
    oylo = MAX(0,oylo);
    oyhi = MIN(ny-1,oyhi);
  }

  if (!periodicity[2]) {
    ozlo = MAX(0,ozlo);
    ozhi = MIN(nz-1,ozhi);
  }

  // error check on size of grid stored by this proc

  bigint total = (bigint)
    (oxhi - oxlo + 1) * (oyhi - oylo + 1) * (ozhi - ozlo + 1);
  if (total > MAXSMALLINT)
    error->one(FLERR, "Too many owned+ghost grid3d points");

  // store grid bounds and proc neighs

  if (layout != Comm::LAYOUT_TILED) {
    int (*procneigh)[2] = comm->procneigh;
    store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          procneigh[0][0],procneigh[0][1],
          procneigh[1][0],procneigh[1][1],
          procneigh[2][0],procneigh[2][1]);
  } else {
    store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          0,0,0,0,0,0);
  }
}

/* ----------------------------------------------------------------------
   constructor called by PPPM classes
   gcomm = world communicator
   gnx, gny, gnz = size of global grid
   i xyz lohi = portion of global grid this proc owns, 0 <= index < N
   o xyz lohi = owned grid portion + ghost grid cells needed in all directions
   if o indices are < 0 or hi indices are >= N,
     then grid is treated as periodic in that dimension,
     communication is done across the periodic boundaries
------------------------------------------------------------------------- */

Grid3d::Grid3d(LAMMPS *lmp, MPI_Comm gcomm,
               int gnx, int gny, int gnz,
               int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
               int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi)
  : Pointers(lmp)
{
  // store commnicator and global grid size
  // set layout mode

  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;
  nz = gnz;

  ngrid[0] = nx; ngrid[1] = ny; ngrid[2] = nz;
  layout = comm->layout;

  // store grid bounds and proc neighs

  if (layout != Comm::LAYOUT_TILED) {
    int (*procneigh)[2] = comm->procneigh;
    store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          procneigh[0][0],procneigh[0][1],
          procneigh[1][0],procneigh[1][1],
          procneigh[2][0],procneigh[2][1]);
  } else {
    store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
          0,0,0,0,0,0);
  }
}

/* ----------------------------------------------------------------------
   constructor called by MSM
   gcomm = world communicator or sub-communicator for a hierarchical grid
   flag = 1 if e xyz lohi values = larger grid stored by caller in gcomm = world
   flag = 2 if e xyz lohi values = 6 neighbor procs in gcomm
   gnx, gny, gnz = size of global grid
   i xyz lohi = portion of global grid this proc owns, 0 <= index < N
   o xyz lohi = owned grid portion + ghost grid cells needed in all directions
   e xyz lohi for flag = 1: extent of larger grid stored by caller
   e xyz lohi for flag = 2: 6 neighbor procs
------------------------------------------------------------------------- */

Grid3d::Grid3d(LAMMPS *lmp, MPI_Comm gcomm, int flag,
               int gnx, int gny, int gnz,
               int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
               int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi,
               int exlo, int exhi, int eylo, int eyhi, int ezlo, int ezhi)
  : Pointers(lmp)
{
  // store commnicator and global grid size
  // set layout mode

  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;
  nz = gnz;

  ngrid[0] = nx; ngrid[1] = ny; ngrid[2] = nz;
  layout = comm->layout;

  // store grid bounds and proc neighs

  if (flag == 1) {
    if (layout != Comm::LAYOUT_TILED) {
      // this assumes gcomm = world
      int (*procneigh)[2] = comm->procneigh;
      store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
            oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
            exlo,exhi,eylo,eyhi,ezlo,ezhi,
            procneigh[0][0],procneigh[0][1],
            procneigh[1][0],procneigh[1][1],
            procneigh[2][0],procneigh[2][1]);
    } else {
      store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
            oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
            exlo,exhi,eylo,eyhi,ezlo,ezhi,
            0,0,0,0,0,0);
    }

  } else if (flag == 2) {
    if (layout != Comm::LAYOUT_TILED) {
      store(ixlo,ixhi,iylo,iyhi,izlo,izhi,
            oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
            oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
            exlo,exhi,eylo,eyhi,ezlo,ezhi);
    } else {
      error->all(FLERR,"Grid3d does not support tiled layout with neighbor procs");
    }
  }
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
// store and access Grid parameters
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   store grid bounds and proc neighs in local variables
------------------------------------------------------------------------- */

void Grid3d::store(int ixlo, int ixhi, int iylo, int iyhi,
                   int izlo, int izhi,
                   int oxlo, int oxhi, int oylo, int oyhi,
                   int ozlo, int ozhi,
                   int fxlo, int fxhi, int fylo, int fyhi,
                   int fzlo, int fzhi,
                   int pxlo, int pxhi, int pylo, int pyhi,
                   int pzlo, int pzhi)
{
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

  fullxlo = fxlo;
  fullxhi = fxhi;
  fullylo = fylo;
  fullyhi = fyhi;
  fullzlo = fzlo;
  fullzhi = fzhi;

  // internal data initializations

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
  
  // for non TILED layout:
  // proc xyz lohi = my 6 neighbor procs in this MPI_Comm
  // xyz split = copy of 1d vectors in Comm
  // grid2proc = copy of 3d array in Comm

  if (layout != Comm::LAYOUT_TILED) {
    procxlo = pxlo;
    procxhi = pxhi;
    procylo = pylo;
    procyhi = pyhi;
    proczlo = pzlo;
    proczhi = pzhi;

    xsplit = new double[comm->procgrid[0]+1];
    ysplit = new double[comm->procgrid[1]+1];
    zsplit = new double[comm->procgrid[2]+1];
    memcpy(xsplit,comm->xsplit,(comm->procgrid[0]+1)*sizeof(double));
    memcpy(ysplit,comm->ysplit,(comm->procgrid[1]+1)*sizeof(double));
    memcpy(zsplit,comm->zsplit,(comm->procgrid[2]+1)*sizeof(double));
    
    memory->create(grid2proc,comm->procgrid[0],comm->procgrid[1],comm->procgrid[2],
                   "grid3d:grid2proc");
    memcpy(&grid2proc[0][0][0],&comm->grid2proc[0][0][0],
           comm->procgrid[0]*comm->procgrid[1]*comm->procgrid[2]*sizeof(int));
  }

  // for TILED layout:
  // create RCB tree of cut info for grid decomp
  // access CommTiled to get cut dimension
  // cut = this proc's inlo in that dim
  // dim is -1 for proc 0, but never accessed

  if (layout == Comm::LAYOUT_TILED) {
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

/* ---------------------------------------------------------------------- */

int Grid3d::identical(Grid3d *grid2)
{
  int inxlo2,inxhi2,inylo2,inyhi2,inzlo2,inzhi2;
  int outxlo2,outxhi2,outylo2,outyhi2,outzlo2,outzhi2;
  
  grid2->get_bounds(inxlo2,inxhi2,inylo2,inyhi2,inzlo2,inzhi2);
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

void Grid3d::get_bounds(int &xlo, int &xhi, int &ylo, int &yhi,
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
// setup of local owned/ghost grid comm
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   setup owned/ghost commmunication
   return sizes of two buffers needed for communication
   either for brick decomp or tiling decomp
   nbuf1 = largest pack or unpack in any Send or Recv or Copy
   nbuf2 = larget of sum of all packs or unpacks in Send or Recv
   for brick comm, nbuf1 = nbuf2
   for tiling comm, nbuf2 >= nbuf2
   nbuf1,nbuf2 are counts of grid points
     caller converts them to message sizes for grid data it stores
------------------------------------------------------------------------- */

void Grid3d::setup(int &nbuf1, int &nbuf2)
{
  if (layout != Comm::LAYOUT_TILED) setup_brick(nbuf1,nbuf2);
  else setup_tiled(nbuf1,nbuf2);
}

/* ----------------------------------------------------------------------
   setup owned/ghost comm for brick comm
   each proc has 6 neighbors
   comm pattern = series of swaps with one of those 6 procs
   can be multiple swaps with same proc if ghost extent is large
   swap may not be symmetric if both procs do not need same layers of ghosts
   all procs perform same # of swaps in a direction, even if some don't need it
------------------------------------------------------------------------- */

void Grid3d::setup_brick(int &nbuf1, int &nbuf2)
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

void Grid3d::setup_tiled(int &nbuf1, int &nbuf2)
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
// query locality of forwrd/reverse grid comm
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   check if all procs only need ghost info from adjacent procs
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int Grid3d::ghost_adjacent()
{
  if (layout != Comm::LAYOUT_TILED) return ghost_adjacent_brick();
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

void Grid3d::forward_comm(int caller, void *ptr, int nper, int nbyte, int which,
                            void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (layout != Comm::LAYOUT_TILED) {
    if (caller == KSPACE)
      forward_comm_brick<KSpace>((KSpace *) ptr,nper,nbyte,which,
				 buf1,buf2,datatype);
    else if (caller == PAIR)
      forward_comm_brick<Pair>((Pair *) ptr,nper,nbyte,which,
			       buf1,buf2,datatype);
    else if (caller == FIX)
      forward_comm_brick<Fix>((Fix *) ptr,nper,nbyte,which,
			      buf1,buf2,datatype);
  } else {
    if (caller == KSPACE)
      forward_comm_tiled<KSpace>((KSpace *) ptr,nper,nbyte,which,
                                 buf1,buf2,datatype);
    else if (caller == PAIR)
      forward_comm_tiled<Pair>((Pair *) ptr,nper,nbyte,which,
                                 buf1,buf2,datatype);
    else if (caller == FIX)
      forward_comm_tiled<Fix>((Fix *) ptr,nper,nbyte,
                              which,buf1,buf2,datatype);
  }
}

/* ----------------------------------------------------------------------
   forward comm for brick decomp via list of swaps with 6 neighbor procs
------------------------------------------------------------------------- */

template < class T >
void Grid3d::
forward_comm_brick(T *ptr, int nper, int /*nbyte*/, int which,
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
forward_comm_tiled(T *ptr, int nper, int nbyte, int which,
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

void Grid3d::reverse_comm(int caller, void *ptr, int nper, int nbyte, int which,
                            void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (layout != Comm::LAYOUT_TILED) {
    if (caller == KSPACE)
      reverse_comm_brick<KSpace>((KSpace *) ptr,nper,nbyte,which,
				buf1,buf2,datatype);
    else if (caller == PAIR)
      reverse_comm_brick<Pair>((Pair *) ptr,nper,nbyte,which,
			       buf1,buf2,datatype);
    else if (caller == FIX)
      reverse_comm_brick<Fix>((Fix *) ptr,nper,nbyte,which,
			      buf1,buf2,datatype);
  } else {
    if (caller == KSPACE)
      reverse_comm_tiled<KSpace>((KSpace *) ptr,nper,nbyte,which,
                                 buf1,buf2,datatype);
    else if (caller == PAIR)
      reverse_comm_tiled<Pair>((Pair *) ptr,nper,nbyte,which,
                                 buf1,buf2,datatype);
    else if (caller == FIX)
      reverse_comm_tiled<Fix>((Fix *) ptr,nper,nbyte,which,
                              buf1,buf2,datatype);
  }
}

/* ----------------------------------------------------------------------
   reverse comm for brick decomp via list of swaps with 6 neighbor procs
------------------------------------------------------------------------- */

template < class T >
void Grid3d::
reverse_comm_brick(T *ptr, int nper, int /*nbyte*/, int which,
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
reverse_comm_tiled(T *ptr, int nper, int nbyte, int which,
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

  // set layout to current Comm layout

  layout = comm->layout;
  
  // overlaps of my old decomp owned box with all owned boxes in new decomp
  // noverlap_old = # of overlaps, including self
  // overlap_old = vector of overlap info in Overlap data struct

  int oldbox[6];
  old->get_bounds(oldbox[0],oldbox[1],oldbox[2],oldbox[3],oldbox[4],oldbox[5]);
  pbc[0] = pbc[1] = pbc[2] = 0;

  Overlap *overlap_old;
  int noverlap_old = compute_overlap(0,oldbox,pbc,overlap_old);

  // use overlap_old to construct send and copy lists
  
  self_remap = 0;

  nsend_remap = 0;
  for (m = 0; m < noverlap_old; m++) {
    if (overlap_old[m].proc == me) self_remap = 1;
    else nsend_remap++;
  }

  send_remap = new Send[nsend_remap];
  
  nsend_remap = 0;
  for (m = 0; m < noverlap_old; m++) {
    box = overlap_old[m].box;
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
  get_bounds(newbox[0],newbox[1],newbox[2],newbox[3],newbox[4],newbox[5]);
  pbc[0] = pbc[1] = pbc[2] = 0;

  Overlap *overlap_new;
  int noverlap_new = old->compute_overlap(0,newbox,pbc,overlap_new);

  // use overlap_new to construct recv and copy lists
  // set offsets for Recv data

  nrecv_remap = 0;
  for (m = 0; m < noverlap_new; m++)
    if (overlap_new[m].proc != me) nrecv_remap++;

  recv_remap = new Recv[nrecv_remap];

  nrecv_remap = 0;
  for (m = 0; m < noverlap_new; m++) {
    box = overlap_new[m].box;
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

void Grid3d::remap(int caller, void *ptr, int nper, int nbyte,
		   void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (caller == FIX) remap_style<Fix>((Fix *) ptr,nper,nbyte,buf1,buf2,datatype);
}

/* ------------------------------------------------------------------------- */

template < class T >
void Grid3d::remap_style(T *ptr, int nper, int nbyte,
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
    ptr->pack_remap_grid(buf1,send_remap[m].npack,send_remap[m].packlist);
    MPI_Send(buf1,nper*send_remap[m].npack,datatype,send_remap[m].proc,0,gridcomm);
  }

  // perform remap to self if defined

  if (self_remap) {
    ptr->pack_remap_grid(buf1,copy_remap.npack,copy_remap.packlist);
    ptr->unpack_remap_grid(buf1,copy_remap.nunpack,copy_remap.unpacklist);
  }

  // unpack all received data

  for (i = 0; i < nrecv_remap; i++) {
    MPI_Waitany(nrecv_remap,requests_remap,&m,MPI_STATUS_IGNORE);
    offset = nper * recv_remap[m].offset * nbyte;
    ptr->unpack_remap_grid((void *) &buf2[offset],
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
  bigint ntotal = (bigint) ngrid[0] * ngrid[1] * ngrid[2];
  bigint nread = 0;

  while (nread < ntotal) {
    int nchunk = MIN(ntotal - nread, nchunk);
    int eof = utils::read_lines_from_file(fp, nchunk, maxline, buffer, me, world);
    if (eof) error->all(FLERR, "Unexpected end of grid data file");

    nread += ptr->unpack_read_grid(buffer);
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

  if (layout != Comm::LAYOUT_TILED) {

    // find comm->procgrid indices in each dim for box bounds
    
    int iproclo = proc_index_uniform(box[0],ngrid[0],0,xsplit);
    int iprochi = proc_index_uniform(box[1],ngrid[0],0,xsplit);
    int jproclo = proc_index_uniform(box[2],ngrid[1],1,ysplit);
    int jprochi = proc_index_uniform(box[3],ngrid[1],1,ysplit);
    int kproclo = proc_index_uniform(box[4],ngrid[2],2,zsplit);
    int kprochi = proc_index_uniform(box[5],ngrid[2],2,zsplit);

    // compute extent of overlap of box with with each proc's obox
    
    for (int k = kproclo; k <= kprochi; k++)
      for (int j = jproclo; j <= jprochi; j++)
	for (int i = iproclo; i <= iprochi; i++) {
	  proc_box_uniform(i,j,k,obox);

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
        obox[1] = ngrid[0]-1;
        obox[2] = 0;
        obox[3] = ngrid[1]-1;
        obox[4] = 0;
        obox[5] = ngrid[2]-1;
        
        proc_box_tiled(overlap_list[m].proc,0,nprocs-1,obox);
        
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
   assume 3d array is allocated as
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
   dim = which dimension (0,1,2)
   split = comm->x/y/z split for fractional bounds of each proc domain
------------------------------------------------------------------------- */

int Grid3d::proc_index_uniform(int igrid, int n, int dim, double *split)
{
  int gridlo,gridhi;
  double fraclo,frachi;

  // loop over # of procs in this dime
  // compute the grid bounds for that proc, same as comm->partition_grid()
  // if igrid falls within those bounds, return m = proc index

  int m;
  for (m = 0; m < comm->procgrid[dim]; m++) {
    fraclo = split[m];
    frachi = split[m+1];
    gridlo = static_cast<int> (fraclo * n);
    if (1.0*gridlo != fraclo*n) gridlo++;
    gridhi = static_cast<int> (frachi * n);
    if (1.0*gridhi == frachi*n) gridhi--;

    if (igrid >= gridlo && igrid <= gridhi) break;
  }

  return m;
}

/* ----------------------------------------------------------------------
   compute the grid box for proc with grid indices i,j,k
   i,j,k = grid index (0 to N-1) in each dim
   return lo/hi bounds of box in 3 dims
   computation is same as Comm::partition_grid()
------------------------------------------------------------------------- */

void Grid3d::proc_box_uniform(int i, int j, int k, int *box)
{
  int lo,hi;
  double fraclo,frachi;
  
  fraclo = xsplit[i];
  frachi = xsplit[i+1];
  lo = static_cast<int> (fraclo * ngrid[0]);
  if (1.0*lo != fraclo*ngrid[0]) lo++;
  hi = static_cast<int> (frachi * ngrid[0]);
  if (1.0*hi == frachi*ngrid[0]) hi--;
  box[0] = lo;
  box[1] = hi;
  
  fraclo = ysplit[j];
  frachi = ysplit[j+1];
  lo = static_cast<int> (fraclo * ngrid[1]);
  if (1.0*lo != fraclo*ngrid[1]) lo++;
  hi = static_cast<int> (frachi * ngrid[1]);
  if (1.0*hi == frachi*ngrid[1]) hi--;
  box[2] = lo;
  box[3] = hi;

  fraclo = zsplit[k];
  frachi = zsplit[k+1];
  lo = static_cast<int> (fraclo * ngrid[2]);
  if (1.0*lo != fraclo*ngrid[2]) lo++;
  hi = static_cast<int> (frachi * ngrid[2]);
  if (1.0*hi == frachi*ngrid[2]) hi--;
  box[4] = lo;
  box[5] = hi;
}

/* ----------------------------------------------------------------------
   compute the grid box for proc within tiled decomposition
   performed recursively until proclower = procupper = proc
   return box = lo/hi bounds of proc's box in 3 dims
------------------------------------------------------------------------- */

void Grid3d::proc_box_tiled(int proc, int proclower, int procupper, int *box)
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
    proc_box_tiled(proc,proclower,procmid-1,box);
  } else {
    box[2*dim] = cut;
    proc_box_tiled(proc,procmid,procupper,box);
  }
}
