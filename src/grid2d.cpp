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

#include "grid2d.h"

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

enum{REGULAR,TILED};

#define DELTA 16

static constexpr int OFFSET = 16384;

/* ----------------------------------------------------------------------
   NOTES
   tiled implementation only currently works for RCB, not general tiled
   b/c RCB tree is used to find neighboring tiles
   if o indices for ghosts are < 0 or hi indices are >= N,
     then grid is treated as periodic in that dimension,
     communication is done across the periodic boundaries
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   constructor called by all classes except PPPM and MSM
   gcomm = world communicator
   gnx, gny = size of global grid
   maxdist = max distance outside of proc domain a particle will be
   extra = additional ghost grid pts needed in each dim, e.g. for stencil
   shift = 0.0 for grid pt in lower-left corner of grid cell, 0.5 for center
   return:
     i xy lohi = portion of global grid this proc owns, 0 <= index < N
     o xy lohi = owned + ghost grid cells needed in all directions
   for non-periodic dims, o indices will not be < 0 or >= N,
     since no grid communication is done across non-periodic boundaries
------------------------------------------------------------------------- */

Grid2d::Grid2d(LAMMPS *lmp, MPI_Comm gcomm,
               int gnx, int gny,
               double maxdist, int extra, double shift,
               int &ixlo, int &ixhi, int &iylo, int &iyhi,
               int &oxlo, int &oxhi, int &oylo, int &oyhi)
  : Pointers(lmp)
{
  // store commnicator and global grid size
  // set layout mode

  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;

  ngrid[0] = nx; ngrid[1] = ny;

  if (comm->layout == Comm::LAYOUT_TILED) layout = TILED;
  else layout = REGULAR;

  // partition global grid across procs
  // i xyz lo/hi = lower/upper bounds of global grid this proc owns
  // indices range from 0 to N-1 inclusive in each dim

  int tmp1,tmp2;
  comm->partition_grid(nx, ny, 1, 0.0, ixlo, ixhi, iylo, iyhi, tmp1, tmp2);

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

  // error check on size of grid stored by this proc

  bigint total = (bigint) (oxhi - oxlo + 1) * (oyhi - oylo + 1);
  if (total > MAXSMALLINT)
    error->one(FLERR, "Too many owned+ghost grid2d points");

  // store grid bounds and proc neighs

  if (layout == REGULAR) {
    int (*procneigh)[2] = comm->procneigh;
    store(ixlo,ixhi,iylo,iyhi,
          oxlo,oxhi,oylo,oyhi,
          oxlo,oxhi,oylo,oyhi,
          procneigh[0][0],procneigh[0][1],
          procneigh[1][0],procneigh[1][1]);
  } else {
    store(ixlo,ixhi,iylo,iyhi,
          oxlo,oxhi,oylo,oyhi,
          oxlo,oxhi,oylo,oyhi,
          0,0,0,0);
  }
}

/* ----------------------------------------------------------------------
   constructor called by PPPM classes
   gcomm = world communicator
   gnx, gny = size of global grid
   i xy lohi = portion of global grid this proc owns, 0 <= index < N
   o xy lohi = owned grid portion + ghost grid cells needed in all directions
   if o indices are < 0 or hi indices are >= N,
     then grid is treated as periodic in that dimension,
     communication is done across the periodic boundaries
------------------------------------------------------------------------- */

Grid2d::Grid2d(LAMMPS *lmp, MPI_Comm gcomm,
               int gnx, int gny,
               int ixlo, int ixhi, int iylo, int iyhi,
               int oxlo, int oxhi, int oylo, int oyhi)
  : Pointers(lmp)
{
  // store commnicator and global grid size
  // set layout mode

  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;

  ngrid[0] = nx; ngrid[1] = ny;

  if (comm->layout == Comm::LAYOUT_TILED) layout = TILED;
  else layout = REGULAR;

  // store grid bounds and proc neighs

  if (layout == REGULAR) {
    int (*procneigh)[2] = comm->procneigh;
    store(ixlo,ixhi,iylo,iyhi,
          oxlo,oxhi,oylo,oyhi,
          oxlo,oxhi,oylo,oyhi,
          procneigh[0][0],procneigh[0][1],
          procneigh[1][0],procneigh[1][1]);
  } else {
    store(ixlo,ixhi,iylo,iyhi,
          oxlo,oxhi,oylo,oyhi,
          oxlo,oxhi,oylo,oyhi,
          0,0,0,0);
  }
}

/* ----------------------------------------------------------------------
   constructor called by MSM
   gcomm = world communicator or sub-communicator for a hierarchical grid
   flag = 1 if e xy lohi values = larger grid stored by caller in gcomm = world
   flag = 2 if e xy lohi values = 6 neighbor procs in gcomm
   gnx, gny = size of global grid
   i xy lohi = portion of global grid this proc owns, 0 <= index < N
   o xy lohi = owned grid portion + ghost grid cells needed in all directions
   e xy lohi for flag = 1: extent of larger grid stored by caller
   e xy lohi for flag = 2: 4 neighbor procs
------------------------------------------------------------------------- */

Grid2d::Grid2d(LAMMPS *lmp, MPI_Comm gcomm, int flag,
               int gnx, int gny,
               int ixlo, int ixhi, int iylo, int iyhi,
               int oxlo, int oxhi, int oylo, int oyhi,
               int exlo, int exhi, int eylo, int eyhi)
  : Pointers(lmp)
{
  // store commnicator and global grid size
  // set layout mode

  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;

  ngrid[0] = nx; ngrid[1] = ny;

  if (comm->layout == Comm::LAYOUT_TILED) layout = TILED;
  else layout = REGULAR;

  // store grid bounds and proc neighs

  if (flag == 1) {
    if (layout == REGULAR) {
      // this assumes gcomm = world
      int (*procneigh)[2] = comm->procneigh;
      store(ixlo,ixhi,iylo,iyhi,
            oxlo,oxhi,oylo,oyhi,
            exlo,exhi,eylo,eyhi,
            procneigh[0][0],procneigh[0][1],
            procneigh[1][0],procneigh[1][1]);
    } else {
      store(ixlo,ixhi,iylo,iyhi,
            oxlo,oxhi,oylo,oyhi,
            exlo,exhi,eylo,eyhi,
            0,0,0,0);
    }

  } else if (flag == 2) {
    if (layout == REGULAR) {
      store(ixlo,ixhi,iylo,iyhi,
            oxlo,oxhi,oylo,oyhi,
            oxlo,oxhi,oylo,oyhi,
            exlo,exhi,eylo,eyhi);
    } else {
      error->all(FLERR,"Grid2d does not support tiled layout with neighbor procs");
    }
  }
}

/* ---------------------------------------------------------------------- */

Grid2d::~Grid2d()
{
  // regular comm data struct

  for (int i = 0; i < nswap; i++) {
    memory->destroy(swap[i].packlist);
    memory->destroy(swap[i].unpacklist);
  }
  memory->sfree(swap);

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
}

/* ----------------------------------------------------------------------
   store constructor args in local variables
------------------------------------------------------------------------- */

void Grid2d::store(int ixlo, int ixhi, int iylo, int iyhi,
                   int oxlo, int oxhi, int oylo, int oyhi,
                   int fxlo, int fxhi, int fylo, int fyhi,
                   int pxlo, int pxhi, int pylo, int pyhi)
{
  inxlo = ixlo;
  inxhi = ixhi;
  inylo = iylo;
  inyhi = iyhi;

  outxlo = oxlo;
  outxhi = oxhi;
  outylo = oylo;
  outyhi = oyhi;

  fullxlo = fxlo;
  fullxhi = fxhi;
  fullylo = fylo;
  fullyhi = fyhi;

  // for REGULAR layout, proc xy lohi = my 4 neighbor procs in this MPI_Comm

  if (layout == REGULAR) {
    procxlo = pxlo;
    procxhi = pxhi;
    procylo = pylo;
    procyhi = pyhi;
  }

  // internal data initializations

  nswap = maxswap = 0;
  swap = nullptr;

  nsend = nrecv = ncopy = 0;
  send = nullptr;
  recv = nullptr;
  copy = nullptr;
  requests = nullptr;
}

/* ---------------------------------------------------------------------- */

void Grid2d::get_size(int &nxgrid, int &nygrid)
{
  nxgrid = nx;
  nygrid = ny;
}

/* ---------------------------------------------------------------------- */

void Grid2d::get_bounds(int &xlo, int &xhi, int &ylo, int &yhi)
{
  xlo = inxlo;
  xhi = inxhi;
  ylo = inylo;
  yhi = inyhi;
}

/* ----------------------------------------------------------------------
   return sizes of two buffers needed for communication
   either on regular grid or procs or irregular tiling
   nbuf1 = largest pack or unpack in any Send or Recv or Copy
   nbuf2 = larget of sum of all packs or unpacks in Send or Recv
   for regular comm, nbuf1 = nbuf2
   for irregular comm, nbuf2 >= nbuf2
   nbuf1,nbuf2 are just count of grid points
     caller converts them to message size for grid data it stores
------------------------------------------------------------------------- */

void Grid2d::setup(int &nbuf1, int &nbuf2)
{
  if (layout == REGULAR) setup_regular(nbuf1,nbuf2);
  else setup_tiled(nbuf1,nbuf2);
}

/* ----------------------------------------------------------------------
   setup comm for a regular grid of procs
   each proc has 6 neighbors
   comm pattern = series of swaps with one of those 6 procs
   can be multiple swaps with same proc if ghost extent is large
   swap may not be symmetric if both procs do not need same layers of ghosts
   all procs perform same # of swaps in a direction, even if some don't need it
------------------------------------------------------------------------- */

void Grid2d::setup_regular(int &nbuf1, int &nbuf2)
{
  int nsent,sendfirst,sendlast,recvfirst,recvlast;
  int sendplanes,recvplanes;
  int notdoneme,notdone;

  // notify 6 neighbor procs how many ghost grid planes I need from them
  // ghost xy lo = # of my lower grid planes that proc xy lo needs as its ghosts
  // ghost xy hi = # of my upper grid planes that proc xy hi needs as its ghosts
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
              sendfirst,sendfirst+sendplanes-1,inylo,inyhi);

    if (procxlo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procxlo,0,
                   &recvplanes,1,MPI_INT,procxhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              recvfirst,recvfirst+recvplanes-1,inylo,inyhi);

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
              sendlast-sendplanes+1,sendlast,inylo,inyhi);

    if (procxhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procxhi,0,
                   &recvplanes,1,MPI_INT,procxlo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              recvlast-recvplanes+1,recvlast,inylo,inyhi);

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
              outxlo,outxhi,sendfirst,sendfirst+sendplanes-1);

    if (procylo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procylo,0,
                   &recvplanes,1,MPI_INT,procyhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              outxlo,outxhi,recvfirst,recvfirst+recvplanes-1);

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
              outxlo,outxhi,sendlast-sendplanes+1,sendlast);

    if (procyhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procyhi,0,
                   &recvplanes,1,MPI_INT,procylo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(swap[nswap].unpacklist,
              outxlo,outxhi,recvlast-recvplanes+1,recvlast);

    nsent += sendplanes;
    sendfirst -= recvplanes;
    sendlast -= sendplanes;
    recvlast -= recvplanes;
    nswap++;

    if (nsent < ghostyhi) notdoneme = 1;
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
   setup comm for RCB tiled proc domains
   each proc has arbitrary # of neighbors that overlap its ghost extent
   identify which procs will send me ghost cells, and vice versa
   may not be symmetric if both procs do not need same layers of ghosts
   comm pattern = post recvs for all my ghosts, send my owned, wait on recvs
     no exchanges by dimension, unlike CommTiled forward/reverse comm of particles
------------------------------------------------------------------------- */

void Grid2d::setup_tiled(int &nbuf1, int &nbuf2)
{
  int i,m;
  double xlo,xhi,ylo,yhi;
  int ghostbox[4],pbc[2];

  // setup RCB tree of cut info for grid
  // access CommTiled to get cut dimension
  // cut = this proc's inlo in that dim
  // dim is -1 for proc 0, but never accessed

  rcbinfo = (RCBinfo *)
    memory->smalloc(nprocs*sizeof(RCBinfo),"grid2d:rcbinfo");
  RCBinfo rcbone;
  rcbone.dim = comm->rcbcutdim;
  if (rcbone.dim <= 0) rcbone.cut = inxlo;
  else if (rcbone.dim == 1) rcbone.cut = inylo;
  MPI_Allgather(&rcbone,sizeof(RCBinfo),MPI_CHAR,
                rcbinfo,sizeof(RCBinfo),MPI_CHAR,gridcomm);

  // find overlaps of my extended ghost box with all other procs
  // accounts for crossings of periodic boundaries
  // noverlap = # of overlaps, including self
  // overlap = vector of overlap info using Overlap data struct

  ghostbox[0] = outxlo;
  ghostbox[1] = outxhi;
  ghostbox[2] = outylo;
  ghostbox[3] = outyhi;

  pbc[0] = pbc[1] = 0;

  memory->create(overlap_procs,nprocs,"grid2d:overlap_procs");
  noverlap = maxoverlap = 0;
  overlap = nullptr;

  ghost_box_drop(ghostbox,pbc);

  // send each proc an overlap message
  // content: me, index of my overlap, box that overlaps with its owned cells
  // ncopy = # of overlaps with myself, across a periodic boundary

  int *proclist;
  memory->create(proclist,noverlap,"grid2d:proclist");
  srequest = (Request *)
    memory->smalloc(noverlap*sizeof(Request),"grid2d:srequest");

  int nsend_request = 0;
  ncopy = 0;

  for (m = 0; m < noverlap; m++) {
    if (overlap[m].proc == me) ncopy++;
    else {
      proclist[nsend_request] = overlap[m].proc;
      srequest[nsend_request].sender = me;
      srequest[nsend_request].index = m;
      for (i = 0; i < 4; i++)
        srequest[nsend_request].box[i] = overlap[m].box[i];
      nsend_request++;
    }
  }

  auto irregular = new Irregular(lmp);
  int nrecv_request = irregular->create_data(nsend_request,proclist,1);
  auto rrequest = (Request *) memory->smalloc(nrecv_request*sizeof(Request),"grid2d:rrequest");
  irregular->exchange_data((char *) srequest,sizeof(Request),(char *) rrequest);
  irregular->destroy_data();

  // compute overlaps between received ghost boxes and my owned box
  // overlap box used to setup my Send data struct and respond to requests

  send = (Send *) memory->smalloc(nrecv_request*sizeof(Send),"grid2d:send");
  sresponse = (Response *) memory->smalloc(nrecv_request*sizeof(Response),"grid2d:sresponse");
  memory->destroy(proclist);
  memory->create(proclist,nrecv_request,"grid2d:proclist");

  for (m = 0; m < nrecv_request; m++) {
    send[m].proc = rrequest[m].sender;
    xlo = MAX(rrequest[m].box[0],inxlo);
    xhi = MIN(rrequest[m].box[1],inxhi);
    ylo = MAX(rrequest[m].box[2],inylo);
    yhi = MIN(rrequest[m].box[3],inyhi);
    send[m].npack = indices(send[m].packlist,xlo,xhi,ylo,yhi);

    proclist[m] = rrequest[m].sender;
    sresponse[m].index = rrequest[m].index;
    sresponse[m].box[0] = xlo;
    sresponse[m].box[1] = xhi;
    sresponse[m].box[2] = ylo;
    sresponse[m].box[3] = yhi;
  }

  nsend = nrecv_request;

  // reply to each Request message with a Response message
  // content: index for the overlap on requestor, overlap box on my owned grid

  int nsend_response = nrecv_request;
  int nrecv_response = irregular->create_data(nsend_response,proclist,1);
  auto rresponse = (Response *) memory->smalloc(nrecv_response*sizeof(Response),"grid2d:rresponse");
  irregular->exchange_data((char *) sresponse,sizeof(Response),(char *) rresponse);
  irregular->destroy_data();
  delete irregular;

  // process received responses
  // box used to setup my Recv data struct after unwrapping via PBC
  // adjacent = 0 if any box of ghost cells does not adjoin my owned cells

  recv = (Recv *) memory->smalloc(nrecv_response*sizeof(Recv),"grid2d:recv");
  adjacent = 1;

  for (i = 0; i < nrecv_response; i++) {
    m = rresponse[i].index;
    recv[i].proc = overlap[m].proc;
    xlo = rresponse[i].box[0] + overlap[m].pbc[0] * nx;
    xhi = rresponse[i].box[1] + overlap[m].pbc[0] * nx;
    ylo = rresponse[i].box[2] + overlap[m].pbc[1] * ny;
    yhi = rresponse[i].box[3] + overlap[m].pbc[1] * ny;
    recv[i].nunpack = indices(recv[i].unpacklist,xlo,xhi,ylo,yhi);

    if (xlo != inxhi+1 && xhi != inxlo-1 &&
        ylo != inyhi+1 && yhi != inylo-1) adjacent = 0;
  }

  nrecv = nrecv_response;

  // create Copy data struct from overlaps with self

  copy = (Copy *) memory->smalloc(ncopy*sizeof(Copy),"grid2d:copy");

  ncopy = 0;
  for (m = 0; m < noverlap; m++) {
    if (overlap[m].proc != me) continue;
    xlo = overlap[m].box[0];
    xhi = overlap[m].box[1];
    ylo = overlap[m].box[2];
    yhi = overlap[m].box[3];
    copy[ncopy].npack = indices(copy[ncopy].packlist,xlo,xhi,ylo,yhi);
    xlo = overlap[m].box[0] + overlap[m].pbc[0] * nx;
    xhi = overlap[m].box[1] + overlap[m].pbc[0] * nx;
    ylo = overlap[m].box[2] + overlap[m].pbc[1] * ny;
    yhi = overlap[m].box[3] + overlap[m].pbc[1] * ny;
    copy[ncopy].nunpack = indices(copy[ncopy].unpacklist,xlo,xhi,ylo,yhi);
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

  memory->sfree(rcbinfo);
  memory->destroy(proclist);
  memory->destroy(overlap_procs);
  memory->sfree(overlap);
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

/* ----------------------------------------------------------------------
   recursively split a box until it doesn't overlap any periodic boundaries
   box = 4 integers = (xlo,xhi,ylo,yhi)
     each lo/hi value may extend beyonw 0 to N-1 into another periodic image
   pbc = flags in each dim of which periodic image the caller box was in
   when a box straddles a periodic bounadry, split it in two
   when a box does not straddle, drop it down RCB tree
     add all the procs it overlaps with to Overlap list
------------------------------------------------------------------------- */

void Grid2d::ghost_box_drop(int *box, int *pbc)
{
  int i,m;

  // newbox12 and newpbc are initially copies of caller box and pbc

  int newbox1[4],newbox2[4],newpbc[2];

  for (i = 0; i < 4; i++) newbox1[i] = newbox2[i] = box[i];
  for (i = 0; i < 2; i++) newpbc[i] = pbc[i];

  // 4 if tests to see if box needs to be split across a periodic boundary
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

  // box is not split, drop on RCB tree
  // returns nprocs = # of procs it overlaps, including self
  // returns proc_overlap = list of proc IDs it overlaps
  // skip self overlap if no crossing of periodic boundaries
  // do not skip self if overlap is in another periodic image

  } else {
    splitflag = 0;
    int np = 0;
    box_drop_grid(box,0,nprocs-1,np,overlap_procs);
    for (m = 0; m < np; m++) {
      if (noverlap == maxoverlap) grow_overlap();
      if (overlap_procs[m] == me &&
          pbc[0] == 0 && pbc[1] == 0 && pbc[2] == 0) continue;
      overlap[noverlap].proc = overlap_procs[m];
      for (i = 0; i < 4; i++) overlap[noverlap].box[i] = box[i];
      for (i = 0; i < 2; i++) overlap[noverlap].pbc[i] = pbc[i];
      noverlap++;
    }
  }

  // recurse with 2 split boxes

  if (splitflag) {
    ghost_box_drop(newbox1,pbc);
    ghost_box_drop(newbox2,newpbc);
  }
}

/* ----------------------------------------------------------------------
   recursively drop a box down the RCB tree to find all procs it overlaps with
   box = 4 integers = (xlo,xhi,ylo,yhi)
     each lo/hi value ranges from 0 to N-1 in a dim, N = grid size in that dim
     box is guaranteed to be wholly within the global domain
   return Np = # of procs, plist = proc IDs
------------------------------------------------------------------------- */

void Grid2d::box_drop_grid(int *box, int proclower, int procupper,
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

/* ----------------------------------------------------------------------
   check if all procs only need ghost info from adjacent procs
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int Grid2d::ghost_adjacent()
{
  if (layout == REGULAR) return ghost_adjacent_regular();
  return ghost_adjacent_tiled();
}

/* ----------------------------------------------------------------------
   adjacent = 0 if a proc's ghost xy lohi values exceed its subdomain size
   return 0 if adjacent=0 for any proc, else 1
------------------------------------------------------------------------- */

int Grid2d::ghost_adjacent_regular()
{
  adjacent = 1;
  if (ghostxlo > inxhi-inxlo+1) adjacent = 0;
  if (ghostxhi > inxhi-inxlo+1) adjacent = 0;
  if (ghostylo > inyhi-inylo+1) adjacent = 0;
  if (ghostyhi > inyhi-inylo+1) adjacent = 0;

  int adjacent_all;
  MPI_Allreduce(&adjacent,&adjacent_all,1,MPI_INT,MPI_MIN,gridcomm);
  return adjacent_all;
}

/* ----------------------------------------------------------------------
   adjacent = 0 if a proc's received ghosts were flagged
     as non-adjacent in setup_tiled()
   return 0 if adjacent=0 for any proc, else 1
------------------------------------------------------------------------- */

int Grid2d::ghost_adjacent_tiled()
{
  int adjacent_all;
  MPI_Allreduce(&adjacent,&adjacent_all,1,MPI_INT,MPI_MIN,gridcomm);
  return adjacent_all;
}

/* ----------------------------------------------------------------------
   forward comm of my owned cells to other's ghost cells
------------------------------------------------------------------------- */

void Grid2d::forward_comm(int caller, void *ptr, int nper, int nbyte, int which,
                          void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (layout == REGULAR) {
    if (caller == KSPACE)
      forward_comm_regular<KSpace>((KSpace *) ptr,nper,nbyte,which,
                                   buf1,buf2,datatype);
    else if (caller == PAIR)
      forward_comm_regular<Pair>((Pair *) ptr,nper,nbyte,which,
                                   buf1,buf2,datatype);
    else if (caller == FIX)
      forward_comm_regular<Fix>((Fix *) ptr,nper,nbyte,which,
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
   forward comm on regular grid of procs via list of swaps with 6 neighbor procs
------------------------------------------------------------------------- */

template < class T >
void Grid2d::
forward_comm_regular(T *ptr, int nper, int /*nbyte*/, int which,
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
   forward comm on tiled grid decomp via Send/Recv lists of each neighbor proc
------------------------------------------------------------------------- */

template < class T >
void Grid2d::
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

void Grid2d::reverse_comm(int caller, void *ptr, int nper, int nbyte, int which,
                          void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (layout == REGULAR) {
    if (caller == KSPACE)
      reverse_comm_regular<KSpace>((KSpace *) ptr,nper,nbyte,which,
                                   buf1,buf2,datatype);
    else if (caller == PAIR)
      reverse_comm_regular<Pair>((Pair *) ptr,nper,nbyte,which,
                                   buf1,buf2,datatype);
    else if (caller == FIX)
      reverse_comm_regular<Fix>((Fix *) ptr,nper,nbyte,which,
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
   reverse comm on regular grid of procs via list of swaps with 6 neighbor procs
------------------------------------------------------------------------- */

template < class T >
void Grid2d::
reverse_comm_regular(T *ptr, int nper, int /*nbyte*/, int which,
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
   reverse comm on tiled grid decomp via Send/Recv lists of each neighbor proc
------------------------------------------------------------------------- */

template < class T >
void Grid2d::
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

/* ----------------------------------------------------------------------
   setup remap from old grid decomposition to this grid decomposition
   pack/unpack operations are performed by caller via callbacks
------------------------------------------------------------------------- */

void Grid2d::remap_setup(Grid2d *old, int &ngrid1_buf, int &ngrid2_buf)
{
  if (layout == REGULAR) remap_setup_regular(old,ngrid1_buf,ngrid2_buf);
  else remap_setup_tiled(old,ngrid1_buf,ngrid2_buf);
}

void Grid2d::remap_setup_regular(Grid2d *old, int &ngrid1_buf, int &ngrid2_buf)
{
  ngrid1_buf = 0;
  ngrid2_buf = 0;
}

void Grid2d::remap_setup_tiled(Grid2d *old, int &ngrid1_buf, int &ngrid2_buf)
{
  ngrid1_buf = 0;
  ngrid2_buf = 0;
}

/* ----------------------------------------------------------------------
   perform remap from old grid decomposition to this grid decomposition
   pack/unpack operations are performed by caller via callbacks
------------------------------------------------------------------------- */

void Grid2d::remap(int caller, void *ptr, int nper, int nbyte,
		   void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (caller == FIX) remap_style<Fix>((Fix *) ptr,nper,nbyte,buf1,buf2,datatype);
}

template < class T >
void Grid2d::remap_style(T *ptr, int nper, int nbyte,
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

/* ----------------------------------------------------------------------
   gather global grid values to proc 0, one grid chunk at a time
   proc 0 pings each proc for its grid chunk
   pack/unpack operations are performed by caller via callbacks
   caller can decide whether to store chunks, output them, etc
------------------------------------------------------------------------- */

void Grid2d::gather(int /*caller*/, void *ptr, int nper, int nbyte,
                      int which, void *buf, MPI_Datatype datatype)
{
  int me = comm->me;
  Fix *fptr = (Fix *) ptr;

  // maxsize = max grid data owned by any proc

  int mysize = (inxhi-inxlo+1) * (inyhi-inylo+1);
  mysize *= nper;
  int maxsize;
  MPI_Allreduce(&mysize,&maxsize,1,MPI_INT,MPI_MAX,world);

  // pack my data via callback to caller

  char *mybuf;
  if (me == 0) memory->create(mybuf,maxsize*nbyte,"grid2d:mybuf");
  else memory->create(mybuf,mysize*nbyte,"grid2d:mybuf");
  fptr->pack_gather_grid(which,mybuf);

  // ping each proc for its data
  // unpack into full buffer via callback to caller

  int xlo,xhi,ylo,yhi,tmp;
  int bounds[4];

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(mybuf,maxsize,datatype,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Recv(bounds,4,MPI_INT,iproc,0,world,&status);
        xlo = bounds[0];
        xhi = bounds[1];
        ylo = bounds[2];
        yhi = bounds[3];
      } else {
        xlo = inxlo;
        xhi = inxhi;
        ylo = inylo;
        yhi = inyhi;
      }

      fptr->unpack_gather_grid(which,mybuf,buf,xlo,xhi,ylo,yhi,0,0);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(mybuf,mysize,datatype,0,0,world);
    bounds[0] = inxlo;
    bounds[1] = inxhi;
    bounds[2] = inylo;
    bounds[3] = inyhi;
    MPI_Send(bounds,4,MPI_INT,0,0,world);
  }

  memory->destroy(mybuf);
}

/* ----------------------------------------------------------------------
   create swap stencil for grid own/ghost communication
   swaps covers all 2 dimensions and both directions
   swaps cover multiple iterations in a direction if need grid pts
     from further away than nearest-neighbor proc
   same swap list used by forward and reverse communication
------------------------------------------------------------------------- */

void Grid2d::grow_swap()
{
  maxswap += DELTA;
  swap = (Swap *) memory->srealloc(swap,maxswap*sizeof(Swap),"grid2d:swap");
}

/* ----------------------------------------------------------------------
   create swap stencil for grid own/ghost communication
   swaps covers all 3 dimensions and both directions
   swaps cover multiple iterations in a direction if need grid pts
     from further away than nearest-neighbor proc
   same swap list used by forward and reverse communication
------------------------------------------------------------------------- */

void Grid2d::grow_overlap()
{
  maxoverlap += DELTA;
  overlap = (Overlap *)
    memory->srealloc(overlap,maxoverlap*sizeof(Overlap),"grid2d:overlap");
}

/* ----------------------------------------------------------------------
   create 1d list of offsets into 2d array section (xlo:xhi,ylo:yhi)
   assume 2d array is allocated as
     (fullxlo:fullxhi,fullylo:fullyhi)
------------------------------------------------------------------------- */

int Grid2d::indices(int *&list, int xlo, int xhi, int ylo, int yhi)
{
  int nmax = (xhi-xlo+1) * (yhi-ylo+1);
  memory->create(list,nmax,"grid2d:indices");
  if (nmax == 0) return 0;

  int nx = (fullxhi-fullxlo+1);

  int n = 0;
  int ix,iy;
  for (iy = ylo; iy <= yhi; iy++)
    for (ix = xlo; ix <= xhi; ix++)
      list[n++] = (iy-fullylo)*nx + (ix-fullxlo);

  return nmax;
}
