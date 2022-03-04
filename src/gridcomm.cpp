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

#include "gridcomm.h"

#include "comm.h"
#include "error.h"
#include "irregular.h"
#include "kspace.h"
#include "fix.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum{REGULAR,TILED};

#define DELTA 16

/* ----------------------------------------------------------------------
   NOTES
   tiled implementation only currently works for RCB, not general tiled
   b/c RCB tree is used to find neighboring tiles
   if o indices for ghosts are < 0 or hi indices are >= N,
     then grid is treated as periodic in that dimension,
     communication is done across the periodic boundaries
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   constructor called by all classes except MSM
   gcomm = world communicator
   gn xyz = size of global grid
   i xyz lohi = portion of global grid this proc owns, 0 <= index < N
   o xyz lohi = owned grid portion + ghost grid cells needed in all directions
   if o indices are < 0 or hi indices are >= N,
     then grid is treated as periodic in that dimension,
     communication is done across the periodic boundaries
------------------------------------------------------------------------- */

GridComm::GridComm(LAMMPS *lmp, MPI_Comm gcomm,
                   int gnx, int gny, int gnz,
                   int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
                   int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi)
  : Pointers(lmp)
{
  if (comm->layout == Comm::LAYOUT_TILED) layout = TILED;
  else layout = REGULAR;

  if (layout == REGULAR) {
    int (*procneigh)[2] = comm->procneigh;
    initialize(gcomm,gnx,gny,gnz,
               ixlo,ixhi,iylo,iyhi,izlo,izhi,
               oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
               oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
               procneigh[0][0],procneigh[0][1],
               procneigh[1][0],procneigh[1][1],
               procneigh[2][0],procneigh[2][1]);
  } else {
    initialize(gcomm,gnx,gny,gnz,
               ixlo,ixhi,iylo,iyhi,izlo,izhi,
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
   gn xyz = size of global grid
   i xyz lohi = portion of global grid this proc owns, 0 <= index < N
   o xyz lohi = owned grid portion + ghost grid cells needed in all directions
   e xyz lohi for flag = 1: extent of larger grid stored by caller
   e xyz lohi for flag = 2: 6 neighbor procs
------------------------------------------------------------------------- */

GridComm::GridComm(LAMMPS *lmp, MPI_Comm gcomm, int flag,
                   int gnx, int gny, int gnz,
                   int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
                   int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi,
                   int exlo, int exhi, int eylo, int eyhi, int ezlo, int ezhi)
  : Pointers(lmp)
{
  if (comm->layout == Comm::LAYOUT_TILED) layout = TILED;
  else layout = REGULAR;

  if (flag == 1) {
    if (layout == REGULAR) {
      // this assumes gcomm = world
      int (*procneigh)[2] = comm->procneigh;
      initialize(gcomm,gnx,gny,gnz,
                 ixlo,ixhi,iylo,iyhi,izlo,izhi,
                 oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
                 exlo,exhi,eylo,eyhi,ezlo,ezhi,
                 procneigh[0][0],procneigh[0][1],
                 procneigh[1][0],procneigh[1][1],
                 procneigh[2][0],procneigh[2][1]);
    } else {
      initialize(gcomm,gnx,gny,gnz,
                 ixlo,ixhi,iylo,iyhi,izlo,izhi,
                 oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
                 exlo,exhi,eylo,eyhi,ezlo,ezhi,
                 0,0,0,0,0,0);
    }

  } else if (flag == 2) {
    if (layout == REGULAR) {
      initialize(gcomm,gnx,gny,gnz,
                 ixlo,ixhi,iylo,iyhi,izlo,izhi,
                 oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
                 oxlo,oxhi,oylo,oyhi,ozlo,ozhi,
                 exlo,exhi,eylo,eyhi,ezlo,ezhi);
    } else {
      error->all(FLERR,"GridComm does not support tiled layout with neighbor procs");
    }
  }
}

/* ---------------------------------------------------------------------- */

GridComm::~GridComm()
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

void GridComm::initialize(MPI_Comm gcomm,
                          int gnx, int gny, int gnz,
                          int ixlo, int ixhi, int iylo, int iyhi,
                          int izlo, int izhi,
                          int oxlo, int oxhi, int oylo, int oyhi,
                          int ozlo, int ozhi,
                          int fxlo, int fxhi, int fylo, int fyhi,
                          int fzlo, int fzhi,
                          int pxlo, int pxhi, int pylo, int pyhi,
                          int pzlo, int pzhi)
{
  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);
  MPI_Comm_size(gridcomm,&nprocs);

  nx = gnx;
  ny = gny;
  nz = gnz;

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

  // for REGULAR layout, proc xyz lohi = my 6 neighbor procs in this MPI_Comm

  if (layout == REGULAR) {
    procxlo = pxlo;
    procxhi = pxhi;
    procylo = pylo;
    procyhi = pyhi;
    proczlo = pzlo;
    proczhi = pzhi;
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

void GridComm::setup(int &nbuf1, int &nbuf2)
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

void GridComm::setup_regular(int &nbuf1, int &nbuf2)
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
   setup comm for RCB tiled proc domains
   each proc has arbitrary # of neighbors that overlap its ghost extent
   identify which procs will send me ghost cells, and vice versa
   may not be symmetric if both procs do not need same layers of ghosts
   comm pattern = post recvs for all my ghosts, send my owned, wait on recvs
     no exchanges by dimension, unlike CommTiled forward/reverse comm of particles
------------------------------------------------------------------------- */

void GridComm::setup_tiled(int &nbuf1, int &nbuf2)
{
  int i,m;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int ghostbox[6],pbc[3];

  // setup RCB tree of cut info for grid
  // access CommTiled to get cut dimension
  // cut = this proc's inlo in that dim
  // dim is -1 for proc 0, but never accessed

  rcbinfo = (RCBinfo *)
    memory->smalloc(nprocs*sizeof(RCBinfo),"GridComm:rcbinfo");
  RCBinfo rcbone;
  rcbone.dim = comm->rcbcutdim;
  if (rcbone.dim <= 0) rcbone.cut = inxlo;
  else if (rcbone.dim == 1) rcbone.cut = inylo;
  else if (rcbone.dim == 2) rcbone.cut = inzlo;
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
  ghostbox[4] = outzlo;
  ghostbox[5] = outzhi;

  pbc[0] = pbc[1] = pbc[2] = 0;

  memory->create(overlap_procs,nprocs,"GridComm:overlap_procs");
  noverlap = maxoverlap = 0;
  overlap = nullptr;

  ghost_box_drop(ghostbox,pbc);

  // send each proc an overlap message
  // content: me, index of my overlap, box that overlaps with its owned cells
  // ncopy = # of overlaps with myself, across a periodic boundary

  int *proclist;
  memory->create(proclist,noverlap,"GridComm:proclist");
  srequest = (Request *)
    memory->smalloc(noverlap*sizeof(Request),"GridComm:srequest");

  int nsend_request = 0;
  ncopy = 0;

  for (m = 0; m < noverlap; m++) {
    if (overlap[m].proc == me) ncopy++;
    else {
      proclist[nsend_request] = overlap[m].proc;
      srequest[nsend_request].sender = me;
      srequest[nsend_request].index = m;
      for (i = 0; i < 6; i++)
        srequest[nsend_request].box[i] = overlap[m].box[i];
      nsend_request++;
    }
  }

  Irregular *irregular = new Irregular(lmp);
  int nrecv_request = irregular->create_data(nsend_request,proclist,1);
  Request *rrequest =
    (Request *) memory->smalloc(nrecv_request*sizeof(Request),"GridComm:rrequest");
  irregular->exchange_data((char *) srequest,sizeof(Request),(char *) rrequest);
  irregular->destroy_data();

  // compute overlaps between received ghost boxes and my owned box
  // overlap box used to setup my Send data struct and respond to requests

  send = (Send *) memory->smalloc(nrecv_request*sizeof(Send),"GridComm:send");
  sresponse = (Response *)
    memory->smalloc(nrecv_request*sizeof(Response),"GridComm:sresponse");
  memory->destroy(proclist);
  memory->create(proclist,nrecv_request,"GridComm:proclist");

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
  Response *rresponse =
    (Response *) memory->smalloc(nrecv_response*sizeof(Response),"GridComm:rresponse");
  irregular->exchange_data((char *) sresponse,sizeof(Response),(char *) rresponse);
  irregular->destroy_data();
  delete irregular;

  // process received responses
  // box used to setup my Recv data struct after unwrapping via PBC
  // adjacent = 0 if any box of ghost cells does not adjoin my owned cells

  recv = (Recv *) memory->smalloc(nrecv_response*sizeof(Recv),"GridComm:recv");
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

  copy = (Copy *) memory->smalloc(ncopy*sizeof(Copy),"GridComm:copy");

  ncopy = 0;
  for (m = 0; m < noverlap; m++) {
    if (overlap[m].proc != me) continue;
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
   box = 6 integers = (xlo,xhi,ylo,yhi,zlo,zhi)
     each lo/hi value may extend beyonw 0 to N-1 into another periodic image
   pbc = flags in each dim of which periodic image the caller box was in
   when a box straddles a periodic bounadry, split it in two
   when a box does not straddle, drop it down RCB tree
     add all the procs it overlaps with to Overlap list
------------------------------------------------------------------------- */

void GridComm::ghost_box_drop(int *box, int *pbc)
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
      for (i = 0; i < 6; i++) overlap[noverlap].box[i] = box[i];
      for (i = 0; i < 3; i++) overlap[noverlap].pbc[i] = pbc[i];
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
   box = 6 integers = (xlo,xhi,ylo,yhi,zlo,zhi)
     each lo/hi value ranges from 0 to N-1 in a dim, N = grid size in that dim
     box is guaranteed to be wholly within the global domain
   return Np = # of procs, plist = proc IDs
------------------------------------------------------------------------- */

void GridComm::box_drop_grid(int *box, int proclower, int procupper,
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

int GridComm::ghost_adjacent()
{
  if (layout == REGULAR) return ghost_adjacent_regular();
  return ghost_adjacent_tiled();
}

/* ----------------------------------------------------------------------
   adjacent = 0 if a proc's ghost xyz lohi values exceed its subdomain size
   return 0 if adjacent=0 for any proc, else 1
------------------------------------------------------------------------- */

int GridComm::ghost_adjacent_regular()
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

int GridComm::ghost_adjacent_tiled()
{
  int adjacent_all;
  MPI_Allreduce(&adjacent,&adjacent_all,1,MPI_INT,MPI_MIN,gridcomm);
  return adjacent_all;
}

/* ----------------------------------------------------------------------
   forward comm of my owned cells to other's ghost cells
------------------------------------------------------------------------- */

void GridComm::forward_comm(int caller, void *ptr, int nper, int nbyte, int which,
                            void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (layout == REGULAR) {
    if (caller == KSPACE)
      forward_comm_regular<KSpace>((KSpace *) ptr,nper,nbyte,which,
                                   buf1,buf2,datatype);
    else if (caller == FIX)
      forward_comm_regular<Fix>((Fix *) ptr,nper,nbyte,which,
                                buf1,buf2,datatype);
  } else {
    if (caller == KSPACE)
      forward_comm_tiled<KSpace>((KSpace *) ptr,nper,nbyte,which,
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
void GridComm::
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
void GridComm::
forward_comm_tiled(T *ptr, int nper, int nbyte, int which,
                   void *buf1, void *vbuf2, MPI_Datatype datatype)
{
  int i,m,offset;

  char *buf2 = (char *) vbuf2;

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

void GridComm::reverse_comm(int caller, void *ptr, int nper, int nbyte, int which,
                            void *buf1, void *buf2, MPI_Datatype datatype)
{
  if (layout == REGULAR) {
    if (caller == KSPACE)
      reverse_comm_regular<KSpace>((KSpace *) ptr,nper,nbyte,which,
                                   buf1,buf2,datatype);
    else if (caller == FIX)
      reverse_comm_regular<Fix>((Fix *) ptr,nper,nbyte,which,
                                buf1,buf2,datatype);
  } else {
    if (caller == KSPACE)
      reverse_comm_tiled<KSpace>((KSpace *) ptr,nper,nbyte,which,
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
void GridComm::
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
void GridComm::
reverse_comm_tiled(T *ptr, int nper, int nbyte, int which,
                   void *buf1, void *vbuf2, MPI_Datatype datatype)
{
  int i,m,offset;

  char *buf2 = (char *) vbuf2;

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
   gather global grid values to proc 0, one grid chunk at a time
   proc 0 pings each proc for its grid chunk
   pack/unpack operations are performed by caller via callbacks
   caller can decide whether to store chunks, output them, etc
------------------------------------------------------------------------- */

void GridComm::gather(int /*caller*/, void *ptr, int nper, int nbyte,
                      int which, void *buf, MPI_Datatype datatype)
{
  int me = comm->me;
  Fix *fptr = (Fix *) ptr;

  // maxsize = max grid data owned by any proc

  int mysize = (inxhi-inxlo+1) * (inyhi-inylo+1) * (inzhi-inzlo+1);
  mysize *= nper;
  int maxsize;
  MPI_Allreduce(&mysize,&maxsize,1,MPI_INT,MPI_MAX,world);

  // pack my data via callback to caller

  char *mybuf;
  if (me == 0) memory->create(mybuf,maxsize*nbyte,"GridComm:mybuf");
  else memory->create(mybuf,mysize*nbyte,"GridComm:mybuf");
  fptr->pack_gather_grid(which,mybuf);

  // ping each proc for its data
  // unpack into full buffer via callback to caller

  int xlo,xhi,ylo,yhi,zlo,zhi,tmp;
  int bounds[6];

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(mybuf,maxsize,datatype,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Recv(bounds,6,MPI_INT,iproc,0,world,&status);
        xlo = bounds[0];
        xhi = bounds[1];
        ylo = bounds[2];
        yhi = bounds[3];
        zlo = bounds[4];
        zhi = bounds[5];
      } else {
        xlo = inxlo;
        xhi = inxhi;
        ylo = inylo;
        yhi = inyhi;
        zlo = inzlo;
        zhi = inzhi;
      }

      fptr->unpack_gather_grid(which,mybuf,buf,xlo,xhi,ylo,yhi,zlo,zhi);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(mybuf,mysize,datatype,0,0,world);
    bounds[0] = inxlo;
    bounds[1] = inxhi;
    bounds[2] = inylo;
    bounds[3] = inyhi;
    bounds[4] = inzlo;
    bounds[5] = inzhi;
    MPI_Send(bounds,6,MPI_INT,0,0,world);
  }

  memory->destroy(mybuf);
}

/* ----------------------------------------------------------------------
   create swap stencil for grid own/ghost communication
   swaps covers all 3 dimensions and both directions
   swaps cover multiple iterations in a direction if need grid pts
     from further away than nearest-neighbor proc
   same swap list used by forward and reverse communication
------------------------------------------------------------------------- */

void GridComm::grow_swap()
{
  maxswap += DELTA;
  swap = (Swap *) memory->srealloc(swap,maxswap*sizeof(Swap),"GridComm:swap");
}

/* ----------------------------------------------------------------------
   create swap stencil for grid own/ghost communication
   swaps covers all 3 dimensions and both directions
   swaps cover multiple iterations in a direction if need grid pts
     from further away than nearest-neighbor proc
   same swap list used by forward and reverse communication
------------------------------------------------------------------------- */

void GridComm::grow_overlap()
{
  maxoverlap += DELTA;
  overlap = (Overlap *)
    memory->srealloc(overlap,maxoverlap*sizeof(Overlap),"GridComm:overlap");
}

/* ----------------------------------------------------------------------
   create 1d list of offsets into 3d array section (xlo:xhi,ylo:yhi,zlo:zhi)
   assume 3d array is allocated as
     (fullxlo:fullxhi,fullylo:fullyhi,fullzlo:fullzhi)
------------------------------------------------------------------------- */

int GridComm::indices(int *&list,
                       int xlo, int xhi, int ylo, int yhi, int zlo, int zhi)
{
  int nmax = (xhi-xlo+1) * (yhi-ylo+1) * (zhi-zlo+1);
  memory->create(list,nmax,"GridComm:indices");
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
