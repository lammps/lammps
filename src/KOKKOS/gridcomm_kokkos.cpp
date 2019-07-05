/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gridcomm_kokkos.h"
#include <mpi.h>
#include "comm.h"
#include "kspace.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos_base.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

#define SWAPDELTA 8

/* ---------------------------------------------------------------------- */

template<class DeviceType>
GridCommKokkos<DeviceType>::GridCommKokkos(LAMMPS *lmp, MPI_Comm gcomm, int forward, int reverse,
                   int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
                   int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi,
                   int pxlo, int pxhi, int pylo, int pyhi, int pzlo, int pzhi)
  : Pointers(lmp)
{
  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);

  nforward = forward;
  nreverse = reverse;

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

  outxlo_max = oxlo;
  outxhi_max = oxhi;
  outylo_max = oylo;
  outyhi_max = oyhi;
  outzlo_max = ozlo;
  outzhi_max = ozhi;

  procxlo = pxlo;
  procxhi = pxhi;
  procylo = pylo;
  procyhi = pyhi;
  proczlo = pzlo;
  proczhi = pzhi;

  nswap = 0;
  swap = NULL;
  //buf1 = buf2 = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
GridCommKokkos<DeviceType>::GridCommKokkos(LAMMPS *lmp, MPI_Comm gcomm, int forward, int reverse,
                   int ixlo, int ixhi, int iylo, int iyhi, int izlo, int izhi,
                   int oxlo, int oxhi, int oylo, int oyhi, int ozlo, int ozhi,
                   int oxlo_max, int oxhi_max, int oylo_max, int oyhi_max,
                   int ozlo_max, int ozhi_max,
                   int pxlo, int pxhi, int pylo, int pyhi, int pzlo, int pzhi)
  : Pointers(lmp)
{
  gridcomm = gcomm;
  MPI_Comm_rank(gridcomm,&me);

  nforward = forward;
  nreverse = reverse;

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

  outxlo_max = oxlo_max;
  outxhi_max = oxhi_max;
  outylo_max = oylo_max;
  outyhi_max = oyhi_max;
  outzlo_max = ozlo_max;
  outzhi_max = ozhi_max;

  procxlo = pxlo;
  procxhi = pxhi;
  procylo = pylo;
  procyhi = pyhi;
  proczlo = pzlo;
  proczhi = pzhi;

  nswap = 0;
  swap = NULL;
  //buf1 = buf2 = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
GridCommKokkos<DeviceType>::~GridCommKokkos()
{
  for (int i = 0; i < nswap; i++) {
    //memoryKK->destroy_kokkos(swap[i].k_packlist,swap[i].packlist);
    //memoryKK->destroy_kokkos(swap[i].k_unpacklist,swap[i].unpacklist);
  }
  memory->sfree(swap);

  //memory->destroy(buf1);
  //memory->destroy(buf2);
}

/* ----------------------------------------------------------------------
   notify 6 neighbor procs how many ghost grid planes I need from them
   ghostxlo = # of lower grid planes I own that are needed from me
              by procxlo to become its upper ghost planes
   ghostxhi = # of upper grid planes I own that are needed from me
              by procxhi to become its lower ghost planes
   if no neighbor proc, value is from self
------------------------------------------------------------------------- */

template<class DeviceType>
void GridCommKokkos<DeviceType>::ghost_notify()
{
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
}

/* ----------------------------------------------------------------------
   check if all ghost grid comm needs overlap into non nearest-neighbor proc
   if yes, return 1, else return 0
------------------------------------------------------------------------- */

template<class DeviceType>
int GridCommKokkos<DeviceType>::ghost_overlap()
{
  int nearest = 0;
  if (ghostxlo > inxhi-inxlo+1) nearest = 1;
  if (ghostxhi > inxhi-inxlo+1) nearest = 1;
  if (ghostylo > inyhi-inylo+1) nearest = 1;
  if (ghostyhi > inyhi-inylo+1) nearest = 1;
  if (ghostzlo > inzhi-inzlo+1) nearest = 1;
  if (ghostzhi > inzhi-inzlo+1) nearest = 1;

  int nearest_all;
  MPI_Allreduce(&nearest,&nearest_all,1,MPI_INT,MPI_MIN,gridcomm);

  return nearest_all;
}

/* ----------------------------------------------------------------------
   create swap stencil for grid own/ghost communication
   swaps covers all 3 dimensions and both directions
   swaps cover multiple iterations in a direction if need grid pts
     from further away than nearest-neighbor proc
   same swap list used by forward and reverse communication
------------------------------------------------------------------------- */

template<class DeviceType>
void GridCommKokkos<DeviceType>::setup()
{
  int nsent,sendfirst,sendlast,recvfirst,recvlast;
  int sendplanes,recvplanes;
  int notdoneme,notdone;

  int maxswap = 6;
  swap = (Swap *) memory->smalloc(maxswap*sizeof(Swap),"Commgrid:swap");
  k_packlist = DAT::tdual_int_2d("Commgrid:packlist",maxswap,1);
  k_unpacklist = DAT::tdual_int_2d("Commgrid:unpacklist",maxswap,1);
  nswap = 0;

  // send own grid pts to -x processor, recv ghost grid pts from +x processor

  nsent = 0;
  sendfirst = inxlo;
  sendlast = inxhi;
  recvfirst = inxhi+1;
  notdone = 1;

  while (notdone) {
    if (nswap == maxswap) {
      maxswap += SWAPDELTA;
      swap = (Swap *)
        memory->srealloc(swap,maxswap*sizeof(Swap),"Commgrid:swap");
      k_packlist.resize(maxswap,k_packlist.extent(1));
      k_unpacklist.resize(maxswap,k_unpacklist.extent(1));
    }

    swap[nswap].sendproc = procxlo;
    swap[nswap].recvproc = procxhi;
    sendplanes = MIN(sendlast-sendfirst+1,ghostxlo-nsent);
    swap[nswap].npack =
      indices(k_packlist,nswap,
              sendfirst,sendfirst+sendplanes-1,inylo,inyhi,inzlo,inzhi);

    if (procxlo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procxlo,0,
                   &recvplanes,1,MPI_INT,procxhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(k_unpacklist,nswap,
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
    if (nswap == maxswap) {
      maxswap += 1;
      swap = (Swap *)
        memory->srealloc(swap,maxswap*sizeof(Swap),"Commgrid:swap");
      k_packlist.resize(maxswap,k_packlist.extent(1));
      k_unpacklist.resize(maxswap,k_unpacklist.extent(1));
    }

    swap[nswap].sendproc = procxhi;
    swap[nswap].recvproc = procxlo;
    sendplanes = MIN(sendlast-sendfirst+1,ghostxhi-nsent);
    swap[nswap].npack =
      indices(k_packlist,nswap,
              sendlast-sendplanes+1,sendlast,inylo,inyhi,inzlo,inzhi);

    if (procxhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procxhi,0,
                   &recvplanes,1,MPI_INT,procxlo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(k_unpacklist,nswap,
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
    if (nswap == maxswap) {
      maxswap += SWAPDELTA;
      swap = (Swap *)
        memory->srealloc(swap,maxswap*sizeof(Swap),"Commgrid:swap");
      k_packlist.resize(maxswap,k_packlist.extent(1));
      k_unpacklist.resize(maxswap,k_unpacklist.extent(1));
    }

    swap[nswap].sendproc = procylo;
    swap[nswap].recvproc = procyhi;
    sendplanes = MIN(sendlast-sendfirst+1,ghostylo-nsent);
    swap[nswap].npack =
      indices(k_packlist,nswap,
              outxlo,outxhi,sendfirst,sendfirst+sendplanes-1,inzlo,inzhi);

    if (procylo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procylo,0,
                   &recvplanes,1,MPI_INT,procyhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(k_unpacklist,nswap,
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
    if (nswap == maxswap) {
      maxswap += 1;
      swap = (Swap *)
        memory->srealloc(swap,maxswap*sizeof(Swap),"Commgrid:swap");
      k_packlist.resize(maxswap,k_packlist.extent(1));
      k_unpacklist.resize(maxswap,k_unpacklist.extent(1));
    }

    swap[nswap].sendproc = procyhi;
    swap[nswap].recvproc = procylo;
    sendplanes = MIN(sendlast-sendfirst+1,ghostyhi-nsent);
    swap[nswap].npack =
      indices(k_packlist,nswap,
              outxlo,outxhi,sendlast-sendplanes+1,sendlast,inzlo,inzhi);

    if (procyhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,procyhi,0,
                   &recvplanes,1,MPI_INT,procylo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(k_unpacklist,nswap,
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
    if (nswap == maxswap) {
      maxswap += SWAPDELTA;
      swap = (Swap *)
        memory->srealloc(swap,maxswap*sizeof(Swap),"Commgrid:swap");
      k_packlist.resize(maxswap,k_packlist.extent(1));
      k_unpacklist.resize(maxswap,k_unpacklist.extent(1));
    }

    swap[nswap].sendproc = proczlo;
    swap[nswap].recvproc = proczhi;
    sendplanes = MIN(sendlast-sendfirst+1,ghostzlo-nsent);
    swap[nswap].npack =
      indices(k_packlist,nswap,
              outxlo,outxhi,outylo,outyhi,sendfirst,sendfirst+sendplanes-1);

    if (proczlo != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,proczlo,0,
                   &recvplanes,1,MPI_INT,proczhi,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(k_unpacklist,nswap,
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
    if (nswap == maxswap) {
      maxswap += 1;
      swap = (Swap *)
        memory->srealloc(swap,maxswap*sizeof(Swap),"Commgrid:swap");
      k_packlist.resize(maxswap,k_packlist.extent(1));
      k_unpacklist.resize(maxswap,k_unpacklist.extent(1));
    }

    swap[nswap].sendproc = proczhi;
    swap[nswap].recvproc = proczlo;
    sendplanes = MIN(sendlast-sendfirst+1,ghostzhi-nsent);
    swap[nswap].npack =
      indices(k_packlist,nswap,
              outxlo,outxhi,outylo,outyhi,sendlast-sendplanes+1,sendlast);

    if (proczhi != me)
      MPI_Sendrecv(&sendplanes,1,MPI_INT,proczhi,0,
                   &recvplanes,1,MPI_INT,proczlo,0,gridcomm,MPI_STATUS_IGNORE);
    else recvplanes = sendplanes;

    swap[nswap].nunpack =
      indices(k_unpacklist,nswap,
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

  // nbuf = max of any forward/reverse pack/unpack

  nbuf = 0;
  for (int i = 0; i < nswap; i++) {
    nbuf = MAX(nbuf,swap[i].npack);
    nbuf = MAX(nbuf,swap[i].nunpack);
  }
  nbuf *= MAX(nforward,nreverse);
  //memory->create(buf1,nbuf,"Commgrid:buf1");
  k_buf1 = DAT::tdual_FFT_SCALAR_1d("Commgrid:buf1",nbuf);
  //memory->create(buf2,nbuf,"Commgrid:buf2");
  k_buf2 = DAT::tdual_FFT_SCALAR_1d("Commgrid:buf2",nbuf);
}

/* ----------------------------------------------------------------------
   use swap list in forward order to acquire copy of all needed ghost grid pts
------------------------------------------------------------------------- */

template<class DeviceType>
void GridCommKokkos<DeviceType>::forward_comm(KSpace *kspace, int which)
{
  k_packlist.sync<DeviceType>();
  k_unpacklist.sync<DeviceType>();

  KokkosBase* kspaceKKBase = dynamic_cast<KokkosBase*>(kspace);

  for (int m = 0; m < nswap; m++) {
    if (swap[m].sendproc == me)
      kspaceKKBase->pack_forward_kspace_kokkos(which,k_buf2,swap[m].npack,k_packlist,m);
    else
      kspaceKKBase->pack_forward_kspace_kokkos(which,k_buf1,swap[m].npack,k_packlist,m);
    DeviceType::fence();

    if (swap[m].sendproc != me) {
      FFT_SCALAR* buf1;
      FFT_SCALAR* buf2;
      if (lmp->kokkos->gpu_direct_flag) {
        buf1 = k_buf1.view<DeviceType>().data();
        buf2 = k_buf2.view<DeviceType>().data();
      } else {
        k_buf1.modify<DeviceType>();
        k_buf1.sync<LMPHostType>();
        buf1 = k_buf1.h_view.data();
        buf2 = k_buf2.h_view.data();
      }

      MPI_Irecv(buf2,nforward*swap[m].nunpack,MPI_FFT_SCALAR,
                swap[m].recvproc,0,gridcomm,&request);
      MPI_Send(buf1,nforward*swap[m].npack,MPI_FFT_SCALAR,
               swap[m].sendproc,0,gridcomm);
      MPI_Wait(&request,MPI_STATUS_IGNORE);

      if (!lmp->kokkos->gpu_direct_flag) {
        k_buf2.modify<LMPHostType>();
        k_buf2.sync<DeviceType>();
      }
    }

    kspaceKKBase->unpack_forward_kspace_kokkos(which,k_buf2,swap[m].nunpack,k_unpacklist,m);
    DeviceType::fence();
  }
}

/* ----------------------------------------------------------------------
   use swap list in reverse order to compute fully summed value
   for each owned grid pt that some other proc has copy of as a ghost grid pt
------------------------------------------------------------------------- */

template<class DeviceType>
void GridCommKokkos<DeviceType>::reverse_comm(KSpace *kspace, int which)
{
  k_packlist.sync<DeviceType>();
  k_unpacklist.sync<DeviceType>();

  KokkosBase* kspaceKKBase = dynamic_cast<KokkosBase*>(kspace);

  for (int m = nswap-1; m >= 0; m--) {
    if (swap[m].recvproc == me)
      kspaceKKBase->pack_reverse_kspace_kokkos(which,k_buf2,swap[m].nunpack,k_unpacklist,m);
    else
      kspaceKKBase->pack_reverse_kspace_kokkos(which,k_buf1,swap[m].nunpack,k_unpacklist,m);
    DeviceType::fence();

    if (swap[m].recvproc != me) {
      FFT_SCALAR* buf1;
      FFT_SCALAR* buf2;
      if (lmp->kokkos->gpu_direct_flag) {
        buf1 = k_buf1.view<DeviceType>().data();
        buf2 = k_buf2.view<DeviceType>().data();
      } else {
        k_buf1.modify<DeviceType>();
        k_buf1.sync<LMPHostType>();
        buf1 = k_buf1.h_view.data();
        buf2 = k_buf2.h_view.data();
      }

      MPI_Irecv(buf2,nreverse*swap[m].npack,MPI_FFT_SCALAR,
                swap[m].sendproc,0,gridcomm,&request);
      MPI_Send(buf1,nreverse*swap[m].nunpack,MPI_FFT_SCALAR,
               swap[m].recvproc,0,gridcomm);
      MPI_Wait(&request,MPI_STATUS_IGNORE);

      if (!lmp->kokkos->gpu_direct_flag) {
        k_buf2.modify<LMPHostType>();
        k_buf2.sync<DeviceType>();
      }
    }

    kspaceKKBase->unpack_reverse_kspace_kokkos(which,k_buf2,swap[m].npack,k_packlist,m);
    DeviceType::fence();
  }
}

/* ----------------------------------------------------------------------
   create 1d list of offsets into 3d array section (xlo:xhi,ylo:yhi,zlo:zhi)
   assume 3d array is allocated as (0:outxhi_max-outxlo_max+1,0:outyhi_max-outylo_max+1,
     0:outzhi_max-outzlo_max+1)
------------------------------------------------------------------------- */

template<class DeviceType>
int GridCommKokkos<DeviceType>::indices(DAT::tdual_int_2d &k_list, int index,
                       int xlo, int xhi, int ylo, int yhi, int zlo, int zhi)
{
  int nmax = (xhi-xlo+1) * (yhi-ylo+1) * (zhi-zlo+1);
  if (k_list.extent(1) < nmax)
    k_list.resize(k_list.extent(0),nmax);

  int nx = (outxhi_max-outxlo_max+1);
  int ny = (outyhi_max-outylo_max+1);

  k_list.sync<LMPHostType>();

  int n = 0;
  int ix,iy,iz;
  for (iz = zlo; iz <= zhi; iz++)
    for (iy = ylo; iy <= yhi; iy++)
      for (ix = xlo; ix <= xhi; ix++)
        k_list.h_view(index,n++) = (iz-outzlo_max)*ny*nx + (iy-outylo_max)*nx + (ix-outxlo_max);

  k_list.modify<LMPHostType>();
  k_list.sync<DeviceType>();

  return nmax;
}


/* ----------------------------------------------------------------------
   memory usage of send/recv bufs
------------------------------------------------------------------------- */

template<class DeviceType>
double GridCommKokkos<DeviceType>::memory_usage()
{
  double bytes = 2*nbuf * sizeof(double);
  return bytes;
}

namespace LAMMPS_NS {
template class GridCommKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class GridCommKokkos<LMPHostType>;
#endif
}
