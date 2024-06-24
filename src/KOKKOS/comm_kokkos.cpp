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

#include "comm_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "atom_vec_kokkos.h"
#include "compute.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "kokkos.h"
#include "kokkos_base.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "output.h"
#include "pair.h"

#include <Kokkos_Sort.hpp>

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 10000;
static constexpr int BUFEXTRA = 1000;

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space
------------------------------------------------------------------------- */

CommKokkos::CommKokkos(LAMMPS *lmp) : CommBrick(lmp)
{
  if (sendlist) for (int i = 0; i < maxswap; i++) memory->destroy(sendlist[i]);
  memory->sfree(sendlist);
  sendlist = nullptr;
  k_sendlist = DAT::tdual_int_2d();
  k_total_send = DAT::tdual_int_scalar("comm::k_total_send");

  // error check for disallow of OpenMP threads?

  // initialize comm buffers & exchange memory

  memory->destroy(buf_send);
  buf_send = nullptr;
  memory->destroy(buf_recv);
  buf_recv = nullptr;

  k_exchange_sendlist = DAT::tdual_int_1d("comm:k_exchange_sendlist",100);
  k_exchange_copylist = DAT::tdual_int_1d("comm:k_exchange_copylist",100);
  k_count = DAT::tdual_int_scalar("comm:k_count");

  memory->destroy(maxsendlist);
  maxsendlist = nullptr;
  memory->create(maxsendlist,maxswap,"comm:maxsendlist");
  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
  }
  memoryKK->create_kokkos(k_sendlist,sendlist,maxswap,BUFMIN,"comm:sendlist");

  max_buf_pair = 0;
  k_buf_send_pair = DAT::tdual_xfloat_1d("comm:k_buf_send_pair",1);
  k_buf_recv_pair = DAT::tdual_xfloat_1d("comm:k_recv_send_pair",1);

  max_buf_fix = 0;
  k_buf_send_fix = DAT::tdual_xfloat_1d("comm:k_buf_send_fix",1);
  k_buf_recv_fix = DAT::tdual_xfloat_1d("comm:k_recv_send_fix",1);
}

/* ---------------------------------------------------------------------- */

CommKokkos::~CommKokkos()
{
  memoryKK->destroy_kokkos(k_sendlist,sendlist);
  sendlist = nullptr;
  memoryKK->destroy_kokkos(k_buf_send,buf_send);
  buf_send = nullptr;
  memoryKK->destroy_kokkos(k_buf_recv,buf_recv);
  buf_recv = nullptr;
}

/* ---------------------------------------------------------------------- */

void CommKokkos::init()
{
  maxsend = BUFMIN;
  maxrecv = BUFMIN;

  grow_send_kokkos(maxsend+bufextra,0,Host);
  grow_recv_kokkos(maxrecv,Host);

  atomKK = (AtomKokkos *) atom;
  exchange_comm_classic = lmp->kokkos->exchange_comm_classic;
  forward_comm_classic = lmp->kokkos->forward_comm_classic;
  forward_pair_comm_classic = lmp->kokkos->forward_pair_comm_classic;
  reverse_pair_comm_classic = lmp->kokkos->reverse_pair_comm_classic;
  forward_fix_comm_classic = lmp->kokkos->forward_fix_comm_classic;
  reverse_comm_classic = lmp->kokkos->reverse_comm_classic;
  exchange_comm_on_host = lmp->kokkos->exchange_comm_on_host;
  forward_comm_on_host = lmp->kokkos->forward_comm_on_host;
  reverse_comm_on_host = lmp->kokkos->reverse_comm_on_host;

  CommBrick::init();

  int check_forward = 0;
  int check_reverse = 0;
  if (force->pair && (force->pair->execution_space == Host))
    check_forward += force->pair->comm_forward;
  if (force->pair && (force->pair->execution_space == Host))
    check_reverse += force->pair->comm_reverse;

  for (const auto &fix : modify->get_fix_list()) {
    check_forward += fix->comm_forward;
    check_reverse += fix->comm_reverse;
  }

  for (const auto &compute : modify->get_compute_list()) {
    check_forward += compute->comm_forward;
    check_reverse += compute->comm_reverse;
  }

  for (const auto &dump : output->get_dump_list()) {
    check_forward += dump->comm_forward;
    check_reverse += dump->comm_reverse;
  }

  if (force->newton == 0) check_reverse = 0;
  if (force->pair) check_reverse += force->pair->comm_reverse_off;

  if (!comm_f_only) {// not all Kokkos atom_vec styles have reverse pack/unpack routines yet
    reverse_comm_classic = true;
    lmp->kokkos->reverse_comm_classic = 1;
  }

  if (ghost_velocity && atomKK->avecKK->no_comm_vel_flag) { // not all Kokkos atom_vec styles have comm vel pack/unpack routines yet
    forward_comm_classic = true;
    lmp->kokkos->forward_comm_classic = 1;
  }
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommKokkos::forward_comm(int dummy)
{
  if (!forward_comm_classic) {
    if (forward_comm_on_host) forward_comm_device<LMPHostType>();
    else forward_comm_device<LMPDeviceType>();
    return;
  }

  k_sendlist.sync<LMPHostType>();

  if (comm_x_only) {
    atomKK->sync(Host,X_MASK);
    atomKK->modified(Host,X_MASK);
  } else if (ghost_velocity) {
    atomKK->sync(Host,X_MASK | V_MASK);
    atomKK->modified(Host,X_MASK | V_MASK);
  } else {
    atomKK->sync(Host,ALL_MASK);
    atomKK->modified(Host,ALL_MASK);
  }

  CommBrick::forward_comm(dummy);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::forward_comm_device()
{
  int n;
  MPI_Request request;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  k_sendlist.sync<DeviceType>();

  if (comm->nprocs == 1 && !ghost_velocity) {
    k_swap.sync<DeviceType>();
    k_swap2.sync<DeviceType>();
    k_pbc.sync<DeviceType>();
    n = atomKK->avecKK->pack_comm_self_fused(totalsend,k_sendlist,k_sendnum_scan,
                    k_firstrecv,k_pbc_flag,k_pbc,k_g2l);
  } else {

    for (int iswap = 0; iswap < nswap; iswap++) {
      if (sendproc[iswap] != me) {
        if (comm_x_only) {
          if (size_forward_recv[iswap]) {
            buf = atomKK->k_x.view<DeviceType>().data() +
              firstrecv[iswap]*atomKK->k_x.view<DeviceType>().extent(1);
            MPI_Irecv(buf,size_forward_recv[iswap],MPI_DOUBLE,
                      recvproc[iswap],0,world,&request);
          }
          auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_kokkos(sendnum[iswap],k_sendlist_iswap,
                                     k_buf_send,pbc_flag[iswap],pbc[iswap]);
          DeviceType().fence();
          if (n) {
            MPI_Send(k_buf_send.view<DeviceType>().data(),
                     n,MPI_DOUBLE,sendproc[iswap],0,world);
          }

          if (size_forward_recv[iswap])
            MPI_Wait(&request,MPI_STATUS_IGNORE);

        } else if (ghost_velocity) {
          if (size_forward_recv[iswap]) {
            MPI_Irecv(k_buf_recv.view<DeviceType>().data(),
                      size_forward_recv[iswap],MPI_DOUBLE,
                      recvproc[iswap],0,world,&request);
          }
          auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_vel_kokkos(sendnum[iswap],k_sendlist_iswap,
                                         k_buf_send,pbc_flag[iswap],pbc[iswap]);
          DeviceType().fence();
          if (n) {
            MPI_Send(k_buf_send.view<DeviceType>().data(),n,
                     MPI_DOUBLE,sendproc[iswap],0,world);
          }
          if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          atomKK->avecKK->unpack_comm_vel_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_recv);
          DeviceType().fence();
        } else {
          if (size_forward_recv[iswap])
            MPI_Irecv(k_buf_recv.view<DeviceType>().data(),
                      size_forward_recv[iswap],MPI_DOUBLE,
                      recvproc[iswap],0,world,&request);
          auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_kokkos(sendnum[iswap],k_sendlist_iswap,
                                     k_buf_send,pbc_flag[iswap],pbc[iswap]);
          DeviceType().fence();
          if (n)
            MPI_Send(k_buf_send.view<DeviceType>().data(),n,
                     MPI_DOUBLE,sendproc[iswap],0,world);
          if (size_forward_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
          atomKK->avecKK->unpack_comm_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_recv);
          DeviceType().fence();
        }
      } else {
        if (!ghost_velocity) {
          if (sendnum[iswap]) {
            auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
            n = atomKK->avecKK->pack_comm_self(sendnum[iswap],k_sendlist_iswap,
                                     firstrecv[iswap],pbc_flag[iswap],pbc[iswap]);
            DeviceType().fence();
          }
        } else {
          auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_vel_kokkos(sendnum[iswap],k_sendlist_iswap,
                                         k_buf_send,pbc_flag[iswap],pbc[iswap]);
          DeviceType().fence();
          atomKK->avecKK->unpack_comm_vel_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_send);
          DeviceType().fence();
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommKokkos::reverse_comm()
{
  if (!reverse_comm_classic) {
    if (reverse_comm_on_host) reverse_comm_device<LMPHostType>();
    else reverse_comm_device<LMPDeviceType>();
    return;
  }

  k_sendlist.sync<LMPHostType>();

  if (comm_f_only)
    atomKK->sync(Host,F_MASK);
  else
    atomKK->sync(Host,ALL_MASK);

  CommBrick::reverse_comm();

  if (comm_f_only)
    atomKK->modified(Host,F_MASK);
  else
    atomKK->modified(Host,ALL_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::reverse_comm_device()
{
  int n;
  MPI_Request request;
  double *buf;

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  k_sendlist.sync<DeviceType>();

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap] != me) {
      if (comm_f_only) {
        if (size_reverse_recv[iswap])
            MPI_Irecv(k_buf_recv.view<DeviceType>().data(),size_reverse_recv[iswap],MPI_DOUBLE,
                    sendproc[iswap],0,world,&request);
        if (size_reverse_send[iswap]) {
          buf = atomKK->k_f.view<DeviceType>().data() +
            firstrecv[iswap]*atomKK->k_f.view<DeviceType>().extent(1);

          MPI_Send(buf,size_reverse_send[iswap],MPI_DOUBLE,
                   recvproc[iswap],0,world);
        }
        if (size_reverse_recv[iswap])
          MPI_Wait(&request,MPI_STATUS_IGNORE);

      } else {
        if (size_reverse_recv[iswap])
          MPI_Irecv(k_buf_recv.view<DeviceType>().data(),
                    size_reverse_recv[iswap],MPI_DOUBLE,
                    sendproc[iswap],0,world,&request);
        n = atomKK->avecKK->pack_reverse_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_send);
        DeviceType().fence();
        if (n)
          MPI_Send(k_buf_send.view<DeviceType>().data(),n,
                   MPI_DOUBLE,recvproc[iswap],0,world);
        if (size_reverse_recv[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);
      }
      auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
      atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap],k_sendlist_iswap,
                                k_buf_recv);
      DeviceType().fence();
    } else {
      if (sendnum[iswap]) {
        auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
        n = atomKK->avecKK->pack_reverse_self(sendnum[iswap],k_sendlist_iswap,
                                 firstrecv[iswap]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommKokkos::forward_comm(Fix *fix, int size)
{
  if (fix->execution_space == Host || !fix->forward_comm_device || forward_fix_comm_classic) {
    k_sendlist.sync<LMPHostType>();
    CommBrick::forward_comm(fix,size);
  } else {
    k_sendlist.sync<LMPDeviceType>();
    forward_comm_device<LMPDeviceType>(fix,size);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::forward_comm_device(Fix *fix, int size)
{
  int iswap,n,nsize;
  MPI_Request request;
  DAT::tdual_xfloat_1d k_buf_tmp;

  if (size) nsize = size;
  else nsize = fix->comm_forward;
  KokkosBase* fixKKBase = dynamic_cast<KokkosBase*>(fix);

  for (iswap = 0; iswap < nswap; iswap++) {
    int n = MAX(max_buf_fix,nsize*sendnum[iswap]);
    n = MAX(n,nsize*recvnum[iswap]);
    if (n > max_buf_fix)
      grow_buf_fix(n);
  }

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
    n = fixKKBase->pack_forward_comm_kokkos(sendnum[iswap],k_sendlist_iswap,
                                      k_buf_send_fix,pbc_flag[iswap],pbc[iswap]);
    DeviceType().fence();

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      double* buf_send_fix;
      double* buf_recv_fix;
      if (lmp->kokkos->gpu_aware_flag) {
        buf_send_fix = k_buf_send_fix.view<DeviceType>().data();
        buf_recv_fix = k_buf_recv_fix.view<DeviceType>().data();
      } else {
        k_buf_send_fix.modify<DeviceType>();
        k_buf_send_fix.sync<LMPHostType>();
        buf_send_fix = k_buf_send_fix.h_view.data();
        buf_recv_fix = k_buf_recv_fix.h_view.data();
      }

      if (recvnum[iswap]) {
        MPI_Irecv(buf_recv_fix,nsize*recvnum[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);
      }
      if (sendnum[iswap])
        MPI_Send(buf_send_fix,n,MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);

      if (!lmp->kokkos->gpu_aware_flag) {
        k_buf_recv_fix.modify<LMPHostType>();
        k_buf_recv_fix.sync<DeviceType>();
      }
      k_buf_tmp = k_buf_recv_fix;
    } else k_buf_tmp = k_buf_send_fix;

    // unpack buffer

    fixKKBase->unpack_forward_comm_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_tmp);
    DeviceType().fence();
  }
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommKokkos::reverse_comm(Fix *fix, int size)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::reverse_comm(fix, size);
}


/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for pack size to ensure buf_send is big enough
   handshake sizes before each Irecv/Send to ensure buf_recv is big enough
------------------------------------------------------------------------- */

void CommKokkos::reverse_comm_variable(Fix *fix)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::reverse_comm_variable(fix);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommKokkos::forward_comm(Compute *compute)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::forward_comm(compute);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommKokkos::forward_comm(Bond *bond)
{
  CommBrick::forward_comm(bond);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommKokkos::reverse_comm(Bond *bond)
{
  CommBrick::reverse_comm(bond);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommKokkos::reverse_comm(Compute *compute)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::reverse_comm(compute);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommKokkos::forward_comm(Pair *pair)
{
  if (pair->execution_space == Host || forward_pair_comm_classic) {
    k_sendlist.sync<LMPHostType>();
    CommBrick::forward_comm(pair);
  } else {
    k_sendlist.sync<LMPDeviceType>();
    forward_comm_device<LMPDeviceType>(pair);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::forward_comm_device(Pair *pair)
{
  int iswap,n;
  MPI_Request request;
  DAT::tdual_xfloat_1d k_buf_tmp;

  int nsize = pair->comm_forward;
  KokkosBase* pairKKBase = dynamic_cast<KokkosBase*>(pair);

  int nmax = max_buf_pair;
  for (iswap = 0; iswap < nswap; iswap++) {
    nmax = MAX(nmax,nsize*sendnum[iswap]);
    nmax = MAX(nmax,nsize*recvnum[iswap]);
  }
  if (nmax > max_buf_pair)
    grow_buf_pair(nmax);

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
    n = pairKKBase->pack_forward_comm_kokkos(sendnum[iswap],k_sendlist_iswap,
                                       k_buf_send_pair,pbc_flag[iswap],pbc[iswap]);
    DeviceType().fence();

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      double* buf_send_pair;
      double* buf_recv_pair;
      if (lmp->kokkos->gpu_aware_flag) {
        buf_send_pair = k_buf_send_pair.view<DeviceType>().data();
        buf_recv_pair = k_buf_recv_pair.view<DeviceType>().data();
      } else {
        k_buf_send_pair.modify<DeviceType>();
        k_buf_send_pair.sync<LMPHostType>();
        buf_send_pair = k_buf_send_pair.h_view.data();
        buf_recv_pair = k_buf_recv_pair.h_view.data();
      }

      if (recvnum[iswap]) {
        MPI_Irecv(buf_recv_pair,nsize*recvnum[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);
      }
      if (sendnum[iswap])
        MPI_Send(buf_send_pair,n,MPI_DOUBLE,sendproc[iswap],0,world);
      if (recvnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);

      if (!lmp->kokkos->gpu_aware_flag) {
        k_buf_recv_pair.modify<LMPHostType>();
        k_buf_recv_pair.sync<DeviceType>();
      }
      k_buf_tmp = k_buf_recv_pair;
    } else k_buf_tmp = k_buf_send_pair;

    // unpack buffer

    pairKKBase->unpack_forward_comm_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_tmp);
    DeviceType().fence();
  }
}

/* ---------------------------------------------------------------------- */

void CommKokkos::grow_buf_pair(int n) {
  max_buf_pair = n * BUFFACTOR;
  k_buf_send_pair.resize(max_buf_pair);
  k_buf_recv_pair.resize(max_buf_pair);
}

/* ---------------------------------------------------------------------- */

void CommKokkos::grow_buf_fix(int n) {
  max_buf_fix = n * BUFFACTOR;
  k_buf_send_fix.resize(max_buf_fix);
  k_buf_recv_fix.resize(max_buf_fix);
}

/* ---------------------------------------------------------------------- */

void CommKokkos::reverse_comm(Pair *pair)
{
  if (pair->execution_space == Host || !pair->reverse_comm_device || reverse_pair_comm_classic) {
    k_sendlist.sync<LMPHostType>();
    CommBrick::reverse_comm(pair);
  } else {
    k_sendlist.sync<LMPDeviceType>();
    reverse_comm_device<LMPDeviceType>(pair);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::reverse_comm_device(Pair *pair)
{
  int iswap,n;
  MPI_Request request;
  DAT::tdual_xfloat_1d k_buf_tmp;

  KokkosBase* pairKKBase = dynamic_cast<KokkosBase*>(pair);

  int nsize = MAX(pair->comm_reverse,pair->comm_reverse_off);

  int nmax = max_buf_pair;
  for (iswap = 0; iswap < nswap; iswap++) {
    nmax = MAX(nmax,nsize*sendnum[iswap]);
    nmax = MAX(nmax,nsize*recvnum[iswap]);
  }
  if (nmax > max_buf_pair)
    grow_buf_pair(nmax);

  for (iswap = nswap-1; iswap >= 0; iswap--) {

    // pack buffer

    n = pairKKBase->pack_reverse_comm_kokkos(recvnum[iswap],firstrecv[iswap],k_buf_send_pair);
    DeviceType().fence();

    // exchange with another proc
    // if self, set recv buffer to send buffer

    double* buf_send_pair;
    double* buf_recv_pair;
    if (lmp->kokkos->gpu_aware_flag) {
      buf_send_pair = k_buf_send_pair.view<DeviceType>().data();
      buf_recv_pair = k_buf_recv_pair.view<DeviceType>().data();
    } else {
      k_buf_send_pair.modify<DeviceType>();
      k_buf_send_pair.sync<LMPHostType>();
      buf_send_pair = k_buf_send_pair.h_view.data();
      buf_recv_pair = k_buf_recv_pair.h_view.data();
    }

    if (sendproc[iswap] != me) {
      if (sendnum[iswap])
        MPI_Irecv(buf_recv_pair,nsize*sendnum[iswap],MPI_DOUBLE,sendproc[iswap],0,world,&request);
      if (recvnum[iswap])
        MPI_Send(buf_send_pair,n,MPI_DOUBLE,recvproc[iswap],0,world);
      if (sendnum[iswap]) MPI_Wait(&request,MPI_STATUS_IGNORE);

      if (!lmp->kokkos->gpu_aware_flag) {
        k_buf_recv_pair.modify<LMPHostType>();
        k_buf_recv_pair.sync<DeviceType>();
      }
      k_buf_tmp = k_buf_recv_pair;
    } else k_buf_tmp = k_buf_send_pair;

    // unpack buffer

    auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
    pairKKBase->unpack_reverse_comm_kokkos(sendnum[iswap],k_sendlist_iswap,
                                       k_buf_tmp);
    DeviceType().fence();
  }
}

/* ---------------------------------------------------------------------- */

void CommKokkos::forward_comm(Dump *dump)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::forward_comm(dump);
}

/* ---------------------------------------------------------------------- */

void CommKokkos::reverse_comm(Dump *dump)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::reverse_comm(dump);
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside some proc's box
     can happen if atom moves outside of non-periodic boundary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommKokkos::exchange()
{
  if (!exchange_comm_classic) {
    if (atom->nextra_grow) {

      // check if all fixes with atom-based arrays support exchange on device

      int flag = 1;
      for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
        auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
        if (!fix_iextra->exchange_comm_device) {
          flag = 0;
          break;
        }
      }

      if (!atomKK->avecKK->unpack_exchange_indices_flag || !flag) {
        if (!atomKK->avecKK->unpack_exchange_indices_flag) {
          if (comm->me == 0) {
            error->warning(FLERR,"Atom style not compatible with fix sending data in Kokkos communication, "
                           "switching to classic exchange/border communication");
          }
        } else if (!flag) {
          if (comm->me == 0) {
            error->warning(FLERR,"Fix with atom-based arrays not compatible with sending data in Kokkos communication, "
                           "switching to classic exchange/border communication");
          }
        }
        exchange_comm_classic = true;
        lmp->kokkos->exchange_comm_classic = 1;
      }
    }
  }

  if (!exchange_comm_classic) {
    if (exchange_comm_on_host) exchange_device<LMPHostType>();
    else exchange_device<LMPDeviceType>();
    return;
  }

  atomKK->sync(Host,ALL_MASK);
  CommBrick::exchange();
  atomKK->modified(Host,ALL_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct BuildExchangeListFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  X_FLOAT _lo,_hi;
  typename AT::t_x_array _x;

  int _nlocal,_dim;
  typename AT::t_int_scalar _nsend;
  typename AT::t_int_1d _sendlist;


  BuildExchangeListFunctor(
      const typename AT::tdual_x_array x,
      const typename AT::tdual_int_1d sendlist,
      typename AT::tdual_int_scalar nsend,
      int nlocal, int dim,
      X_FLOAT lo, X_FLOAT hi):
                _lo(lo),_hi(hi),
                _x(x.template view<DeviceType>()),
                _nlocal(nlocal),_dim(dim),
                _nsend(nsend.template view<DeviceType>()),
                _sendlist(sendlist.template view<DeviceType>()) { }

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    if (_x(i,_dim) < _lo || _x(i,_dim) >= _hi) {
      const int mysend = Kokkos::atomic_fetch_add(&_nsend(),1);
      if (mysend < (int)_sendlist.extent(0))
        _sendlist(mysend) = i;
    }
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::exchange_device()
{
  int nsend,nrecv,nrecv1,nrecv2,nlocal;
  double *sublo,*subhi;
  double lo,hi;
  MPI_Request request;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  //   new ghosts are created in borders()
  // map_set() is done at end of borders()

  if (lmp->kokkos->atom_map_classic)
    if (map_style != Atom::MAP_NONE) atom->map_clear();

  // clear ghost count and any ghost bonus data internal to AtomVec

  atom->nghost = 0;
  atom->avec->clear_bonus();

  if (comm->nprocs > 1) { // otherwise no-op

    // subbox bounds for orthogonal or triclinic

    if (triclinic == 0) {
      sublo = domain->sublo;
      subhi = domain->subhi;
    } else {
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    }

    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,ALL_MASK);

    // loop over dimensions
    for (int dim = 0; dim < 3; dim++) {

      lo = sublo[dim];
      hi = subhi[dim];
      nlocal = atom->nlocal;
      nsend = 0;

      // fill buffer with atoms leaving my box, using < and >=

      k_count.h_view() = k_exchange_sendlist.h_view.extent(0);
      while (k_count.h_view() >= (int)k_exchange_sendlist.h_view.extent(0)) {
        k_count.h_view() = 0;
        k_count.modify<LMPHostType>();
        k_count.sync<DeviceType>();

        BuildExchangeListFunctor<DeviceType>
          f(atomKK->k_x,k_exchange_sendlist,k_count,
            nlocal,dim,lo,hi);
        Kokkos::parallel_for(nlocal,f);
        k_exchange_sendlist.modify<DeviceType>();
        k_count.modify<DeviceType>();

        k_count.sync<LMPHostType>();
        int count = k_count.h_view();
        if (count >= (int)k_exchange_sendlist.h_view.extent(0)) {
          MemKK::realloc_kokkos(k_exchange_sendlist,"comm:k_exchange_sendlist",count*1.1);
          MemKK::realloc_kokkos(k_exchange_copylist,"comm:k_exchange_copylist",count*1.1);
          k_count.h_view() = k_exchange_sendlist.h_view.extent(0);
        }
      }
      int count = k_count.h_view();

      // sort exchange_sendlist

      auto d_exchange_sendlist = Kokkos::subview(k_exchange_sendlist.view<DeviceType>(),std::make_pair(0,count));
      Kokkos::sort(DeviceType(), d_exchange_sendlist);
      k_exchange_sendlist.sync<LMPHostType>();

      // when atom is deleted, fill it in with last atom

      int sendpos = count-1;
      int icopy = nlocal-1;
      nlocal -= count;
      for (int recvpos = 0; recvpos < count; recvpos++) {
        int irecv = k_exchange_sendlist.h_view(recvpos);
        if (irecv < nlocal) {
          if (icopy == k_exchange_sendlist.h_view(sendpos)) icopy--;
          while (sendpos > 0 && icopy <= k_exchange_sendlist.h_view(sendpos-1)) {
            sendpos--;
            icopy = k_exchange_sendlist.h_view(sendpos) - 1;
          }
          k_exchange_copylist.h_view(recvpos) = icopy;
          icopy--;
        } else
          k_exchange_copylist.h_view(recvpos) = -1;
      }

      k_exchange_copylist.modify<LMPHostType>();
      k_exchange_copylist.sync<DeviceType>();
      nsend = count;
      if (nsend > maxsend) grow_send_kokkos(nsend,0);
      nsend =
        atomKK->avecKK->pack_exchange_kokkos(count,k_buf_send,
                                   k_exchange_sendlist,k_exchange_copylist,
                                   ExecutionSpaceFromDevice<DeviceType>::space);
      DeviceType().fence();
      atom->nlocal = nlocal;

      // send/recv atoms in both directions
      // send size of message first so receiver can realloc buf_recv if needed
      // if 1 proc in dimension, no send/recv
      //   set nrecv = 0 so buf_send atoms will be lost
      // if 2 procs in dimension, single send/recv
      // if more than 2 procs in dimension, send/recv to both neighbors

      const int data_size = atomKK->avecKK->size_exchange;

      if (procgrid[dim] == 1) nrecv = 0;
      else {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
                     &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,MPI_STATUS_IGNORE);
        nrecv = nrecv1;
        if (procgrid[dim] > 2) {
          MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
                       &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,MPI_STATUS_IGNORE);
          nrecv += nrecv2;
        }
        if (nrecv > maxrecv) grow_recv_kokkos(nrecv);

        MPI_Irecv(k_buf_recv.view<DeviceType>().data(),nrecv1,
                  MPI_DOUBLE,procneigh[dim][1],0,
                  world,&request);
        MPI_Send(k_buf_send.view<DeviceType>().data(),nsend,
                 MPI_DOUBLE,procneigh[dim][0],0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        if (procgrid[dim] > 2) {
          MPI_Irecv(k_buf_recv.view<DeviceType>().data()+nrecv1,
                    nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
                    world,&request);
          MPI_Send(k_buf_send.view<DeviceType>().data(),nsend,
                   MPI_DOUBLE,procneigh[dim][1],0,world);
          MPI_Wait(&request,MPI_STATUS_IGNORE);
        }

        if (nrecv) {

          if (atom->nextra_grow) {
            if ((int) k_indices.extent(0) < nrecv/data_size)
              MemoryKokkos::realloc_kokkos(k_indices,"comm:indices",nrecv/data_size);
          } else if (k_indices.h_view.data())
           k_indices = DAT::tdual_int_1d();


          atom->nlocal = atomKK->avecKK->
            unpack_exchange_kokkos(k_buf_recv,nrecv,atom->nlocal,dim,lo,hi,
                                     ExecutionSpaceFromDevice<DeviceType>::space,k_indices);

          DeviceType().fence();
        }
      }

      if (atom->nextra_grow) {
        for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
          auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
          KokkosBase *kkbase = dynamic_cast<KokkosBase*>(fix_iextra);
          int nextrasend = 0;
          nsend = count;
          if (nsend) {
            if (nsend*fix_iextra->maxexchange > maxsend)
              grow_send_kokkos(nsend*fix_iextra->maxexchange,0);
            nextrasend = kkbase->pack_exchange_kokkos(
              count,k_buf_send,k_exchange_sendlist,k_exchange_copylist,
              ExecutionSpaceFromDevice<DeviceType>::space);
            DeviceType().fence();
          }

          int nextrarecv,nextrarecv1,nextrarecv2;
          if (procgrid[dim] == 1) nextrarecv = 0;
          else {
            MPI_Sendrecv(&nextrasend,1,MPI_INT,procneigh[dim][0],0,
                         &nextrarecv1,1,MPI_INT,procneigh[dim][1],0,
                         world,MPI_STATUS_IGNORE);

            nextrarecv = nextrarecv1;

            if (procgrid[dim] > 2) {
              MPI_Sendrecv(&nextrasend,1,MPI_INT,procneigh[dim][1],0,
                           &nextrarecv2,1,MPI_INT,procneigh[dim][0],0,
                           world,MPI_STATUS_IGNORE);

              nextrarecv += nextrarecv2;
            }

            if (nextrarecv > maxrecv) grow_recv_kokkos(nextrarecv);

            MPI_Irecv(k_buf_recv.view<DeviceType>().data(),nextrarecv1,
                      MPI_DOUBLE,procneigh[dim][1],0,
                      world,&request);
            MPI_Send(k_buf_send.view<DeviceType>().data(),nextrasend,
                     MPI_DOUBLE,procneigh[dim][0],0,world);
            MPI_Wait(&request,MPI_STATUS_IGNORE);

            if (procgrid[dim] > 2) {
              MPI_Irecv(k_buf_recv.view<DeviceType>().data()+nextrarecv1,
                        nextrarecv2,MPI_DOUBLE,procneigh[dim][0],0,
                        world,&request);
              MPI_Send(k_buf_send.view<DeviceType>().data(),nextrasend,
                       MPI_DOUBLE,procneigh[dim][1],0,world);
              MPI_Wait(&request,MPI_STATUS_IGNORE);
            }

            if (nextrarecv) {
              kkbase->unpack_exchange_kokkos(
                k_buf_recv,k_indices,nrecv/data_size,
                nrecv1/data_size,nextrarecv1,
                ExecutionSpaceFromDevice<DeviceType>::space);
              DeviceType().fence();
            }
          }
        }
      }
    }
    atomKK->modified(ExecutionSpaceFromDevice<DeviceType>::space,ALL_MASK);
  }

  if (atom->firstgroupname) {
    /* this is not yet implemented with Kokkos */
    atomKK->sync(Host,ALL_MASK);
    atom->first_reorder();
    atomKK->modified(Host,ALL_MASK);
  }
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate, so don't need to explicitly
     call communicate routine on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommKokkos::borders()
{
  if (!exchange_comm_classic) {

    if (atom->nextra_border || mode != Comm::SINGLE || bordergroup ||
         (ghost_velocity && atomKK->avecKK->no_border_vel_flag)) {

      if (comm->me == 0) {
        error->warning(FLERR,"Required border comm not yet implemented in Kokkos communication, "
                      "switching to classic exchange/border communication");
      }
      exchange_comm_classic = true;
      lmp->kokkos->exchange_comm_classic = 1;
    }
  }

  if (!exchange_comm_classic) {
    if (exchange_comm_on_host) borders_device<LMPHostType>();
    else borders_device<LMPDeviceType>();
  } else {
    atomKK->sync(Host,ALL_MASK);
    k_sendlist.sync<LMPHostType>();
    int prev_auto_sync = lmp->kokkos->auto_sync;
    lmp->kokkos->auto_sync = 1;
    CommBrick::borders();
    lmp->kokkos->auto_sync = prev_auto_sync;
    k_sendlist.modify<LMPHostType>();
    atomKK->modified(Host,ALL_MASK);
  }

  if (comm->nprocs == 1 && !ghost_velocity && !forward_comm_classic)
    copy_swap_info();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct BuildBorderListFunctor {
        typedef DeviceType device_type;
        typedef ArrayTypes<DeviceType> AT;
  X_FLOAT lo,hi;
  typename AT::t_x_array x;
  int iswap,maxsendlist;
  int nfirst,nlast,dim;
  typename AT::t_int_2d sendlist;
  typename AT::t_int_scalar nsend;

  BuildBorderListFunctor(typename AT::tdual_x_array _x,
                         typename AT::tdual_int_2d _sendlist,
                         typename AT::tdual_int_scalar _nsend,int _nfirst,
                         int _nlast, int _dim,
                         X_FLOAT _lo, X_FLOAT _hi, int _iswap,
                         int _maxsendlist):
    lo(_lo),hi(_hi),x(_x.template view<DeviceType>()),iswap(_iswap),
    maxsendlist(_maxsendlist),nfirst(_nfirst),nlast(_nlast),dim(_dim),
    sendlist(_sendlist.template view<DeviceType>()),
    nsend(_nsend.template view<DeviceType>()) {}


  KOKKOS_INLINE_FUNCTION
  void operator() (typename Kokkos::TeamPolicy<DeviceType>::member_type dev) const {
    const int chunk = ((nlast - nfirst + dev.league_size() - 1 ) /
                       dev.league_size());
    const int teamstart = chunk*dev.league_rank() + nfirst;
    const int teamend = (teamstart + chunk) < nlast?(teamstart + chunk):nlast;
    int mysend = 0;
    for (int i=teamstart + dev.team_rank(); i<teamend; i+=dev.team_size()) {
      if (x(i,dim) >= lo && x(i,dim) <= hi) mysend++;
    }
    const int my_store_pos = dev.team_scan(mysend,&nsend());

    if (my_store_pos+mysend < maxsendlist) {
    mysend = my_store_pos;
      for (int i=teamstart + dev.team_rank(); i<teamend; i+=dev.team_size()) {
        if (x(i,dim) >= lo && x(i,dim) <= hi) {
          sendlist(iswap,mysend++) = i;
        }
      }
    }
  }

  size_t shmem_size(const int team_size) const { (void) team_size; return 1000u;}
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommKokkos::borders_device() {
  int i,n,itype,iswap,dim,ineed,twoneed,smax,rmax;
  int nsend,nrecv,sendflag,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *mlo,*mhi;
  MPI_Request request;

  ExecutionSpace exec_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(exec_space,ALL_MASK);
  k_sendlist.sync<DeviceType>();

  int team_size = 1;
  if (exec_space == Device)
    team_size = 128;

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in list for use in future timesteps

      x = atom->x;
      if (mode == Comm::SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;

      // sendflag = 0 if I do not send on this swap
      // sendneed test indicates receiver no longer requires data
      // e.g. due to non-PBC or non-uniform sub-domains

      if (ineed/2 >= sendneed[dim][ineed % 2]) sendflag = 0;
      else sendflag = 1;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus atoms in bordergroup
      // only need to limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost

      if (sendflag) {
        if (!bordergroup || ineed >= 2) {
          if (mode == Comm::SINGLE) {
            k_total_send.h_view() = 0;
            k_total_send.template modify<LMPHostType>();
            k_total_send.template sync<LMPDeviceType>();

            BuildBorderListFunctor<DeviceType> f(atomKK->k_x,k_sendlist,
                k_total_send,nfirst,nlast,dim,lo,hi,iswap,maxsendlist[iswap]);
            Kokkos::TeamPolicy<DeviceType> config((nlast-nfirst+team_size-1)/team_size,team_size);
            Kokkos::parallel_for(config,f);

            k_total_send.template modify<DeviceType>();
            k_total_send.template sync<LMPHostType>();

            k_sendlist.modify<DeviceType>();

            if (k_total_send.h_view() >= maxsendlist[iswap]) {
              grow_list(iswap,k_total_send.h_view());

              k_total_send.h_view() = 0;
              k_total_send.template modify<LMPHostType>();
              k_total_send.template sync<LMPDeviceType>();

              BuildBorderListFunctor<DeviceType> f(atomKK->k_x,k_sendlist,
                  k_total_send,nfirst,nlast,dim,lo,hi,iswap,maxsendlist[iswap]);
              Kokkos::TeamPolicy<DeviceType> config((nlast-nfirst+team_size-1)/team_size,team_size);
              Kokkos::parallel_for(config,f);

              k_total_send.template modify<DeviceType>();
              k_total_send.template sync<LMPHostType>();

              k_sendlist.modify<DeviceType>();
            }
            nsend = k_total_send.h_view();
          } else {
            error->all(FLERR,"Required border comm not yet "
                       "implemented with Kokkos");
            for (i = nfirst; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }

        } else {
          error->all(FLERR,"Required border comm not yet "
                     "implemented with Kokkos");
          if (mode == Comm::SINGLE) {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            for (i = atom->nlocal; i < nlast; i++)
              if (x[i][dim] >= lo && x[i][dim] <= hi) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
          } else {
            ngroup = atom->nfirst;
            for (i = 0; i < ngroup; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
            for (i = atom->nlocal; i < nlast; i++) {
              itype = type[i];
              if (x[i][dim] >= mlo[itype] && x[i][dim] <= mhi[itype]) {
                if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
                sendlist[iswap][nsend++] = i;
              }
            }
          }
        }
      }

      // pack up list of border atoms

      if (nsend*size_border > maxsend)
        grow_send_kokkos(nsend*size_border,0);
      if (ghost_velocity) {
        auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
        n = atomKK->avecKK->
          pack_border_vel_kokkos(nsend,k_sendlist_iswap,k_buf_send,
                                 pbc_flag[iswap],pbc[iswap],exec_space);
        DeviceType().fence();
      } else {
        auto k_sendlist_iswap = Kokkos::subview(k_sendlist,iswap,Kokkos::ALL);
        n = atomKK->avecKK->
          pack_border_kokkos(nsend,k_sendlist_iswap,k_buf_send,
                             pbc_flag[iswap],pbc[iswap],exec_space);
        DeviceType().fence();
      }

      // swap atoms with other proc
      // no MPI calls except SendRecv if nsend/nrecv = 0
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,MPI_STATUS_IGNORE);
        if (nrecv*size_border > maxrecv) grow_recv_kokkos(nrecv*size_border);
        if (nrecv) MPI_Irecv(k_buf_recv.view<DeviceType>().data(),
                             nrecv*size_border,MPI_DOUBLE,
                             recvproc[iswap],0,world,&request);
        if (n) MPI_Send(k_buf_send.view<DeviceType>().data(),n,
                        MPI_DOUBLE,sendproc[iswap],0,world);
        if (nrecv) MPI_Wait(&request,MPI_STATUS_IGNORE);
      } else {
        nrecv = nsend;
      }

      // unpack buffer

      if (ghost_velocity) {
        if (sendproc[iswap] != me) {
          atomKK->avecKK->unpack_border_vel_kokkos(nrecv,atom->nlocal+atom->nghost,
                                         k_buf_recv,exec_space);
          DeviceType().fence();
        } else {
          atomKK->avecKK->unpack_border_vel_kokkos(nrecv,atom->nlocal+atom->nghost,
                                         k_buf_send,exec_space);
          DeviceType().fence();
        }
      } else {
        if (sendproc[iswap] != me) {
          atomKK->avecKK->unpack_border_kokkos(nrecv,atom->nlocal+atom->nghost,
                                     k_buf_recv,exec_space);
          DeviceType().fence();
        } else {
          atomKK->avecKK->unpack_border_kokkos(nrecv,atom->nlocal+atom->nghost,
                                     k_buf_send,exec_space);
          DeviceType().fence();
        }
      }
      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }

  // ensure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send_kokkos(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv_kokkos(max);

  atomKK->modified(exec_space,ALL_MASK);

  // reset global->local map

  if (map_style != Atom::MAP_NONE)
    atom->map_set();
}

/* ----------------------------------------------------------------------
   copy swap info
------------------------------------------------------------------------- */

void CommKokkos::copy_swap_info()
{
  if (nswap > (int)k_swap.extent(1)) {
    k_swap = DAT::tdual_int_2d("comm:swap",2,nswap);
    k_firstrecv    = Kokkos::subview(k_swap,0,Kokkos::ALL);
    k_sendnum_scan = Kokkos::subview(k_swap,1,Kokkos::ALL);
  }
  int scan = 0;
  for (int iswap = 0; iswap < nswap; iswap++) {
    scan += sendnum[iswap];
    k_sendnum_scan.h_view[iswap] = scan;
    k_firstrecv.h_view[iswap] = firstrecv[iswap];
  }
  totalsend = scan;

  // create map of ghost to local atom id
  // store periodic boundary transform from local to ghost

  k_sendlist.sync<LMPHostType>();

  if (totalsend > (int)k_pbc.extent(0)) {
    k_pbc = DAT::tdual_int_2d("comm:pbc",totalsend,6);
    k_swap2 = DAT::tdual_int_2d("comm:swap2",2,totalsend);
    k_pbc_flag = Kokkos::subview(k_swap2,0,Kokkos::ALL);
    k_g2l = Kokkos::subview(k_swap2,1,Kokkos::ALL);
  }

  for (int iswap = 0; iswap < nswap; iswap++) {
    for (int i = 0; i < sendnum[iswap]; i++) {
      int source = sendlist[iswap][i] - atom->nlocal;
      int dest = firstrecv[iswap] + i - atom->nlocal;
      k_pbc_flag.h_view(dest) = pbc_flag[iswap];
      k_pbc.h_view(dest,0) = pbc[iswap][0];
      k_pbc.h_view(dest,1) = pbc[iswap][1];
      k_pbc.h_view(dest,2) = pbc[iswap][2];
      k_pbc.h_view(dest,3) = pbc[iswap][3];
      k_pbc.h_view(dest,4) = pbc[iswap][4];
      k_pbc.h_view(dest,5) = pbc[iswap][5];
      k_g2l.h_view(dest) = atom->nlocal + source;

      if (source >= 0) {
        k_pbc_flag.h_view(dest) = k_pbc_flag.h_view(dest) || k_pbc_flag.h_view(source);
        k_pbc.h_view(dest,0) += k_pbc.h_view(source,0);
        k_pbc.h_view(dest,1) += k_pbc.h_view(source,1);
        k_pbc.h_view(dest,2) += k_pbc.h_view(source,2);
        k_pbc.h_view(dest,3) += k_pbc.h_view(source,3);
        k_pbc.h_view(dest,4) += k_pbc.h_view(source,4);
        k_pbc.h_view(dest,5) += k_pbc.h_view(source,5);
        k_g2l.h_view(dest) = k_g2l.h_view(source);
      }
    }
  }

  k_swap.modify<LMPHostType>();
  k_swap2.modify<LMPHostType>();
  k_pbc.modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommKokkos::grow_send(int n, int flag)
{
  grow_send_kokkos(n,flag,Host);
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommKokkos::grow_recv(int n)
{
  grow_recv_kokkos(n,Host);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommKokkos::grow_send_kokkos(int n, int flag, ExecutionSpace space)
{

  maxsend = static_cast<int> (BUFFACTOR * n);
  int maxsend_border = (maxsend+BUFEXTRA)/atomKK->avecKK->size_border;
  if (flag) {
    if (space == Device)
      k_buf_send.modify<LMPDeviceType>();
    else
      k_buf_send.modify<LMPHostType>();

    if (ghost_velocity)
      k_buf_send.resize(maxsend_border,
                        atomKK->avecKK->size_border + atomKK->avecKK->size_velocity);
    else
      k_buf_send.resize(maxsend_border,atomKK->avecKK->size_border);
    buf_send = k_buf_send.view<LMPHostType>().data();
  } else {
    if (ghost_velocity)
      MemoryKokkos::realloc_kokkos(k_buf_send,"comm:k_buf_send",maxsend_border,
                        atomKK->avecKK->size_border + atomKK->avecKK->size_velocity);
    else
      MemoryKokkos::realloc_kokkos(k_buf_send,"comm:k_buf_send",maxsend_border,
                        atomKK->avecKK->size_border);
    buf_send = k_buf_send.view<LMPHostType>().data();
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommKokkos::grow_recv_kokkos(int n, ExecutionSpace /*space*/)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  int maxrecv_border = (maxrecv+BUFEXTRA)/atomKK->avecKK->size_border;

  MemoryKokkos::realloc_kokkos(k_buf_recv,"comm:k_buf_recv",maxrecv_border,
    atomKK->avecKK->size_border);
  buf_recv = k_buf_recv.view<LMPHostType>().data();
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommKokkos::grow_list(int /*iswap*/, int n)
{
  int size = static_cast<int> (BUFFACTOR * n);

  if (exchange_comm_classic) { // force realloc on Host
    k_sendlist.sync<LMPHostType>();
    k_sendlist.modify<LMPHostType>();
  }

  memoryKK->grow_kokkos(k_sendlist,sendlist,maxswap,size,"comm:sendlist");

  for (int i=0;i<maxswap;i++) {
    maxsendlist[i]=size; sendlist[i]=&k_sendlist.view<LMPHostType>()(i,0);
  }
}

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
------------------------------------------------------------------------- */

void CommKokkos::grow_swap(int n)
{
  free_swap();
  allocate_swap(n);
  if (mode == Comm::MULTI) {
    free_multi();
    allocate_multi(n);
  }

  maxswap = n;
  int size = MAX(k_sendlist.d_view.extent(1),BUFMIN);

  if (exchange_comm_classic) { // force realloc on Host
    k_sendlist.sync<LMPHostType>();
    k_sendlist.modify<LMPHostType>();
  }

  memoryKK->grow_kokkos(k_sendlist,sendlist,maxswap,size,"comm:sendlist");

  memory->grow(maxsendlist,n,"comm:maxsendlist");
  for (int i=0;i<maxswap;i++) maxsendlist[i]=size;
}

/* ----------------------------------------------------------------------
   forward communication of N values in per-atom array
------------------------------------------------------------------------- */

void CommKokkos::forward_comm_array(int nsize, double **array)
{
  k_sendlist.sync<LMPHostType>();
  CommBrick::forward_comm_array(nsize,array);
}
