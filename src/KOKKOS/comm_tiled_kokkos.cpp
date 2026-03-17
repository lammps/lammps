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

#include "comm_tiled_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "atom_vec_ellipsoid.h"
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
#include "neighbor.h"
#include "output.h"

#include <Kokkos_Sort.hpp>
#include <vector>

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;

/* ---------------------------------------------------------------------- */

CommTiledKokkos::CommTiledKokkos(LAMMPS *_lmp) : CommTiled(_lmp)
{
  sendlist = nullptr;
  maxsendlist = nullptr;
  nprocmaxtot = 0;

  k_exchange_sendlist       = DAT::tdual_int_1d("comm:k_exchange_sendlist",       100);
  k_exchange_copylist       = DAT::tdual_int_1d("comm:k_exchange_copylist",       100);
  k_exchange_sendlist_bonus = DAT::tdual_int_1d("comm:k_exchange_sendlist_bonus", 100);
  k_exchange_copylist_bonus = DAT::tdual_int_1d("comm:k_exchange_copylist_bonus", 100);
  k_count      = DAT::tdual_int_1d("comm:k_count", 2);
  k_total_send = DAT::tdual_int_scalar("comm:k_total_send");
  atomKK = nullptr;
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommTiledKokkos::CommTiledKokkos(LAMMPS *_lmp, Comm *oldcomm) : CommTiled(_lmp,oldcomm)
{
  sendlist = nullptr;
  maxsendlist = nullptr;
  nprocmaxtot = 0;

  k_exchange_sendlist       = DAT::tdual_int_1d("comm:k_exchange_sendlist",       100);
  k_exchange_copylist       = DAT::tdual_int_1d("comm:k_exchange_copylist",       100);
  k_exchange_sendlist_bonus = DAT::tdual_int_1d("comm:k_exchange_sendlist_bonus", 100);
  k_exchange_copylist_bonus = DAT::tdual_int_1d("comm:k_exchange_copylist_bonus", 100);
  k_count      = DAT::tdual_int_1d("comm:k_count", 2);
  k_total_send = DAT::tdual_int_scalar("comm:k_total_send");
  atomKK = nullptr;
}

/* ---------------------------------------------------------------------- */

CommTiledKokkos::~CommTiledKokkos()
{
  memoryKK->destroy_kokkos(k_sendlist,sendlist);
  memory->destroy(maxsendlist);
  sendlist = nullptr;
  maxsendlist = nullptr;
  buf_send = nullptr;
  buf_recv = nullptr;
}

/* ---------------------------------------------------------------------- */

void CommTiledKokkos::init()
{
  atomKK = (AtomKokkos *) atom;
  exchange_comm_legacy = lmp->kokkos->exchange_comm_legacy;
  forward_comm_legacy = lmp->kokkos->forward_comm_legacy;
  forward_pair_comm_legacy = lmp->kokkos->forward_pair_comm_legacy;
  reverse_pair_comm_legacy = lmp->kokkos->reverse_pair_comm_legacy;
  forward_fix_comm_legacy = lmp->kokkos->forward_fix_comm_legacy;
  reverse_comm_legacy = lmp->kokkos->reverse_comm_legacy;
  exchange_comm_on_host = lmp->kokkos->exchange_comm_on_host;
  forward_comm_on_host = lmp->kokkos->forward_comm_on_host;
  reverse_comm_on_host = lmp->kokkos->reverse_comm_on_host;

  CommTiled::init();
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(int dummy)
{
  if (!forward_comm_legacy) {
    if (forward_comm_on_host) forward_comm_device<LMPHostType>();
    else forward_comm_device<LMPDeviceType>();
    return;
  }

  k_sendlist.sync_host();

  if (ghost_velocity)
    atomKK->sync(Host,atomKK->avecKK->datamask_comm_vel);
  else
    atomKK->sync(Host,atomKK->avecKK->datamask_comm);

  CommTiled::forward_comm(dummy);

  if (ghost_velocity)
    atomKK->modified(Host,atomKK->avecKK->datamask_comm_vel);
  else
    atomKK->modified(Host,atomKK->avecKK->datamask_comm);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommTiledKokkos::forward_comm_device()
{
  int i,irecv,n,nsend,nrecv;
  double *buf;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  k_sendlist.sync<DeviceType>();

  for (int iswap = 0; iswap < nswap; iswap++) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (comm_x_only && !atomKK->k_x.NEED_TRANSFORM) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = (double*)atomKK->k_x.view<DeviceType>().data() +
            firstrecv[iswap][i]*atomKK->k_x.view<DeviceType>().extent(1);
          DeviceType().fence();
          MPI_Irecv(buf,size_forward_recv[iswap][i],
                    MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,i,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_kokkos(sendnum[iswap][i],k_sendlist_small,
                              k_buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          DeviceType().fence();
          MPI_Send(k_buf_send.view<DeviceType>().data(),n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
        atomKK->avecKK->pack_comm_self_kokkos(sendnum[iswap][nsend],k_sendlist_small,
                        firstrecv[iswap][nrecv],pbc_flag[iswap][nsend],pbc[iswap][nsend]);
      }
      if (recvother[iswap]) {
        MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);
        DeviceType().fence();
      }

    } else if (ghost_velocity) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            forward_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
          DeviceType().fence();
          MPI_Irecv(buf,
                    size_forward_recv[iswap][i],MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,i,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_vel_kokkos(sendnum[iswap][i],k_sendlist_small,
                                  k_buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          DeviceType().fence();
          MPI_Send(k_buf_send.view<DeviceType>().data(),n,
                   MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
        atomKK->avecKK->pack_comm_vel_kokkos(sendnum[iswap][nsend],k_sendlist_small,
                            k_buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        atomKK->avecKK->unpack_comm_vel_kokkos(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],k_buf_send);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          DeviceType().fence();
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(forward_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_comm_vel_kokkos(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                                k_buf_recv_offset);
        }
      }

    } else {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            forward_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
          DeviceType().fence();
          MPI_Irecv(buf,
                    size_forward_recv[iswap][i],MPI_DOUBLE,recvproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,i,Kokkos::ALL);
          n = atomKK->avecKK->pack_comm_kokkos(sendnum[iswap][i],k_sendlist_small,
                              k_buf_send,pbc_flag[iswap][i],pbc[iswap][i]);
          DeviceType().fence();
          MPI_Send(k_buf_send.view<DeviceType>().data(),n,MPI_DOUBLE,sendproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
        n = atomKK->avecKK->pack_comm_kokkos(sendnum[iswap][nsend],k_sendlist_small,
                        k_buf_send,pbc_flag[iswap][nsend],pbc[iswap][nsend]);
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          DeviceType().fence();
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(forward_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_comm_kokkos(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                                   k_buf_recv_offset);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm()
{
  if (!reverse_comm_legacy) {
    if (reverse_comm_on_host) reverse_comm_device<LMPHostType>();
    else reverse_comm_device<LMPDeviceType>();
    return;
  }

  k_sendlist.sync_host();
  atomKK->sync(Host,atomKK->avecKK->datamask_reverse);

  CommTiled::reverse_comm();

  atomKK->modified(Host,atomKK->avecKK->datamask_reverse);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommTiledKokkos::reverse_comm_device()
{
  int i,irecv,n,nsend,nrecv;
  double *buf;

  // exchange data with another set of procs in each swap
  // post recvs from all procs except self
  // send data to all procs except self
  // copy data to self if sendself is set
  // wait on all procs except self and unpack received data
  // if comm_f_only set, exchange or copy directly from f, don't pack

  k_sendlist.sync<DeviceType>();

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (comm_f_only  && !atomKK->k_f.NEED_TRANSFORM) {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            reverse_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
          DeviceType().fence();
          MPI_Irecv(buf,
                    size_reverse_recv[iswap][i],MPI_DOUBLE,sendproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = (double*)atomKK->k_f.view<DeviceType>().data() +
            firstrecv[iswap][i]*atomKK->k_f.view<DeviceType>().extent(1);
          DeviceType().fence();
          MPI_Send(buf,size_reverse_send[iswap][i],
                   MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
        atomKK->avecKK->pack_reverse_self_kokkos(sendnum[iswap][nsend],k_sendlist_small,
                             firstrecv[iswap][nrecv]);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          DeviceType().fence();
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,irecv,Kokkos::ALL);
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(reverse_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap][irecv],k_sendlist_small,
                                      k_buf_recv_offset);
        }
      }

    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            reverse_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
          DeviceType().fence();
          MPI_Irecv(buf,
                    size_reverse_recv[iswap][i],MPI_DOUBLE,sendproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          n = atomKK->avecKK->pack_reverse_kokkos(recvnum[iswap][i],firstrecv[iswap][i],k_buf_send);
          DeviceType().fence();
          MPI_Send(k_buf_send.view<DeviceType>().data(),n,MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
        atomKK->avecKK->pack_reverse_kokkos(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],k_buf_send);
        atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap][nsend],k_sendlist_small,k_buf_send);
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          DeviceType().fence();
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,irecv,Kokkos::ALL);
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(reverse_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap][irecv],k_sendlist_small,
                               k_buf_recv_offset);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with procs that touch sub-box in each of 3 dims
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside a touching proc's box
     can happen if atom moves outside of non-periodic boundary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommTiledKokkos::exchange()
{
  if (!exchange_comm_legacy) {
    // Check all nextra_grow fixes support device exchange; fall back if not
    if (atom->nextra_grow) {
      for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
        auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
        if (!fix_iextra->exchange_comm_device) {
          if (comm->me == 0)
            error->warning(FLERR,
              "Fix {} not compatible with Kokkos comm, "
              "switching to legacy exchange/border communication",
              fix_iextra->style);
          exchange_comm_legacy = true;
          lmp->kokkos->exchange_comm_legacy = 1;
          break;
        }
      }
    }
  }

  if (!exchange_comm_legacy) {
    if (exchange_comm_on_host) exchange_device<LMPHostType>();
    else                       exchange_device<LMPDeviceType>();
    return;
  }

  atomKK->sync(Host,atomKK->avecKK->datamask_exchange);

  int prev_auto_sync = lmp->kokkos->auto_sync;
  lmp->kokkos->auto_sync = 1;
  CommTiled::exchange();
  lmp->kokkos->auto_sync = prev_auto_sync;

  atomKK->modified(Host,atomKK->avecKK->datamask_exchange);
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created per swap/proc that will be made
   as list is made, actually do communication
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */

void CommTiledKokkos::borders()
{
  if (!exchange_comm_legacy) {
    if (atom->nextra_border || mode != Comm::SINGLE || bordergroup) {
      if (comm->me == 0)
        error->warning(FLERR,
          "Required border comm not yet implemented in Kokkos comm/tiled, "
          "switching to legacy exchange/border communication");
      exchange_comm_legacy = true;
      lmp->kokkos->exchange_comm_legacy = 1;
    }
  }

  if (!exchange_comm_legacy) {
    if (exchange_comm_on_host) borders_device<LMPHostType>();
    else                       borders_device<LMPDeviceType>();
    return;
  }

  k_sendlist.sync_host();
  if (ghost_velocity)
    atomKK->sync(Host,atomKK->avecKK->datamask_border_vel);
  else
    atomKK->sync(Host,atomKK->avecKK->datamask_border);

  int prev_auto_sync = lmp->kokkos->auto_sync;
  lmp->kokkos->auto_sync = 1;
  CommTiled::borders();
  lmp->kokkos->auto_sync = prev_auto_sync;

  k_sendlist.modify_host();
  if (ghost_velocity)
    atomKK->modified(Host,atomKK->avecKK->datamask_border_vel);
  else
    atomKK->modified(Host,atomKK->avecKK->datamask_border);


}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Pair *pair, int size)
{
  k_sendlist.sync_host();
  CommTiled::forward_comm(pair, size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_reverse from Pair
   size > 0 -> Pair passes max size per atom
   the latter is only useful if Pair does several comm modes,
     some are smaller than max stored in its comm_reverse
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Pair *pair, int size)
{
  k_sendlist.sync_host();
  CommTiled::reverse_comm(pair, size);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Bond
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Bond
   size > 0 -> Bond passes max size per atom
   the latter is only useful if Bond does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Bond *bond, int size)
{
  k_sendlist.sync_host();
  CommTiled::forward_comm(bond, size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Bond
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_reverse from Bond
   size > 0 -> Bond passes max size per atom
   the latter is only useful if Bond does several comm modes,
     some are smaller than max stored in its comm_reverse
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Bond *bond, int size)
{
  k_sendlist.sync_host();
  CommTiled::reverse_comm(bond, size);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Fix *fix, int size)
{
  k_sendlist.sync_host();
  CommTiled::forward_comm(fix, size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_reverse from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_reverse
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Fix *fix, int size)
{
  k_sendlist.sync_host();
  CommTiled::reverse_comm(fix, size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for all pack sizes to ensure buf_send is big enough
   handshake sizes before irregular comm to ensure buf_recv is big enough
   NOTE: how to setup one big buf recv with correct offsets ??
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_variable(Fix *fix)
{
  k_sendlist.sync_host();
  CommTiled::reverse_comm_variable(fix);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Compute
   size > 0 -> Compute passes max size per atom
   the latter is only useful if Compute does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Compute *compute, int size)
{
  k_sendlist.sync_host();
  CommTiled::forward_comm(compute, size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_reverse from Compute
   size > 0 -> Compute passes max size per atom
   the latter is only useful if Compute does several comm modes,
     some are smaller than max stored in its comm_reverse
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Compute *compute, int size)
{
  k_sendlist.sync_host();
  CommTiled::reverse_comm(compute, size);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Dump
   size > 0 -> Dump passes max size per atom
   the latter is only useful if Dump does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Dump *dump, int size)
{
  k_sendlist.sync_host();
  CommTiled::forward_comm(dump, size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_reverse from Dump
   size > 0 -> Dump passes max size per atom
   the latter is only useful if Dump does several comm modes,
     some are smaller than max stored in its comm_reverse
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Dump *dump, int size)
{
  k_sendlist.sync_host();
  CommTiled::reverse_comm(dump, size);
}

/* ----------------------------------------------------------------------
   forward communication of Nsize values in per-atom array
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_array(int nsize, double **array)
{
  k_sendlist.sync_host();
  CommTiled::forward_comm_array(nsize,array);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR and bufextra
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_send(int n, int flag)
{
  grow_send_kokkos(n,flag,Host);
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_recv(int n, int flag)
{
  grow_recv_kokkos(n,flag,Host);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_send_kokkos(int n, int flag, ExecutionSpace space)
{

  maxsend = static_cast<int> (BUFFACTOR * n);
  int maxsend_border = (maxsend+Comm::BUFEXTRA)/atomKK->avecKK->size_border;
  if (flag) {
    if (space == Device)
      k_buf_send.modify_device();
    else
      k_buf_send.modify_host();

    if (ghost_velocity)
      k_buf_send.resize(maxsend_border,
                        atomKK->avecKK->size_border + atomKK->avecKK->size_velocity);
    else
      k_buf_send.resize(maxsend_border,atomKK->avecKK->size_border);
    buf_send = k_buf_send.view_host().data();
  } else {
    if (ghost_velocity)
      MemoryKokkos::realloc_kokkos(k_buf_send,"comm:k_buf_send",maxsend_border,
                        atomKK->avecKK->size_border + atomKK->avecKK->size_velocity);
    else
      MemoryKokkos::realloc_kokkos(k_buf_send,"comm:k_buf_send",maxsend_border,
                        atomKK->avecKK->size_border);
    buf_send = k_buf_send.view_host().data();
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_recv_kokkos(int n, int flag, ExecutionSpace /*space*/)
{
  if (flag) maxrecv = n;
  else maxrecv = static_cast<int> (BUFFACTOR * n);

  int maxrecv_border = (maxrecv+Comm::BUFEXTRA)/atomKK->avecKK->size_border;

  MemoryKokkos::realloc_kokkos(k_buf_recv,"comm:k_buf_recv",maxrecv_border,
    atomKK->avecKK->size_border);
  buf_recv = k_buf_recv.view_host().data();
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_list(int iswap, int iwhich, int n)
{
  int size = static_cast<int> (BUFFACTOR * n);

  k_sendlist.sync_host();
  k_sendlist.modify_host();

  memoryKK->grow_kokkos(k_sendlist,sendlist,maxswap,nprocmaxtot,size,"comm:sendlist");

  for (int i = 0; i < maxswap; i++)
    for (int j = 0; j < nprocmaxtot; j++)
      maxsendlist[i][j] = size;
}

/* ----------------------------------------------------------------------
   grow info for swap I, to allow for N procs to communicate with
   ditto for complementary recv for swap I+1 or I-1, as invoked by caller
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_swap_send(int i, int n, int /*nold*/)
{
  delete [] sendproc[i];
  sendproc[i] = new int[n];
  delete [] sendnum[i];
  sendnum[i] = new int[n];

  delete [] size_reverse_recv[i];
  size_reverse_recv[i] = new int[n];
  delete [] reverse_recv_offset[i];
  reverse_recv_offset[i] = new int[n];

  delete [] pbc_flag[i];
  pbc_flag[i] = new int[n];
  memory->destroy(pbc[i]);
  memory->create(pbc[i],n,6,"comm:pbc_flag");
  memory->destroy(sendbox[i]);
  memory->create(sendbox[i],n,6,"comm:sendbox");
  grow_swap_send_multi(i,n);

  if (sendlist && !k_sendlist.view_host().data()) {
    delete [] sendlist;
    delete [] maxsendlist;

    sendlist = nullptr;
    maxsendlist = nullptr;
  } else {
    memoryKK->destroy_kokkos(k_sendlist,sendlist);
    memory->destroy(maxsendlist);
  }

  nprocmaxtot = MAX(nprocmaxtot,n);

  memoryKK->create_kokkos(k_sendlist,sendlist,maxswap,nprocmaxtot,BUFMIN,"comm:sendlist");
  memory->create(maxsendlist,maxswap,nprocmaxtot,"comm:maxsendlist");

  for (int i = 0; i < maxswap; i++)
    for (int j = 0; j < nprocmaxtot; j++)
      maxsendlist[i][j] = BUFMIN;
}
/* ======================================================================
   Functor: find atoms leaving [lo, hi) in dimension _dim.
   BONUS_FLAG=1 simultaneously collects departing ellipsoid-bonus indices.
   ====================================================================== */

template<class DeviceType, int BONUS_FLAG>
struct BuildExchangeListFunctorTiled {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  double _lo, _hi;
  typename AT::t_kkfloat_1d_3_lr _x;
  int _nlocal, _dim;
  typename AT::t_int_1d _nsend;
  typename AT::t_int_1d _sendlist;
  typename AT::t_int_1d _sendlist_bonus;
  typename AT::t_int_1d _bonus_flags;

  BuildExchangeListFunctorTiled(
      const DAT::ttransform_kkfloat_1d_3_lr  x,
      const DAT::tdual_int_1d               sendlist,
      DAT::tdual_int_1d                      nsend,
      const typename AT::tdual_int_1d        sendlist_bonus,
      const typename AT::tdual_int_1d        bonus_flags,
      int nlocal, int dim, double lo, double hi)
      : _lo(lo), _hi(hi),
        _x(x.template view<DeviceType>()),
        _nlocal(nlocal), _dim(dim),
        _nsend(nsend.template view<DeviceType>()),
        _sendlist(sendlist.template view<DeviceType>()),
        _sendlist_bonus(sendlist_bonus.template view<DeviceType>()),
        _bonus_flags(bonus_flags.template view<DeviceType>()) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    if (_x(i,_dim) < _lo || _x(i,_dim) >= _hi) {
      const int mysend = Kokkos::atomic_fetch_add(&_nsend(0), 1);
      if (mysend < (int)_sendlist.extent(0))
        _sendlist(mysend) = i;
      if (BONUS_FLAG) {
        if (_bonus_flags(i) >= 0) {
          const int mysend_bonus = Kokkos::atomic_fetch_add(&_nsend(1), 1);
          if (mysend_bonus < (int)_sendlist_bonus.extent(0))
            _sendlist_bonus(mysend_bonus) = _bonus_flags(i);
        }
      }
    }
  }
};

/* ======================================================================
   Functor: find atoms inside a 3-D sendbox for borders().
   Team-parallel scan; team_scan gives a contiguous store offset into
   k_sendlist[iswap][iproc] without per-atom atomics on the write path.
   ====================================================================== */

template<class DeviceType>
struct BuildBorderListFunctorTiled {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  double xlo,ylo,zlo,xhi,yhi,zhi;
  typename AT::t_kkfloat_1d_3_lr x;
  int iswap, iproc;
  int nfirst, nlast;
  int maxsendlist;
  typename AT::t_int_3d_lr sendlist;
  typename AT::t_int_scalar nsend;

  BuildBorderListFunctorTiled(
      const DAT::ttransform_kkfloat_1d_3_lr _x,
      const DAT::tdual_int_3d_lr            _sendlist,
      const DAT::tdual_int_scalar           _nsend,
      int _nfirst, int _nlast,
      int _iswap,  int _iproc, int _maxsendlist,
      double _xlo, double _ylo, double _zlo,
      double _xhi, double _yhi, double _zhi)
      : xlo(_xlo),ylo(_ylo),zlo(_zlo),
        xhi(_xhi),yhi(_yhi),zhi(_zhi),
        x(_x.template view<DeviceType>()),
        iswap(_iswap), iproc(_iproc),
        nfirst(_nfirst), nlast(_nlast),
        maxsendlist(_maxsendlist),
        sendlist(_sendlist.template view<DeviceType>()),
        nsend(_nsend.template view<DeviceType>()) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(typename Kokkos::TeamPolicy<DeviceType>::member_type dev) const {
    const int league = dev.league_size();
    const int chunk  = ((nlast - nfirst + league - 1) / league);
    const int tstart = chunk * dev.league_rank() + nfirst;
    const int tend   = (tstart + chunk) < nlast ? (tstart + chunk) : nlast;

    int mysend = 0;
    for (int i = tstart + dev.team_rank(); i < tend; i += dev.team_size()) {
      if (x(i,0) >= xlo && x(i,0) < xhi &&
          x(i,1) >= ylo && x(i,1) < yhi &&
          x(i,2) >= zlo && x(i,2) < zhi)
        mysend++;
    }

    const int my_store_pos = dev.team_scan(mysend, &nsend());

    if (my_store_pos + mysend < maxsendlist) {
      mysend = my_store_pos;
      for (int i = tstart + dev.team_rank(); i < tend; i += dev.team_size()) {
        if (x(i,0) >= xlo && x(i,0) < xhi &&
            x(i,1) >= ylo && x(i,1) < yhi &&
            x(i,2) >= zlo && x(i,2) < zhi)
          sendlist(iswap,iproc,mysend++) = i;
      }
    }
  }

  size_t shmem_size(const int /*team_size*/) const { return 1000u; }
};

/* ======================================================================
   exchange_device<DeviceType>
   -------------------------------------------------------------------
   Device-native atom migration for the RCB tiled decomposition.

   Unlike CommKokkos (brick), tiled has nexchproc[dim] exchange
   neighbours per dimension.  The packed buffer is broadcast to ALL of
   them; unpack_exchange_kokkos filters by coordinate range [lo, hi).
   Never calls modified(Host,...) or sets auto_sync, so TransformView
   stays single-owner and the concurrent-modification error cannot fire.
   ====================================================================== */

template<class DeviceType>
void CommTiledKokkos::exchange_device()
{
  int nsend, nrecv, nlocal, nlocal_bonus;
  double *sublo_ptr, *subhi_ptr;

  if (lmp->kokkos->atom_map_legacy)
    if (map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  if (triclinic == 0) {
    sublo_ptr = domain->sublo;
    subhi_ptr = domain->subhi;
  } else {
    sublo_ptr = domain->sublo_lamda;
    subhi_ptr = domain->subhi_lamda;
  }

  const int ellipsoid_flag = atom->ellipsoid_flag;
  const int line_flag      = atom->line_flag;
  const int tri_flag       = atom->tri_flag;
  const int body_flag      = atom->body_flag;

  if (line_flag || tri_flag || body_flag)
    error->all(FLERR,
      "Bonus struct for line/tri/body not yet supported by Kokkos comm/tiled");

  const int bonus_flag = ellipsoid_flag ? 1 : 0;
  DAT::tdual_int_1d k_bonus_flags;

  ExecutionSpace exec_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(exec_space, X_MASK);
  if (bonus_flag) {
    atomKK->sync(exec_space, ELLIPSOID_MASK);
    k_bonus_flags = atomKK->k_ellipsoid;
  }

  dimension = domain->dimension;

  for (int dim = 0; dim < dimension; dim++) {
    const double lo = sublo_ptr[dim];
    const double hi = subhi_ptr[dim];
    nlocal = atom->nlocal;
    if (bonus_flag) nlocal_bonus = atomKK->avecKK->get_status_nlocal_bonus();

    // ---- build sendlist (retry on overflow) ----
    k_count.view_host()(0) = (int)k_exchange_sendlist.view_host().extent(0);
    k_count.view_host()(1) = 0;

    while (k_count.view_host()(0) >= (int)k_exchange_sendlist.view_host().extent(0)) {
      auto d_count = k_count.view<DeviceType>();
      Kokkos::deep_copy(d_count, 0);

      if (bonus_flag) {
        const int ext = (int)k_exchange_sendlist.view_host().extent(0);
        if ((int)k_exchange_sendlist_bonus.extent(0) < ext) {
          MemKK::realloc_kokkos(k_exchange_sendlist_bonus,"comm:k_exchange_sendlist_bonus",ext);
          MemKK::realloc_kokkos(k_exchange_copylist_bonus,"comm:k_exchange_copylist_bonus",ext);
        }
        BuildExchangeListFunctorTiled<DeviceType,1>
          f(atomKK->k_x,k_exchange_sendlist,k_count,
            k_exchange_sendlist_bonus,k_bonus_flags,nlocal,dim,lo,hi);
        Kokkos::parallel_for(nlocal,f);
      } else {
        BuildExchangeListFunctorTiled<DeviceType,0>
          f(atomKK->k_x,k_exchange_sendlist,k_count,
            k_exchange_sendlist_bonus,k_bonus_flags,nlocal,dim,lo,hi);
        Kokkos::parallel_for(nlocal,f);
      }

      k_exchange_sendlist.modify<DeviceType>();
      k_count.modify<DeviceType>();
      k_count.sync_host();
      if (bonus_flag) k_exchange_sendlist_bonus.modify<DeviceType>();

      const int cnt  = k_count.view_host()(0);
      const int cntb = k_count.view_host()(1);

      if (cnt >= (int)k_exchange_sendlist.view_host().extent(0)) {
        MemKK::realloc_kokkos(k_exchange_sendlist,"comm:k_exchange_sendlist",(int)(cnt*1.1));
        MemKK::realloc_kokkos(k_exchange_copylist,"comm:k_exchange_copylist",(int)(cnt*1.1));
        k_count.view_host()(0) = (int)k_exchange_sendlist.view_host().extent(0);
      }
      if (bonus_flag && cntb >= (int)k_exchange_sendlist_bonus.view_host().extent(0)) {
        MemKK::realloc_kokkos(k_exchange_sendlist_bonus,"comm:k_exchange_sendlist_bonus",(int)(cntb*1.1));
        MemKK::realloc_kokkos(k_exchange_copylist_bonus,"comm:k_exchange_copylist_bonus",(int)(cntb*1.1));
        k_count.view_host()(1) = (int)k_exchange_sendlist_bonus.view_host().extent(0);
      }
    }

    const int count = k_count.view_host()(0);

    // ---- sort sendlist on device ----
    {
      auto d_sl = Kokkos::subview(k_exchange_sendlist.view<DeviceType>(),
                                  std::make_pair(0,count));
      Kokkos::sort(DeviceType(), d_sl);
    }
    k_exchange_sendlist.sync_host();

    // ---- build copylist on host (backfill compaction) ----
    {
      int sendpos = count - 1;
      int icopy   = nlocal - 1;
      nlocal -= count;
      for (int recvpos = 0; recvpos < count; recvpos++) {
        const int irecv_idx = k_exchange_sendlist.view_host()(recvpos);
        if (irecv_idx < nlocal) {
          if (icopy == k_exchange_sendlist.view_host()(sendpos)) icopy--;
          while (sendpos > 0 &&
                 icopy <= k_exchange_sendlist.view_host()(sendpos-1)) {
            sendpos--;
            icopy = k_exchange_sendlist.view_host()(sendpos) - 1;
          }
          k_exchange_copylist.view_host()(recvpos) = icopy;
          icopy--;
        } else {
          k_exchange_copylist.view_host()(recvpos) = -1;
        }
      }
    }
    k_exchange_copylist.modify_host();
    k_exchange_copylist.sync<DeviceType>();

    // ---- bonus copylist (ellipsoid) ----
    if (bonus_flag) {
      atomKK->sync(Host, BONUS_MASK);
      const int count_bonus = k_count.view_host()(1);
      {
        auto d_sb = Kokkos::subview(k_exchange_sendlist_bonus.view<DeviceType>(),
                                    std::make_pair(0,count_bonus));
        Kokkos::sort(DeviceType(), d_sb);
      }
      k_exchange_sendlist_bonus.sync_host();

      AtomVecEllipsoid *avec_ellipsoid =
        dynamic_cast<AtomVecEllipsoid*>(atom->style_match("ellipsoid"));
      AtomVecEllipsoid::Bonus *ebonus = nullptr;
      if (avec_ellipsoid) ebonus = avec_ellipsoid->bonus;

      HAT::t_int_1d i2recv;
      MemKK::realloc_kokkos(i2recv,"comm:i2recv",atom->nmax);
      for (int rp = 0; rp < count; rp++) {
        const int ii = k_exchange_sendlist.view_host()(rp);
        i2recv[ii] = rp;
      }

      Kokkos::deep_copy(k_exchange_copylist_bonus.view_host(), -1);
      int sendpos = count_bonus - 1;
      int icopy   = nlocal_bonus - 1;
      nlocal_bonus -= count_bonus;
      for (int recvpos = 0; recvpos < count_bonus; recvpos++) {
        const int irecv_idx = k_exchange_sendlist_bonus.view_host()(recvpos);
        if (irecv_idx < nlocal_bonus) {
          if (icopy == k_exchange_sendlist_bonus.view_host()(sendpos)) icopy--;
          while (sendpos > 0 &&
                 icopy <= k_exchange_sendlist_bonus.view_host()(sendpos-1)) {
            sendpos--;
            icopy = k_exchange_sendlist_bonus.view_host()(sendpos) - 1;
          }
          if (ebonus) {
            const int irecv_all = i2recv[ebonus[irecv_idx].ilocal];
            k_exchange_copylist_bonus.view_host()(irecv_all) = icopy;
          }
          icopy--;
        }
      }
      k_exchange_copylist_bonus.modify_host();
      k_exchange_copylist_bonus.sync<DeviceType>();
    }

    // ---- pack (also compacts local array in-place) ----
    if (nsend > maxsend) grow_send_kokkos(nsend, 0);
    nsend = atomKK->avecKK->pack_exchange_kokkos(
              count, k_buf_send,
              k_exchange_sendlist, k_exchange_copylist,
              k_exchange_copylist_bonus, exec_space);
    buf_send = k_buf_send.view_host().data();
    atom->nlocal = nlocal;
    if (bonus_flag) atomKK->avecKK->set_status_nlocal_bonus(nlocal_bonus);

    const int data_size = atomKK->avecKK->size_exchange;
    const int nexch     = nexchproc[dim];
    if (!nexch) continue;

    // ---- all-to-all size handshake with every exchange proc ----
    std::vector<MPI_Request> reqs(nexch);
    for (int m = 0; m < nexch; m++)
      MPI_Irecv(&exchnum[dim][m],1,MPI_INT,exchproc[dim][m],0,world,&reqs[m]);
    for (int m = 0; m < nexch; m++)
      MPI_Send(&nsend,1,MPI_INT,exchproc[dim][m],0,world);
    MPI_Waitall(nexch,reqs.data(),MPI_STATUS_IGNORE);

    nrecv = 0;
    for (int m = 0; m < nexch; m++) nrecv += exchnum[dim][m];
    if (nrecv > maxrecv) grow_recv_kokkos(nrecv,0);

    // ---- post all Irecvs, broadcast send buffer to every exchange proc ----
    int offset = 0;
    for (int m = 0; m < nexch; m++) {
      if (exchnum[dim][m] > 0) {
        DeviceType().fence();
        MPI_Irecv(k_buf_recv.view<DeviceType>().data() + offset,
                  exchnum[dim][m],MPI_DOUBLE,exchproc[dim][m],0,world,&reqs[m]);
      } else {
        reqs[m] = MPI_REQUEST_NULL;
      }
      offset += exchnum[dim][m];
    }
    if (nsend > 0) {
      DeviceType().fence();
      for (int m = 0; m < nexch; m++)
        MPI_Send(k_buf_send.view<DeviceType>().data(),nsend,MPI_DOUBLE,
                 exchproc[dim][m],0,world);
    }
    MPI_Waitall(nexch,reqs.data(),MPI_STATUS_IGNORE);
    DeviceType().fence();

    // ---- unpack: coordinate-range filter accepts only atoms in [lo,hi) ----
    if (nrecv > 0) {
      if (atom->nextra_grow || atomKK->avecKK->size_exchange_bonus) {
        if ((int)k_indices.extent(0) < nrecv/data_size)
          MemoryKokkos::realloc_kokkos(k_indices,"comm:indices",nrecv/data_size);
      } else if (k_indices.view_host().data()) {
        k_indices = DAT::tdual_int_1d();
      }
      atom->nlocal = atomKK->avecKK->
        unpack_exchange_kokkos(k_buf_recv,nrecv,atom->nlocal,
                               dim,lo,hi,exec_space,k_indices);
    }

    // ---- nextra_grow fix loop ----
    if (atom->nextra_grow) {
      for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
        auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
        KokkosBase *kkbase = dynamic_cast<KokkosBase*>(fix_iextra);
        int nextrasend = 0;
        if (count > 0) {
          if (count*fix_iextra->maxexchange > maxsend)
            grow_send_kokkos(count*fix_iextra->maxexchange,0);
          nextrasend = kkbase->pack_exchange_kokkos(
            count,k_buf_send,k_exchange_sendlist,k_exchange_copylist,exec_space);
        }
        std::vector<int> nextrarecv_per(nexch,0);
        for (int m = 0; m < nexch; m++)
          MPI_Irecv(&nextrarecv_per[m],1,MPI_INT,exchproc[dim][m],0,world,&reqs[m]);
        for (int m = 0; m < nexch; m++)
          MPI_Send(&nextrasend,1,MPI_INT,exchproc[dim][m],0,world);
        MPI_Waitall(nexch,reqs.data(),MPI_STATUS_IGNORE);
        int nextrarecv = 0;
        for (int m = 0; m < nexch; m++) nextrarecv += nextrarecv_per[m];
        if (nextrarecv > maxrecv) grow_recv_kokkos(nextrarecv,0);
        offset = 0;
        for (int m = 0; m < nexch; m++) {
          if (nextrarecv_per[m] > 0) {
            DeviceType().fence();
            MPI_Irecv(k_buf_recv.view<DeviceType>().data()+offset,
                      nextrarecv_per[m],MPI_DOUBLE,exchproc[dim][m],0,world,&reqs[m]);
          } else {
            reqs[m] = MPI_REQUEST_NULL;
          }
          offset += nextrarecv_per[m];
        }
        if (nextrasend > 0) {
          DeviceType().fence();
          for (int m = 0; m < nexch; m++)
            MPI_Send(k_buf_send.view<DeviceType>().data(),nextrasend,MPI_DOUBLE,
                     exchproc[dim][m],0,world);
        }
        MPI_Waitall(nexch,reqs.data(),MPI_STATUS_IGNORE);
        DeviceType().fence();
        if (nextrarecv > 0) {
          const int nrecv_atoms = nrecv/data_size;
          kkbase->unpack_exchange_kokkos(k_buf_recv,k_indices,
            nrecv_atoms,nrecv_atoms,nextrarecv,exec_space);
        }
      }
    }
  } // end dim loop

  if (atom->firstgroupname) {
    atomKK->sync(Host,ALL_MASK);
    atom->first_reorder();
    atomKK->modified(Host,ALL_MASK);
  }
}

/* ======================================================================
   borders_device<DeviceType>
   -------------------------------------------------------------------
   Device-native ghost-atom list construction for the RCB tiled
   decomposition.  Uses BuildBorderListFunctorTiled (team-parallel 3-D
   box test) instead of the CPU loops in CommTiled::borders().
   Calls modified(exec_space,...) — never modified(Host,...) — so
   TransformView stays single-owner.
   ====================================================================== */

template<class DeviceType>
void CommTiledKokkos::borders_device()
{
  int n, nsend, nrecv, nlast;
  int smaxone = 0, smaxall = 0, rmaxone = 0, rmaxall = 0;

  ExecutionSpace exec_space = ExecutionSpaceFromDevice<DeviceType>::space;

  if (mode == Comm::MULTI) neighbor->build_collection(0);

  if (ghost_velocity)
    atomKK->sync(exec_space, atomKK->avecKK->datamask_border_vel);
  else
    atomKK->sync(exec_space, atomKK->avecKK->datamask_border);

  k_sendlist.sync<DeviceType>();

  const int team_size = (exec_space == Device) ? 128 : 1;

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (iswap % 2 == 0)
      nlast = atom->nlocal + atom->nghost;

    int ncountall = 0;

    for (int m = 0; m < nsendproc[iswap]; m++) {
      const double *bbox = sendbox[iswap][m];
      const double xlo = bbox[0], ylo = bbox[1], zlo = bbox[2];
      const double xhi = bbox[3], yhi = bbox[4], zhi = bbox[5];

      bool overflow = true;
      while (overflow) {
        k_total_send.view<DeviceType>()() = 0;
        DeviceType().fence();

        BuildBorderListFunctorTiled<DeviceType> f(
          atomKK->k_x, k_sendlist, k_total_send,
          0, nlast, iswap, m, maxsendlist[iswap][m],
          xlo,ylo,zlo,xhi,yhi,zhi);

        const int league_size = (nlast + team_size - 1) / team_size;
        Kokkos::parallel_for(
          Kokkos::TeamPolicy<DeviceType>(league_size,team_size)
            .set_scratch_size(0,Kokkos::PerTeam(1000)), f);
        DeviceType().fence();
        k_total_send.modify<DeviceType>();
        k_total_send.sync_host();

        const int nfound = k_total_send.view_host()();
        if (nfound > maxsendlist[iswap][m] - 1) {
          k_sendlist.sync_host();
          grow_list(iswap, m, nfound);
          k_sendlist.sync<DeviceType>();
        } else {
          sendnum[iswap][m] = nfound;
          smaxone = MAX(smaxone, nfound);
          ncountall += nfound;
          overflow = false;
        }
      }
    }
    smaxall = MAX(smaxall, ncountall);

    nsend = nsendproc[iswap] - sendself[iswap];
    nrecv = nrecvproc[iswap] - sendself[iswap];

    if (recvother[iswap])
      for (int m = 0; m < nrecv; m++)
        MPI_Irecv(&recvnum[iswap][m],1,MPI_INT,recvproc[iswap][m],0,world,&requests[m]);
    if (sendother[iswap])
      for (int m = 0; m < nsend; m++)
        MPI_Send(&sendnum[iswap][m],1,MPI_INT,sendproc[iswap][m],0,world);
    if (sendself[iswap]) recvnum[iswap][nrecv] = sendnum[iswap][nsend];
    if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

    for (int m = 0; m < nsendproc[iswap]; m++) {
      size_reverse_recv[iswap][m] = sendnum[iswap][m]*size_reverse;
      reverse_recv_offset[iswap][m] =
        (m == 0) ? 0 : reverse_recv_offset[iswap][m-1] + sendnum[iswap][m-1];
    }

    int ncountall_recv = 0;
    for (int m = 0; m < nrecvproc[iswap]; m++) {
      const int nc = recvnum[iswap][m];
      rmaxone = MAX(rmaxone, nc);
      ncountall_recv += nc;
      size_forward_recv[iswap][m] = nc*size_forward;
      size_reverse_send[iswap][m] = nc*size_reverse;
      if (m == 0) {
        firstrecv[iswap][0]           = atom->nlocal + atom->nghost;
        forward_recv_offset[iswap][0] = 0;
      } else {
        firstrecv[iswap][m]           = firstrecv[iswap][m-1] + recvnum[iswap][m-1];
        forward_recv_offset[iswap][m] = forward_recv_offset[iswap][m-1] + recvnum[iswap][m-1];
      }
    }
    rmaxall = MAX(rmaxall, ncountall_recv);

    if (smaxone*size_border > maxsend) grow_send_kokkos(smaxone*size_border,0);
    if (rmaxall*size_border > maxrecv) grow_recv_kokkos(rmaxall*size_border,0);

    if (recvother[iswap]) {
      for (int m = 0; m < nrecv; m++) {
        double *buf = k_buf_recv.view<DeviceType>().data() +
          forward_recv_offset[iswap][m]*k_buf_recv.view<DeviceType>().extent(1);
        DeviceType().fence();
        MPI_Irecv(buf,recvnum[iswap][m]*size_border,
                  MPI_DOUBLE,recvproc[iswap][m],0,world,&requests[m]);
      }
    }
    if (sendother[iswap]) {
      for (int m = 0; m < nsend; m++) {
        auto k_sl = Kokkos::subview(k_sendlist,iswap,m,Kokkos::ALL);
        if (ghost_velocity)
          n = atomKK->avecKK->pack_border_vel_kokkos(sendnum[iswap][m],k_sl,k_buf_send,
                                pbc_flag[iswap][m],pbc[iswap][m],exec_space);
        else
          n = atomKK->avecKK->pack_border_kokkos(sendnum[iswap][m],k_sl,k_buf_send,
                                pbc_flag[iswap][m],pbc[iswap][m],exec_space);
        if (n) { DeviceType().fence();
          MPI_Send(k_buf_send.view<DeviceType>().data(),n,MPI_DOUBLE,sendproc[iswap][m],0,world);
        }
      }
    }
    if (sendself[iswap]) {
      auto k_sl = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
      if (ghost_velocity) {
        atomKK->avecKK->pack_border_vel_kokkos(sendnum[iswap][nsend],k_sl,k_buf_send,
                            pbc_flag[iswap][nsend],pbc[iswap][nsend],exec_space);
        atomKK->avecKK->unpack_border_vel_kokkos(recvnum[iswap][nrecv],
                            firstrecv[iswap][nrecv],k_buf_send,exec_space);
      } else {
        atomKK->avecKK->pack_border_kokkos(sendnum[iswap][nsend],k_sl,k_buf_send,
                            pbc_flag[iswap][nsend],pbc[iswap][nsend],exec_space);
        atomKK->avecKK->unpack_border_kokkos(recvnum[iswap][nrecv],
                            firstrecv[iswap][nrecv],k_buf_send,exec_space);
      }
    }
    if (recvother[iswap]) {
      for (int i = 0; i < nrecv; i++) {
        int irecv;
        MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
        DeviceType().fence();
        auto k_buf_offset = Kokkos::subview(k_buf_recv,
          std::pair<int,int>(forward_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),
          Kokkos::ALL);
        if (ghost_velocity)
          atomKK->avecKK->unpack_border_vel_kokkos(recvnum[iswap][irecv],
                               firstrecv[iswap][irecv],k_buf_offset,exec_space);
        else
          atomKK->avecKK->unpack_border_kokkos(recvnum[iswap][irecv],
                               firstrecv[iswap][irecv],k_buf_offset,exec_space);
      }
    }

    if (nrecvproc[iswap]) {
      const int nlast_recv = nrecvproc[iswap] - 1;
      atom->nghost += forward_recv_offset[iswap][nlast_recv]
                    + recvnum[iswap][nlast_recv];
    }
  } // end iswap loop

  int max = MAX(maxforward*smaxone, maxreverse*rmaxone);
  if (max > maxsend) grow_send_kokkos(max,0);
  max = MAX(maxforward*rmaxall, maxreverse*smaxall);
  if (max > maxrecv) grow_recv_kokkos(max,0);

  k_sendlist.modify<DeviceType>();

  // Mark modified on device — NOT Host — TransformView stays single-owner
  if (ghost_velocity)
    atomKK->modified(exec_space, atomKK->avecKK->datamask_border_vel);
  else
    atomKK->modified(exec_space, atomKK->avecKK->datamask_border);

  if (map_style != Atom::MAP_NONE) atom->map_set();
}
