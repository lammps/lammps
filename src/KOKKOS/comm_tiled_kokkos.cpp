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

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec_kokkos.h"
#include "compute.h"
#include "dump.h"
#include "fix.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "output.h"

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;
static constexpr int BUFEXTRA = 1000;

/* ---------------------------------------------------------------------- */

CommTiledKokkos::CommTiledKokkos(LAMMPS *_lmp) : CommTiled(_lmp)
{
  sendlist = nullptr;
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
}

/* ---------------------------------------------------------------------- */

CommTiledKokkos::~CommTiledKokkos()
{
  memoryKK->destroy_kokkos(k_sendlist,sendlist);
  sendlist = nullptr;
  buf_send = nullptr;
  buf_recv = nullptr;
}

/* ---------------------------------------------------------------------- */

void CommTiledKokkos::init()
{
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

  CommTiled::init();

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

  if (!comm_f_only) { // not all Kokkos atom_vec styles have reverse pack/unpack routines yet
    reverse_comm_classic = true;
    lmp->kokkos->reverse_comm_classic = 1;
  }
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(int dummy)
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

  CommTiled::forward_comm(dummy);
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

    if (comm_x_only) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = atomKK->k_x.view<DeviceType>().data() +
            firstrecv[iswap][i]*atomKK->k_x.view<DeviceType>().extent(1);
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
        atomKK->avecKK->pack_comm_self(sendnum[iswap][nsend],k_sendlist_small,
                        firstrecv[iswap][nrecv],pbc_flag[iswap][nsend],pbc[iswap][nsend]);
        DeviceType().fence();
      }
      if (recvother[iswap]) MPI_Waitall(nrecv,requests,MPI_STATUS_IGNORE);

    } else if (ghost_velocity) {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            forward_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
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
        DeviceType().fence();
        atomKK->avecKK->unpack_comm_vel_kokkos(recvnum[iswap][nrecv],firstrecv[iswap][nrecv],k_buf_send);
        DeviceType().fence();
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(forward_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_comm_vel_kokkos(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                                k_buf_recv_offset);
          DeviceType().fence();
        }
      }

    } else {
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            forward_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
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
        DeviceType().fence();
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          MPI_Waitany(nrecv,requests,&irecv,MPI_STATUS_IGNORE);
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(forward_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_comm_kokkos(recvnum[iswap][irecv],firstrecv[iswap][irecv],
                                   k_buf_recv_offset);
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

void CommTiledKokkos::reverse_comm()
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

  CommTiled::reverse_comm();

  if (comm_f_only)
    atomKK->modified(Host,F_MASK);
  else
    atomKK->modified(Host,ALL_MASK);
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

    if (comm_f_only) {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            reverse_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
          MPI_Irecv(buf,
                    size_reverse_recv[iswap][i],MPI_DOUBLE,sendproc[iswap][i],0,world,&requests[i]);
        }
      }
      if (recvother[iswap]) {
        for (i = 0; i < nrecv; i++) {
          buf = atomKK->k_f.view<DeviceType>().data() +
            firstrecv[iswap][i]*atomKK->k_f.view<DeviceType>().extent(1);
          MPI_Send(buf,size_reverse_send[iswap][i],
                   MPI_DOUBLE,recvproc[iswap][i],0,world);
        }
      }
      if (sendself[iswap]) {
        auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,nsend,Kokkos::ALL);
        atomKK->avecKK->pack_reverse_self(sendnum[iswap][nsend],k_sendlist_small,
                             firstrecv[iswap][nrecv]);
        DeviceType().fence();
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,irecv,Kokkos::ALL);
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(reverse_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap][irecv],k_sendlist_small,
                                      k_buf_recv_offset);
          DeviceType().fence();
        }
      }

    } else {
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          buf = k_buf_recv.view<DeviceType>().data() +
            reverse_recv_offset[iswap][i]*k_buf_recv.view<DeviceType>().extent(1);
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
        DeviceType().fence();
        atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap][nsend],k_sendlist_small,k_buf_send);
        DeviceType().fence();
      }
      if (sendother[iswap]) {
        for (i = 0; i < nsend; i++) {
          MPI_Waitany(nsend,requests,&irecv,MPI_STATUS_IGNORE);
          auto k_sendlist_small = Kokkos::subview(k_sendlist,iswap,irecv,Kokkos::ALL);
          auto k_buf_recv_offset = Kokkos::subview(k_buf_recv,std::pair<int,int>(reverse_recv_offset[iswap][irecv],(int)k_buf_recv.extent(0)),Kokkos::ALL);
          atomKK->avecKK->unpack_reverse_kokkos(sendnum[iswap][irecv],k_sendlist_small,
                               k_buf_recv_offset);
          DeviceType().fence();
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
  atomKK->sync(Host,ALL_MASK);
  CommTiled::exchange();
  atomKK->modified(Host,ALL_MASK);
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
  atomKK->sync(Host,ALL_MASK);
  CommTiled::borders();
  atomKK->modified(Host,ALL_MASK);
  k_sendlist.modify_host();
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Pair
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Pair *pair)
{
  CommTiled::forward_comm(pair);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Pair
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Pair *pair)
{
  CommTiled::reverse_comm(pair);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Bond *bond)
{
  CommTiled::forward_comm(bond);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Bond
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Bond *bond)
{
  CommTiled::reverse_comm(bond);
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
  CommTiled::forward_comm(fix,size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix
   size/nsize used only to set recv buffer limit
   size = 0 (default) -> use comm_forward from Fix
   size > 0 -> Fix passes max size per atom
   the latter is only useful if Fix does several comm modes,
     some are smaller than max stored in its comm_forward
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Fix *fix, int size)
{
  CommTiled::reverse_comm(fix,size);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Fix with variable size data
   query fix for all pack sizes to ensure buf_send is big enough
   handshake sizes before irregular comm to ensure buf_recv is big enough
   NOTE: how to setup one big buf recv with correct offsets ??
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm_variable(Fix *fix)
{
  CommTiled::reverse_comm_variable(fix);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Compute
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Compute *compute)
{
  CommTiled::forward_comm(compute);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Compute
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Compute *compute)
{
  CommTiled::reverse_comm(compute);
}

/* ----------------------------------------------------------------------
   forward communication invoked by a Dump
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm(Dump *dump)
{
  CommTiled::forward_comm(dump);
}

/* ----------------------------------------------------------------------
   reverse communication invoked by a Dump
   nsize used only to set recv buffer limit
------------------------------------------------------------------------- */

void CommTiledKokkos::reverse_comm(Dump *dump)
{
  CommTiled::reverse_comm(dump);
}

/* ----------------------------------------------------------------------
   forward communication of Nsize values in per-atom array
------------------------------------------------------------------------- */

void CommTiledKokkos::forward_comm_array(int nsize, double **array)
{
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

void CommTiledKokkos::grow_recv_kokkos(int n, int flag, ExecutionSpace /*space*/)
{
  if (flag) maxrecv = n;
  else maxrecv = static_cast<int> (BUFFACTOR * n);

  int maxrecv_border = (maxrecv+BUFEXTRA)/atomKK->avecKK->size_border;

  MemoryKokkos::realloc_kokkos(k_buf_recv,"comm:k_buf_recv",maxrecv_border,
    atomKK->avecKK->size_border);
  buf_recv = k_buf_recv.view<LMPHostType>().data();
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommTiledKokkos::grow_list(int iswap, int iwhich, int n)
{
  int size = static_cast<int> (BUFFACTOR * n);

  k_sendlist.sync<LMPHostType>();
  k_sendlist.modify<LMPHostType>();

  if (size > (int)k_sendlist.extent(2)) {
    memoryKK->grow_kokkos(k_sendlist,sendlist,maxswap,maxsend,size,"comm:sendlist");

    for (int i = 0; i < maxswap; i++)
      maxsendlist[iswap][iwhich] = size;
  }
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
  memory->destroy(sendbox_multiold[i]);
  memory->create(sendbox_multiold[i],n,atom->ntypes+1,6,"comm:sendbox_multiold");

  delete [] maxsendlist[i];
  maxsendlist[i] = new int[n];

  for (int j = 0; j < n; j++)
    maxsendlist[i][j] = BUFMIN;

  if (sendlist && !k_sendlist.d_view.data()) {
    for (int ii = 0; ii < maxswap; ii++) {
      if (sendlist[ii]) {
        for (int jj = 0; jj < nprocmax[ii]; jj++)
          memory->destroy(sendlist[ii][jj]);
        delete [] sendlist[ii];
      }
    }
    delete [] sendlist;
  } else {
    memoryKK->destroy_kokkos(k_sendlist,sendlist);
  }

  memoryKK->create_kokkos(k_sendlist,sendlist,maxswap,n,BUFMIN,"comm:sendlist");
}
