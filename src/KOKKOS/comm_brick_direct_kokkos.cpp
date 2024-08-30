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

#include "comm_brick_direct_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec_kokkos.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor.h"

// NOTES:
// still need cutoff calculation for nonuniform layout
// need forward_comm_array to test molecular systems
// test msg tags with individual procs as multiple neighbors via big stencil
// test when cutoffs >> box length
// test with triclinic
// doc msg tag logic in code
// doc stencil data structs and logic in code
// CommBrick could use local maxsend in its borders() check for sendlist realloc
//   instead of indexing the swap for each atom

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;

/* ---------------------------------------------------------------------- */

CommBrickDirectKokkos::CommBrickDirectKokkos(LAMMPS *lmp) : CommBrickDirect(lmp)
{
}

/* ---------------------------------------------------------------------- */

CommBrickDirectKokkos::~CommBrickDirectKokkos()
{
  buf_send_direct = nullptr;
  buf_recv_direct = nullptr;
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommBrickDirectKokkos::CommBrickDirectKokkos(LAMMPS *lmp, Comm *oldcomm) : CommBrickDirect(lmp, oldcomm)
{
}

/* ----------------------------------------------------------------------
   create stencil of direct swaps this procs make with each proc in stencil
   direct swap = send and recv
     same proc can appear multiple times in stencil, self proc can also appear
   stencil is used for border and forward and reverse comm
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::setup()
{
  CommBrickDirect::setup();

  MemKK::realloc_kokkos(k_swap2list,"comm_direct:swap2list",ndirect);
  MemKK::realloc_kokkos(k_pbc_flag_direct,"comm_direct:pbc_flag",ndirect);
  MemKK::realloc_kokkos(k_pbc_direct,"comm_direct:pbc",ndirect,6);
  MemKK::realloc_kokkos(k_self_flags,"comm_direct:pbc",ndirect);

  for (int iswap = 0; iswap < ndirect; iswap++) {
    k_swap2list.h_view[iswap] = swap2list[iswap];
    k_pbc_flag_direct.h_view[iswap] = pbc_flag_direct[iswap];
    k_pbc_direct.h_view(iswap,0) = pbc_direct[iswap][0];
    k_pbc_direct.h_view(iswap,1) = pbc_direct[iswap][1];
    k_pbc_direct.h_view(iswap,2) = pbc_direct[iswap][2];
    k_pbc_direct.h_view(iswap,3) = pbc_direct[iswap][3];
    k_pbc_direct.h_view(iswap,4) = pbc_direct[iswap][4];
    k_pbc_direct.h_view(iswap,5) = pbc_direct[iswap][5];
    k_self_flags.h_view(iswap) = proc_direct[iswap]==me;
  }

  k_swap2list.modify_host();
  k_pbc_flag_direct.modify_host();
  k_pbc_direct.modify_host();
  k_self_flags.modify_host();
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
   exchange owned atoms directly with all neighbor procs,
     not via CommBrick 6-way stencil
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::forward_comm(int dummy)
{
  int forward_comm_classic = 0;
  int forward_comm_on_host = 0;

  if (!forward_comm_classic) {
    if (forward_comm_on_host) forward_comm_device<LMPHostType>();
    else forward_comm_device<LMPDeviceType>();
    return;
  }

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

  CommBrickDirect::forward_comm(dummy);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommBrickDirectKokkos::forward_comm_device()
{
  double *buf;

  // post all receives for ghost atoms
  // except for self copies

  int offset;

  int npost = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (size_forward_recv_direct[iswap]) {
      if (comm_x_only) {
        atomKK->k_x.sync<DeviceType>();
        buf = atomKK->k_x.view<DeviceType>().data() + firstrecv_direct[iswap];
      } else {
        offset = recv_offset_forward_direct[iswap];
        buf = k_buf_recv_direct.view<DeviceType>().data() + offset;
      }
      MPI_Irecv(buf,size_forward_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // pack all atom data at once, including copying self data

  k_sendatoms_list.sync<DeviceType>();
  k_swap2list.sync<DeviceType>();
  k_pbc_flag_direct.sync<DeviceType>();
  k_pbc_direct.sync<DeviceType>();
  k_self_flags.sync<DeviceType>();

  if (ghost_velocity) {
    //atomKK->avecKK->pack_comm_vel_direct(totalsend,k_sendatoms_list,
    //                    k_firstrecv,k_pbc_flag_direct,k_pbc_direct,
    //                    k_swap2list,k_buf_send_direct);
  } else {
    atomKK->avecKK->pack_comm_direct(totalsend,k_sendatoms_list,
                        k_sendnum_scan_direct,k_firstrecv_direct,
                        k_pbc_flag_direct,k_pbc_direct,
                        k_swap2list,k_buf_send_direct,k_self_flags);
  }

  // send all owned atoms to receiving procs
  // except for self copies

  offset = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap]) {
      int n = sendnum_direct[iswap]*atomKK->avecKK->size_forward;
      MPI_Send(k_buf_send_direct.view<DeviceType>().data() + offset,n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world);
      offset += n; // TODO: check
    }
  }

  // wait on incoming messages with ghost atoms
  // unpack all messages at once

  if (npost == 0) return;

  MPI_Waitall(npost,requests,MPI_STATUS_IGNORE);

  if (comm_x_only) return;

  if (ghost_velocity) {
    //atomKK->avecKK->unpack_comm_vel_direct(recvnum_direct,firstrecv_direct,buf_recv_direct);
  } else {
    //atomKK->avecKK->unpack_comm_direct(recvnum_direct,firstrecv_direct,buf_recv_direct);
  }
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a forward_comm(), so don't need to explicitly
     call forward_comm() on reneighboring timestep
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
  // loop over conventional 6-way BRICK swaps in 3 dimensions
  // construct BRICK_DIRECT swaps from them
  // unlike borders() in CommBrick, cannot perform borders comm until end
  // this is b/c the swaps take place simultaneously in all dimensions
  //   and thus cannot contain ghost atoms in the forward comm
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::borders()
{
  atomKK->sync(Host,ALL_MASK);
  int prev_auto_sync = lmp->kokkos->auto_sync;
  lmp->kokkos->auto_sync = 1;
  CommBrickDirect::borders();
  lmp->kokkos->auto_sync = prev_auto_sync;
  atomKK->modified(Host,ALL_MASK);

  int maxsend = 0;
  for (int ilist = 0; ilist < maxlist; ilist++)
    maxsend = MAX(maxsend,maxsendatoms_list[ilist]);

  if (k_sendatoms_list.d_view.extent(1) < maxsend)
    MemKK::realloc_kokkos(k_sendatoms_list,"comm_direct:sendatoms_list",maxlist,maxsend);

  if(k_sendnum_scan_direct.extent(0) < nswap) {
    MemKK::realloc_kokkos(k_sendnum_scan_direct,"comm_direct:sendnum_scan",nswap);
    MemKK::realloc_kokkos(k_firstrecv_direct,"comm_direct:firstrecv",nswap);
  }

  for (int ilist = 0; ilist < maxlist; ilist++) {
    const int nsend = sendnum_list[ilist];
    for (int i = 0; i < nsend; i++)
      k_sendatoms_list.h_view(ilist,i) = sendatoms_list[ilist][i];
  }

  int scan = 0;
  for (int iswap = 0; iswap < nswap; iswap++) {
    scan += sendnum_direct[iswap];
    k_sendnum_scan_direct.h_view[iswap] = scan;
    k_firstrecv_direct.h_view[iswap] = firstrecv_direct[iswap];
  }
  totalsend = scan;

  // grow send and recv buffers

  if (totalsend > k_buf_send_direct.d_view.extent(0))
    grow_send_direct(totalsend,0);
}

/* ----------------------------------------------------------------------
   realloc the size of the send_direct buffer as needed with BUFFACTOR
   do not use bufextra as in CommBrick, b/c not using buf_send_direct for exchange()
   flag = 0, don't need to realloc with copy, just free/malloc w/ BUFFACTOR
   flag = 1, realloc with BUFFACTOR
   flag = 2, free/malloc w/out BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::grow_send_direct(int n, int flag)
{
  if (flag == 0) {
    maxsend_direct = static_cast<int> (BUFFACTOR * n);
    MemKK::realloc_kokkos(k_buf_send_direct,"comm:buf_send_direct",maxsend_direct);
  } else if (flag == 1) {
    maxsend_direct = static_cast<int> (BUFFACTOR * n);
    k_buf_send_direct.resize(maxsend_direct);
  } else {
    MemKK::realloc_kokkos(k_buf_send_direct,"comm:buf_send_direct",maxsend_direct);
  }

  buf_send_direct = k_buf_send_direct.h_view.data();
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv_direct buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::grow_recv_direct(int n)
{
  maxrecv_direct = static_cast<int> (BUFFACTOR * n);
  MemKK::realloc_kokkos(k_buf_recv_direct,"comm:buf_recv_direct",maxrecv_direct);
  buf_recv_direct = k_buf_recv_direct.h_view.data();
}

