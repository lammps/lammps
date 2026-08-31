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

using namespace LAMMPS_NS;

static constexpr double BUFFACTOR = 1.5;
static constexpr int BUFMIN = 1024;


/* ----------------------------------------------------------------------
   build one send list on the device
   counts per thread, then a team scan gives each thread a stable slot, so
   the list comes out in increasing atom order exactly as the host build does
------------------------------------------------------------------------- */

template<class DeviceType>
struct CommBrickDirectKokkos_BuildList {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr _x;
  typename AT::t_int_2d_lr _list;
  typename AT::t_int_scalar _nsend;
  int _ilist,_nlocal,_maxlist;
  int _check[3];
  double _lo[3],_hi[3];

  CommBrickDirectKokkos_BuildList(
      const DAT::ttransform_kkfloat_1d_3_lr &x,
      const DAT::tdual_int_2d_lr &list,
      const DAT::tdual_int_scalar &nsend,
      int ilist, int nlocal, int maxlist,
      const int *check, const double *lo, const double *hi):
    _x(x.template view<DeviceType>()),
    _list(list.template view<DeviceType>()),
    _nsend(nsend.template view<DeviceType>()),
    _ilist(ilist),_nlocal(nlocal),_maxlist(maxlist)
  {
    for (int d = 0; d < 3; d++) { _check[d] = check[d]; _lo[d] = lo[d]; _hi[d] = hi[d]; }
  }

  KOKKOS_INLINE_FUNCTION
  bool keep(const int i) const {
    for (int d = 0; d < 3; d++) {
      if (!_check[d]) continue;
      const double v = static_cast<double>(_x(i,d));
      if ((v < _lo[d]) || (v > _hi[d])) return false;
    }
    return true;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (typename Kokkos::TeamPolicy<DeviceType>::member_type dev) const {
    const int chunk = (_nlocal + dev.league_size() - 1) / dev.league_size();
    const int teamstart = chunk*dev.league_rank();
    const int teamend = (teamstart + chunk) < _nlocal ? (teamstart + chunk) : _nlocal;

    int mysend = 0;
    for (int i = teamstart + dev.team_rank(); i < teamend; i += dev.team_size())
      if (keep(i)) mysend++;

    const int my_store_pos = dev.team_scan(mysend,&_nsend());

    if (my_store_pos + mysend < _maxlist) {
      int m = my_store_pos;
      for (int i = teamstart + dev.team_rank(); i < teamend; i += dev.team_size())
        if (keep(i)) _list(_ilist,m++) = i;
    }
  }
};

/* ---------------------------------------------------------------------- */

CommBrickDirectKokkos::CommBrickDirectKokkos(LAMMPS *lmp) : CommBrickDirect(lmp)
{
  totalsend = 0;
  border_device_flag = 0;
  k_total_send = DAT::tdual_int_scalar("comm_direct:total_send");
}

/* ---------------------------------------------------------------------- */

CommBrickDirectKokkos::~CommBrickDirectKokkos()
{
  // the buffers are owned by the dual views, not by the base class

  buf_send_direct = nullptr;
  buf_recv_direct = nullptr;

  memoryKK->destroy_kokkos(k_sendatoms_list,sendatoms_list);
  sendatoms_list = nullptr;
}

/* ---------------------------------------------------------------------- */
//IMPORTANT: we *MUST* pass "*oldcomm" to the Comm initializer here, as
//           the code below *requires* that the (implicit) copy constructor
//           for Comm is run and thus creating a shallow copy of "oldcomm".
//           The call to Comm::copy_arrays() then converts the shallow copy
//           into a deep copy of the class with the new layout.

CommBrickDirectKokkos::CommBrickDirectKokkos(LAMMPS *lmp, Comm *oldcomm) :
  CommBrickDirect(lmp, oldcomm)
{
  totalsend = 0;
  border_device_flag = 0;
  k_total_send = DAT::tdual_int_scalar("comm_direct:total_send");
}

/* ----------------------------------------------------------------------
   replace the base class's per-list arrays with the host side of a dual view
   so the lists built on the device are the same memory the host routines read
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::lists_to_kokkos()
{
  if (sendatoms_list) {
    for (int ilist = 0; ilist < maxlist; ilist++)
      memory->destroy(sendatoms_list[ilist]);
    memory->sfree(sendatoms_list);
    sendatoms_list = nullptr;
  }

  k_sendatoms_list = DAT::tdual_int_2d_lr();
  memoryKK->create_kokkos(k_sendatoms_list,sendatoms_list,maxlist,BUFMIN,
                          "comm_direct:sendatoms_list");

  for (int ilist = 0; ilist < maxlist; ilist++) {
    maxsendatoms_list[ilist] = BUFMIN;
    sendatoms_list[ilist] = &k_sendatoms_list.view_host()(ilist,0);
  }
}

/* ---------------------------------------------------------------------- */

void CommBrickDirectKokkos::allocate_lists()
{
  CommBrickDirect::allocate_lists();
  lists_to_kokkos();
}

/* ---------------------------------------------------------------------- */

void CommBrickDirectKokkos::deallocate_lists(int nlist)
{
  memoryKK->destroy_kokkos(k_sendatoms_list,sendatoms_list);
  sendatoms_list = nullptr;
  CommBrickDirect::deallocate_lists(nlist);
}

/* ----------------------------------------------------------------------
   the dual view is rectangular, so every list grows together
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::grow_list_direct(int /*ilist*/, int n)
{
  const int size = static_cast<int> (BUFFACTOR * n);

  memoryKK->grow_kokkos(k_sendatoms_list,sendatoms_list,maxlist,size,
                        "comm_direct:sendatoms_list");

  for (int i = 0; i < maxlist; i++) {
    maxsendatoms_list[i] = size;
    sendatoms_list[i] = &k_sendatoms_list.view_host()(i,0);
  }
}

/* ----------------------------------------------------------------------
   create stencil of direct swaps this proc makes with each proc in stencil
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::setup()
{
  CommBrickDirect::setup();

  MemKK::realloc_kokkos(k_swap2list,"comm_direct:swap2list",ndirect);
  MemKK::realloc_kokkos(k_pbc_flag_direct,"comm_direct:pbc_flag",ndirect);
  MemKK::realloc_kokkos(k_pbc_direct,"comm_direct:pbc",ndirect,6);
  MemKK::realloc_kokkos(k_self_flag,"comm_direct:self_flag",ndirect);

  for (int iswap = 0; iswap < ndirect; iswap++) {
    k_swap2list.view_host()[iswap] = swap2list[iswap];
    k_pbc_flag_direct.view_host()[iswap] = pbc_flag_direct[iswap];
    for (int m = 0; m < 6; m++)
      k_pbc_direct.view_host()(iswap,m) = pbc_direct[iswap][m];
    k_self_flag.view_host()(iswap) = proc_direct[iswap] == me;
  }

  k_swap2list.modify_host();
  k_pbc_flag_direct.modify_host();
  k_pbc_direct.modify_host();
  k_self_flag.modify_host();
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::forward_comm(int dummy)
{
  // the device path below packs coords only.  ghost velocities and atom
  // styles with extra forward-comm fields both clear comm_x_only and are
  // not handled there, so they use the host path instead

  if (comm_x_only && !atomKK->k_x.NEED_TRANSFORM) {
    if (lmp->kokkos->forward_comm_on_host) forward_comm_device<LMPHostType>();
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
  // post all receives for ghost atoms, except for swaps with self
  // comm_x_only, so receive straight into the ghost region of x

  int npost = 0;
  double *xdata = (double *) atomKK->k_x.view<DeviceType>().data();
  const int xcols = atomKK->k_x.view<DeviceType>().extent(1);

  DeviceType().fence();

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (size_forward_recv_direct[iswap]) {
      MPI_Irecv(xdata + firstrecv_direct[iswap]*xcols,
                size_forward_recv_direct[iswap],MPI_DOUBLE,
                proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
    }
  }

  // pack every swap in one kernel, including copies to self

  k_sendatoms_list.sync<DeviceType>();
  k_swap2list.sync<DeviceType>();
  k_pbc_flag_direct.sync<DeviceType>();
  k_pbc_direct.sync<DeviceType>();
  k_self_flag.sync<DeviceType>();
  k_sendnum_scan_direct.sync<DeviceType>();
  k_firstrecv_direct.sync<DeviceType>();

  if (totalsend)
    atomKK->avecKK->pack_comm_direct(totalsend,k_sendatoms_list,
                                     k_sendnum_scan_direct,k_firstrecv_direct,
                                     k_pbc_flag_direct,k_pbc_direct,
                                     k_swap2list,k_buf_send_direct,k_self_flag);

  DeviceType().fence();

  // each swap already occupies its own rows of the send buffer,
  // so post every message at once with non-blocking sends

  int nsendpost = 0;
  int offset = 0;
  double *sbuf = k_buf_send_direct.view<DeviceType>().data();

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (sendnum_direct[iswap]) {
      int n = sendnum_direct[iswap]*atomKK->avecKK->size_forward;
      if (proc_direct[iswap] != me)
        MPI_Isend(sbuf + offset,n,MPI_DOUBLE,proc_direct[iswap],
                  sendtag[iswap],world,&send_requests[nsendpost++]);
      offset += n;
    }
  }

  if (npost) MPI_Waitall(npost,requests,MPI_STATUS_IGNORE);
  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);

  // MPI wrote the received coords straight into x on this side

  atomKK->modified(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::reverse_comm()
{
  if (lmp->kokkos->reverse_comm_on_host) reverse_comm_device<LMPHostType>();
  else reverse_comm_device<LMPDeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommBrickDirectKokkos::reverse_comm_device()
{
  // buffer offsets are in atoms; reverse comm moves size_reverse per atom
  // recv_offset_reverse_atoms scans sendnum over swaps, so it indexes the
  //   region each swap's owned-atom contributions arrive in

  const int nrev = size_reverse;

  // post all receives for owned atoms, except for swaps with self

  int npost = 0;

  DeviceType().fence();

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;
    MPI_Irecv(k_buf_recv_reverse.view<DeviceType>().data() +
              nrev*recv_offset_reverse_atoms[iswap],
              nrev*sendnum_direct[iswap],MPI_DOUBLE,
              proc_direct[iswap],sendtag[iswap],world,&requests[npost++]);
  }

  // copy/sum to self on device
  // reads the ghost region and sums into owned atoms, so it is disjoint from
  //   the sends below and can run before them

  k_sendatoms_list.sync<DeviceType>();

  for (int iself = 0; iself < nself_direct; iself++) {
    const int iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    auto k_list = Kokkos::subview(k_sendatoms_list,swap2list[iswap],Kokkos::ALL);
    atomKK->avecKK->pack_reverse_self_kokkos(sendnum_direct[iswap],k_list,
                                             firstrecv_direct[iswap]);
  }

  DeviceType().fence();

  // send ghost contributions to the procs that own those atoms
  // comm_f_only sends straight out of the ghost region of f, which the
  //   unpack below never touches, so the sends can stay in flight

  int nsendpost = 0;

  if (comm_f_only && !atomKK->k_f.NEED_TRANSFORM) {

    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,F_MASK);
    DeviceType().fence();

    double *fdata = (double *) atomKK->k_f.view<DeviceType>().data();
    const int fcols = atomKK->k_f.view<DeviceType>().extent(1);

    for (int iswap = 0; iswap < ndirect; iswap++) {
      if (proc_direct[iswap] == me) continue;
      if (recvnum_direct[iswap] == 0) continue;
      MPI_Isend(fdata + firstrecv_direct[iswap]*fcols,
                nrev*recvnum_direct[iswap],MPI_DOUBLE,proc_direct[iswap],
                recvtag[iswap],world,&send_requests[nsendpost++]);
    }

  } else {

    int send_offset = 0;

    for (int iswap = 0; iswap < ndirect; iswap++) {
      if (proc_direct[iswap] == me) continue;
      if (recvnum_direct[iswap] == 0) continue;
      auto k_sub = Kokkos::subview(k_buf_send_reverse,
                                   Kokkos::make_pair(send_offset,
                                     send_offset+recvnum_direct[iswap]),
                                   Kokkos::ALL);
      const int n = atomKK->avecKK->pack_reverse_kokkos(recvnum_direct[iswap],
                                                        firstrecv_direct[iswap],k_sub);
      if (n) {
        DeviceType().fence();
        MPI_Isend(k_buf_send_reverse.view<DeviceType>().data() + nrev*send_offset,
                  n,MPI_DOUBLE,proc_direct[iswap],recvtag[iswap],world,
                  &send_requests[nsendpost++]);
        send_offset += recvnum_direct[iswap];
      }
    }
  }

  // wait on incoming messages, summing each into the owned atoms as it lands

  for (int i = 0; i < npost; i++) {
    int irecv;
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    const int iswap = send_indices_direct[irecv];
    auto k_list = Kokkos::subview(k_sendatoms_list,swap2list[iswap],Kokkos::ALL);
    auto k_sub = Kokkos::subview(k_buf_recv_reverse,
                                 Kokkos::make_pair(recv_offset_reverse_atoms[iswap],
                                   recv_offset_reverse_atoms[iswap]+sendnum_direct[iswap]),
                                 Kokkos::ALL);
    atomKK->avecKK->unpack_reverse_kokkos(sendnum_direct[iswap],k_list,k_sub);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);

  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   CommBrickDirect does not define its own exchange, it uses CommBrick's
     dimension-by-dimension one, and the device version of that lives in
     CommBrickKokkos::exchange_device().  CommBrickKokkos is a sibling of
     this class rather than a parent, so that code cannot be reached from
     here; running it on the device would mean either duplicating it and the
     dual views it owns, or rearranging the class hierarchy so the direct
     stencil and the KOKKOS exchange can be composed.  Until then the atoms
     have to be synced to the host for the inherited routine.
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::exchange()
{
  atomKK->sync(Host,ALL_MASK);
  CommBrickDirect::exchange();
  atomKK->modified(Host,ALL_MASK);
}


/* ----------------------------------------------------------------------
   build the per-swap lists of owned atoms, one kernel per active list
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::build_lists()
{
  if (!border_device_flag) {
    CommBrickDirect::build_lists();
    return;
  }

  if (lmp->kokkos->exchange_comm_on_host) build_lists_device<LMPHostType>();
  else build_lists_device<LMPDeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommBrickDirectKokkos::build_lists_device()
{
  const ExecutionSpace exec = ExecutionSpaceFromDevice<DeviceType>::space;
  const int nlocal = atom->nlocal;
  const int dim = domain->dimension;

  atomKK->sync(exec,X_MASK);
  k_sendatoms_list.sync<DeviceType>();

  const int team_size = (exec == Device) ? 128 : 1;
  const int nteam = (nlocal + team_size - 1) / team_size;

  for (int ilist = 0; ilist < maxlist; ilist++) {
    if (!active_list[ilist]) continue;

    int check[3];
    double lo[3],hi[3];
    for (int d = 0; d < 3; d++) {
      check[d] = check_list[ilist][d];
      lo[d] = bounds_list[ilist][d][0];
      hi[d] = bounds_list[ilist][d][1];
    }
    if (dim == 2) check[2] = 0;

    int nsend = 0;

    for (int attempt = 0; attempt < 2; attempt++) {
      k_total_send.view_host()() = 0;
      k_total_send.modify_host();
      k_total_send.sync<DeviceType>();

      CommBrickDirectKokkos_BuildList<DeviceType>
        f(atomKK->k_x,k_sendatoms_list,k_total_send,ilist,nlocal,
          maxsendatoms_list[ilist],check,lo,hi);
      Kokkos::TeamPolicy<DeviceType> config(nteam > 0 ? nteam : 1,team_size);
      Kokkos::parallel_for(config,f);

      k_total_send.template modify<DeviceType>();
      k_total_send.sync_host();
      nsend = k_total_send.view_host()();

      // the kernel only stores when the whole team fits, so on overflow the
      //   list has to grow and the pass be repeated

      if (nsend < maxsendatoms_list[ilist]) break;
      grow_list_direct(ilist,nsend);
      k_sendatoms_list.sync<DeviceType>();
    }

    sendnum_list[ilist] = nsend;
  }

  k_sendatoms_list.template modify<DeviceType>();

  // the per-object comm routines still run on the host and read
  //   sendatoms_list directly, so the host side has to be current

  k_sendatoms_list.sync_host();
}

/* ----------------------------------------------------------------------
   exchange border data for every swap on the device
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::borders_comm()
{
  if (!border_device_flag) {
    CommBrickDirect::borders_comm();
    return;
  }

  if (lmp->kokkos->exchange_comm_on_host) borders_comm_device<LMPHostType>();
  else borders_comm_device<LMPDeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void CommBrickDirectKokkos::borders_comm_device()
{
  const ExecutionSpace exec = ExecutionSpaceFromDevice<DeviceType>::space;
  const int ncol = ghost_velocity ?
    size_border + atomKK->avecKK->size_velocity : size_border;

  // buffers are indexed in atoms: sends carry sendnum per swap (ssum in all),
  //   receives carry recvnum per swap (rsum in all)

  if ((ssum_direct + 1 > (int) k_buf_send_border.view_device().extent(0)) ||
      (ncol != (int) k_buf_send_border.view_device().extent(1)))
    MemKK::realloc_kokkos(k_buf_send_border,"comm:buf_send_border",ssum_direct+1,ncol);

  if ((rsum_direct + 1 > (int) k_buf_recv_border.view_device().extent(0)) ||
      (ncol != (int) k_buf_recv_border.view_device().extent(1)))
    MemKK::realloc_kokkos(k_buf_recv_border,"comm:buf_recv_border",rsum_direct+1,ncol);

  k_sendatoms_list.sync<DeviceType>();

  // post all receives for ghost atoms, except for swaps with self
  // recv_offset_forward_atoms scans recvnum, so it indexes each swap's region

  int npost = 0;

  DeviceType().fence();

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (recvnum_direct[iswap] == 0) continue;
    MPI_Irecv(k_buf_recv_border.view<DeviceType>().data() +
              ncol*recv_offset_forward_atoms[iswap],
              ncol*recvnum_direct[iswap],MPI_DOUBLE,
              proc_direct[iswap],recvtag[iswap],world,&requests[npost++]);
  }

  // copies to self go through the send buffer, so do them before the sends
  //   claim it

  for (int iself = 0; iself < nself_direct; iself++) {
    const int iswap = self_indices_direct[iself];
    if (sendnum_direct[iswap] == 0) continue;
    auto k_list = Kokkos::subview(k_sendatoms_list,swap2list[iswap],Kokkos::ALL);
    if (ghost_velocity) {
      atomKK->avecKK->pack_border_vel_kokkos(sendnum_direct[iswap],k_list,
                                             k_buf_send_border,pbc_flag_direct[iswap],
                                             pbc_direct[iswap],exec);
      atomKK->avecKK->unpack_border_vel_kokkos(recvnum_direct[iswap],
                                               firstrecv_direct[iswap],
                                               k_buf_send_border,exec);
    } else {
      atomKK->avecKK->pack_border_kokkos(sendnum_direct[iswap],k_list,
                                         k_buf_send_border,pbc_flag_direct[iswap],
                                         pbc_direct[iswap],exec);
      atomKK->avecKK->unpack_border_kokkos(recvnum_direct[iswap],
                                           firstrecv_direct[iswap],
                                           k_buf_send_border,exec);
    }
  }

  // pack each remaining swap into its own region and send it

  int nsendpost = 0;
  int send_offset = 0;

  for (int iswap = 0; iswap < ndirect; iswap++) {
    if (proc_direct[iswap] == me) continue;
    if (sendnum_direct[iswap] == 0) continue;

    auto k_list = Kokkos::subview(k_sendatoms_list,swap2list[iswap],Kokkos::ALL);
    auto k_sub = Kokkos::subview(k_buf_send_border,
                                 Kokkos::make_pair(send_offset,
                                   send_offset+sendnum_direct[iswap]),
                                 Kokkos::ALL);
    int n;
    if (ghost_velocity)
      n = atomKK->avecKK->pack_border_vel_kokkos(sendnum_direct[iswap],k_list,k_sub,
                                                 pbc_flag_direct[iswap],
                                                 pbc_direct[iswap],exec);
    else
      n = atomKK->avecKK->pack_border_kokkos(sendnum_direct[iswap],k_list,k_sub,
                                             pbc_flag_direct[iswap],
                                             pbc_direct[iswap],exec);
    if (n) {
      DeviceType().fence();
      MPI_Isend(k_buf_send_border.view<DeviceType>().data() + ncol*send_offset,
                n,MPI_DOUBLE,proc_direct[iswap],sendtag[iswap],world,
                &send_requests[nsendpost++]);
      send_offset += sendnum_direct[iswap];
    }
  }

  // unpack each message as it arrives

  for (int ipost = 0; ipost < npost; ipost++) {
    int irecv;
    MPI_Waitany(npost,requests,&irecv,MPI_STATUS_IGNORE);
    const int iswap = recv_indices_direct[irecv];
    auto k_sub = Kokkos::subview(k_buf_recv_border,
                                 Kokkos::make_pair(recv_offset_forward_atoms[iswap],
                                   recv_offset_forward_atoms[iswap]+recvnum_direct[iswap]),
                                 Kokkos::ALL);
    if (ghost_velocity)
      atomKK->avecKK->unpack_border_vel_kokkos(recvnum_direct[iswap],
                                               firstrecv_direct[iswap],k_sub,exec);
    else
      atomKK->avecKK->unpack_border_kokkos(recvnum_direct[iswap],
                                           firstrecv_direct[iswap],k_sub,exec);
  }

  if (nsendpost) MPI_Waitall(nsendpost,send_requests,MPI_STATUS_IGNORE);

  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
------------------------------------------------------------------------- */

void CommBrickDirectKokkos::borders()
{
  // decide once whether the device path can handle this configuration
  // extra border data, multi mode and a border group are not implemented there

  if (!lmp->kokkos->exchange_comm_legacy) {
    if (atom->nextra_border || mode != Comm::SINGLE || bordergroup) {
      if (me == 0)
        error->warning(FLERR,"Required border comm not yet implemented in Kokkos "
                       "communication, switching to legacy exchange/border communication");
      lmp->kokkos->exchange_comm_legacy = 1;
    }
  }

  border_device_flag = !lmp->kokkos->exchange_comm_legacy;

  if (border_device_flag) {

    // build_lists() and borders_comm() run on the device

    CommBrickDirect::borders();

  } else {

    if (ghost_velocity)
      atomKK->sync(Host,atomKK->avecKK->datamask_border_vel);
    else
      atomKK->sync(Host,atomKK->avecKK->datamask_border);
    k_sendatoms_list.sync_host();
    int prev_auto_sync = lmp->kokkos->auto_sync;
    lmp->kokkos->auto_sync = 1;
    CommBrickDirect::borders();
    lmp->kokkos->auto_sync = prev_auto_sync;
    k_sendatoms_list.modify_host();
    if (ghost_velocity)
      atomKK->modified(Host,atomKK->avecKK->datamask_border_vel);
    else
      atomKK->modified(Host,atomKK->avecKK->datamask_border);
  }

  // sendatoms_list is the host side of k_sendatoms_list, so the lists the
  //   build above wrote are already the device view's host data; it is grown
  //   only through grow_list_direct(), which keeps the row pointers in step

  if ((int) k_sendnum_scan_direct.extent(0) < ndirect) {
    MemKK::realloc_kokkos(k_sendnum_scan_direct,"comm_direct:sendnum_scan",ndirect);
    MemKK::realloc_kokkos(k_firstrecv_direct,"comm_direct:firstrecv",ndirect);
  }

  int scan = 0;
  for (int iswap = 0; iswap < ndirect; iswap++) {
    scan += sendnum_direct[iswap];
    k_sendnum_scan_direct.view_host()[iswap] = scan;
    k_firstrecv_direct.view_host()[iswap] = firstrecv_direct[iswap];
  }
  totalsend = scan;

  if (totalsend > (int) k_buf_send_direct.view_device().extent(0))
    grow_send_direct(totalsend*size_forward,0);

  // reverse comm buffers, indexed in atoms with size_reverse per atom
  // sends carry the ghosts this proc received (rsum), receives carry the
  //   owned-atom contributions coming back from every swap (ssum)

  if ((rsum_direct > (int) k_buf_send_reverse.view_device().extent(0)) ||
      (size_reverse != (int) k_buf_send_reverse.view_device().extent(1)))
    MemKK::realloc_kokkos(k_buf_send_reverse,"comm:buf_send_reverse",
                          rsum_direct+1,size_reverse);

  if ((ssum_direct > (int) k_buf_recv_reverse.view_device().extent(0)) ||
      (size_reverse != (int) k_buf_recv_reverse.view_device().extent(1)))
    MemKK::realloc_kokkos(k_buf_recv_reverse,"comm:buf_recv_reverse",
                          ssum_direct+1,size_reverse);

  k_sendatoms_list.modify_host();
  k_sendnum_scan_direct.modify_host();
  k_firstrecv_direct.modify_host();
}

/* ---------------------------------------------------------------------- */

void CommBrickDirectKokkos::grow_send_direct(int n, int flag)
{
  const int nrow = static_cast<int> (BUFFACTOR * n) / 3 + 1;

  if (flag == 1) {
    k_buf_send_direct.resize(nrow,3);
    k_buf_send_direct.clear_sync_state();
  } else {
    MemKK::realloc_kokkos(k_buf_send_direct,"comm:buf_send_direct",nrow,3);
  }
  maxsend_direct = nrow*3;
  buf_send_direct = k_buf_send_direct.view_host().data();
}

/* ---------------------------------------------------------------------- */

void CommBrickDirectKokkos::grow_recv_direct(int n)
{
  const int nrow = static_cast<int> (BUFFACTOR * n) / 3 + 1;

  MemKK::realloc_kokkos(k_buf_recv_direct,"comm:buf_recv_direct",nrow,3);
  maxrecv_direct = nrow*3;
  buf_recv_direct = k_buf_recv_direct.view_host().data();
}
