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

/* ----------------------------------------------------------------------
   Contributing authors:
      Yizhong R. Hu
      Marc L. Descoteaux
      Ulrik Unneberg
      William C. Witt
      Affiliation: Harvard University
------------------------------------------------------------------------- */

#include "pair_dispersion_d3_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDispersionD3Kokkos<DeviceType>::PairDispersionD3Kokkos(LAMMPS *lmp) : PairDispersionD3(lmp)
{
  kokkosable = 1;

  // cn and dc6 are reduced with reverse_comm() and are device resident, so
  // let CommKokkos use the device pack/unpack methods below
  reverse_comm_device = 1;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDispersionD3Kokkos<DeviceType>::~PairDispersionD3Kokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    cutsq = nullptr; // prevent base destructor from double-freeing
  }
}

/* ----------------------------------------------------------------------
   Calculate coordination number of atoms
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::calc_coordination_number()
{
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    k_cn = DAT::tdual_kkfloat_1d("pair:cn",nmax);
    k_dc6 = DAT::tdual_kkfloat_1d("pair:dc6",nmax);
    d_cn = k_cn.template view<DeviceType>();
    d_dc6 = k_dc6.template view<DeviceType>();
    h_cn = k_cn.view_host();
    h_dc6 = k_dc6.view_host();
  }

  // cn/dc6 are recomputed from scratch; discard prior sync state
  k_cn.clear_sync_state();
  k_dc6.clear_sync_state();
  k_cn.template modify<DeviceType>();
  k_dc6.template modify<DeviceType>();

  // zero out coordination number and dC6

  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType>(0, newton_pair ? nall : nlocal),
    PairDispersionD3InitializeFunctor<DeviceType>{d_cn, d_dc6});

  // calculate coordination number

  if (neighflag == FULL) {
    dispatch_coordination_kernel<FULL>();
  } else if (neighflag == HALFTHREAD) {
    dispatch_coordination_kernel<HALFTHREAD>();
  } else if (neighflag == HALF) {
    dispatch_coordination_kernel<HALF>();
  } else {
    error->all(FLERR, "Must use half or full neighbor list style with pair dispersion/d3/kk");
  }

  // communicate coordination number
  communicationStage = 1;
  if (newton_pair) comm->reverse_comm(this);
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
void PairDispersionD3Kokkos<DeviceType>::dispatch_coordination_kernel()
{
  if (newton_pair) {
    PairDispersionD3CoordinationNumberKernel<DeviceType,NEIGHFLAG,1> cnkernel(
      x, type, d_rcov, d_cn, d_ilist, d_numneigh, d_neighbors, nlocal, cn_thr);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, inum), cnkernel);
    cnkernel.contribute();
  } else {
    PairDispersionD3CoordinationNumberKernel<DeviceType,NEIGHFLAG,0> cnkernel(
      x, type, d_rcov, d_cn, d_ilist, d_numneigh, d_neighbors, nlocal, cn_thr);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, inum), cnkernel);
    cnkernel.contribute();
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairDispersionD3Kokkos<DeviceType>::init_one(int i, int j)
{
  const double cut = PairDispersionD3::init_one(i, j);
  // Since cutsq is written on host by base class we mark the host view modified
  k_cutsq.modify_host();
  return cut;
}

/* ----------------------------------------------------------------------
   init pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::init_style()
{
  PairDispersionD3::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;

  // a full neighbor list visits every pair twice, so the pairwise forces do
  // not add up to the fdotr virial and it has to be tallied explicitly

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   allocate req. arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::allocate()
{
  PairDispersionD3::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_r2r4 = DAT::tdual_kkfloat_1d("pair:r2r4", n+1);
  k_rcov = DAT::tdual_kkfloat_1d("pair:rcov", n+1);
  k_mxci = DAT::tdual_int_1d("pair:mxci", n+1);
  k_r0ab = DAT::tdual_kkfloat_2d("pair:r0ab", n+1, n+1);

  d_r2r4 = k_r2r4.template view<DeviceType>();
  d_rcov = k_rcov.template view<DeviceType>();
  d_mxci = k_mxci.template view<DeviceType>();
  d_r0ab = k_r0ab.template view<DeviceType>();

  // k_c6ab is created in coeff(), where the reference CN grid size is known
}

/* ----------------------------------------------------------------------
   Coeff: read from pair_coeff
          pair_coeff * * path_r0ab.csv path_c6ab.csv functional element1 element2 ...
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  // the base class parses the arguments and fills the host side tables
  PairDispersionD3::coeff(narg, arg);

  const int ntypes = atom->ntypes;

  // now that max_mxci is known, resize the C6 table to the grid actually used

  const int ngrid = max_mxci + 1;
  if ((k_c6ab.extent_int(2) != ngrid) || (k_c6ab.extent_int(3) != ngrid)) {
    k_c6ab = decltype(k_c6ab)("pair:c6ab", ntypes + 1, ntypes + 1, ngrid, ngrid, 3);
    d_c6ab = k_c6ab.template view<DeviceType>();
  }

  // copy coefficients to device

  auto h_r2r4 = k_r2r4.view_host();
  auto h_rcov = k_rcov.view_host();
  auto h_mxci = k_mxci.view_host();
  auto h_r0ab = k_r0ab.view_host();
  auto h_c6ab = k_c6ab.view_host();

  for (int i = 1; i <= ntypes; i++) {
    h_r2r4(i) = r2r4[i];
    h_rcov(i) = rcov[i];
    h_mxci(i) = mxci[i];
  }

  for (int i = 1; i <= ntypes; i++) {
    for (int j = 1; j <= ntypes; j++) {
      h_r0ab(i, j) = r0ab[i][j];
      for (int ci = 0; ci < ngrid; ci++)
        for (int cj = 0; cj < ngrid; cj++)
          for (int k = 0; k < 3; k++) h_c6ab(i, j, ci, cj, k) = c6ab[i][j][ci][cj][k];
    }
  }

  k_r2r4.template modify<LMPHostType>();
  k_rcov.template modify<LMPHostType>();
  k_mxci.template modify<LMPHostType>();
  k_r0ab.template modify<LMPHostType>();
  k_c6ab.template modify<LMPHostType>();

  k_r2r4.template sync<DeviceType>();
  k_rcov.template sync<DeviceType>();
  k_mxci.template sync<DeviceType>();
  k_r0ab.template sync<DeviceType>();
  k_c6ab.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   Compute : energy, force, and stress (Required)
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_mxci.template sync<DeviceType>();
  k_c6ab.template sync<DeviceType>();
  if (eflag || vflag)
    atomKK->modified(execution_space,datamask_modify);
  else
    atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  // must stay in sync with the DUP alias of the kernels below, or the wrong
  // one of the dup_/ndup_ scatter view pair is handed to the functor

  if (neighflag == FULL) {
    need_dup = std::is_same_v<NeedDup_v<FULL,DeviceType>, Kokkos::Experimental::ScatterDuplicated>;
  } else if (neighflag == HALFTHREAD) {
    need_dup = std::is_same_v<NeedDup_v<HALFTHREAD,DeviceType>, Kokkos::Experimental::ScatterDuplicated>;
  } else {
    need_dup = std::is_same_v<NeedDup_v<HALF,DeviceType>, Kokkos::Experimental::ScatterDuplicated>;
  }

  EV_FLOAT ev_all = {};

  // Calculate the coordination number for each atom
  calc_coordination_number();

  // Since communication impacts host view we refresh device view before the device kernels
  k_cn.template sync<DeviceType>();

  // Could move these scatter view manipulations into a separate function
  dup_eatom = {};
  dup_vatom = {};
  ndup_eatom = {};
  ndup_vatom = {};
  dup_dc6 = {};
  ndup_dc6 = {};
  if (need_dup) {
    dup_f   = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
  }
  ndup_f  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
  if (need_dup) {
    dup_dc6  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_dc6);
  }
  ndup_dc6 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_dc6);

  if (eflag_atom) {
    if (need_dup) {
      dup_eatom  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    }
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
  }
  if (vflag_atom) {
    if (need_dup) {
      dup_vatom  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
    }
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  // first: compute forces, energies, and dC6
  EV_FLOAT ev = {};

  if (neighflag == HALF) {
    dispatch_kernel_A<HALF>(ev);
  } else if (neighflag == HALFTHREAD) {
    dispatch_kernel_A<HALFTHREAD>(ev);
  } else if (neighflag == FULL) {
    dispatch_kernel_A<FULL>(ev);
  } else {
    error->all(FLERR, "Must use half or full neighbor list style with pair dispersion/d3/kk");
  }

  if (evflag) ev_all += ev;

  if (need_dup)
    Kokkos::Experimental::contribute(d_dc6, dup_dc6);
  else
    Kokkos::Experimental::contribute(d_dc6, ndup_dc6);
  k_dc6.template modify<DeviceType>();

  // communicate derivatives of C6
  communicationStage = 2;
  if (newton_pair) comm->reverse_comm(this);
  comm->forward_comm(this);

  // Communication updates host view so we must refresh device view before kernel B
  k_dc6.template sync<DeviceType>();

  // second: compute and apply force contribution from dC6
  ev = {};

  if (neighflag == HALF) {
    dispatch_kernel_B<HALF>(ev);
  } else if (neighflag == HALFTHREAD) {
    dispatch_kernel_B<HALFTHREAD>(ev);
  } else if (neighflag == FULL) {
    dispatch_kernel_B<FULL>(ev);
  } else {
    error->all(FLERR, "Must use half or full neighbor list style with pair dispersion/d3/kk");
  }

  if (evflag) ev_all += ev;

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);
  else
    Kokkos::Experimental::contribute(f, ndup_f);

  if (eflag_global) eng_vdwl += ev_all.evdwl;
  if (vflag_global) {
    virial[0] += ev_all.v[0];
    virial[1] += ev_all.v[1];
    virial[2] += ev_all.v[2];
    virial[3] += ev_all.v[3];
    virial[4] += ev_all.v[4];
    virial[5] += ev_all.v[5];
  }

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    else
      Kokkos::Experimental::contribute(d_eatom, ndup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    else
      Kokkos::Experimental::contribute(d_vatom, ndup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  // free duplicated memory
  if (need_dup) {
    dup_f     = {};
    dup_eatom = {};
    dup_vatom = {};
    dup_dc6   = {};
  }
}

/* ----------------------------------------------------------------------
   Templated kernel methods
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
void PairDispersionD3Kokkos<DeviceType>::launch_kernel_A(EV_FLOAT &ev)
{
  auto functor = PairDispersionD3KernelA<DeviceType, NEIGHFLAG, NEWTON_PAIR, EVFLAG>(
      x, type,
      d_cutsq, d_cn, d_dc6,
      d_r2r4, d_r0ab, d_c6ab, d_mxci,
      d_numneigh, d_neighbors, d_ilist,
      dup_f, ndup_f,
      dup_eatom, ndup_eatom,
      dup_vatom, ndup_vatom,
      dup_dc6, ndup_dc6,
      special_lj,
      nlocal,
      eflag, vflag_either,
      eflag_global, eflag_atom,
      vflag_global, vflag_atom,
      dampingCode,
      s6, s8, rs6, rs8,
      a1, a2, alpha);

  if constexpr (EVFLAG) {
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType>(0, inum),
        functor,
        ev);
  } else {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType>(0, inum),
        functor);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
void PairDispersionD3Kokkos<DeviceType>::dispatch_kernel_A(EV_FLOAT &ev)
{
  if (newton_pair) {
    if (evflag) launch_kernel_A<NEIGHFLAG, 1, 1>(ev);
    else        launch_kernel_A<NEIGHFLAG, 1, 0>(ev);
  } else {
    if (evflag) launch_kernel_A<NEIGHFLAG, 0, 1>(ev);
    else        launch_kernel_A<NEIGHFLAG, 0, 0>(ev);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
void PairDispersionD3Kokkos<DeviceType>::launch_kernel_B(EV_FLOAT &ev)
{
  auto functor = PairDispersionD3KernelB<DeviceType, NEIGHFLAG, NEWTON_PAIR, EVFLAG>(
      x, type,
      d_cutsq, d_dc6, d_rcov,
      d_numneigh, d_neighbors, d_ilist,
      dup_f, ndup_f,
      dup_eatom, ndup_eatom,
      dup_vatom, ndup_vatom,
      special_lj,
      nlocal,
      eflag, vflag_either,
      eflag_global, eflag_atom,
      vflag_global, vflag_atom,
      cn_thr);

  if constexpr (EVFLAG) {
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType>(0, inum),
        functor,
        ev);
  } else {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<DeviceType>(0, inum),
        functor);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
void PairDispersionD3Kokkos<DeviceType>::dispatch_kernel_B(EV_FLOAT &ev)
{
  if (newton_pair) {
    if (evflag) launch_kernel_B<NEIGHFLAG, 1, 1>(ev);
    else        launch_kernel_B<NEIGHFLAG, 1, 0>(ev);
  } else {
    if (evflag) launch_kernel_B<NEIGHFLAG, 0, 1>(ev);
    else        launch_kernel_B<NEIGHFLAG, 0, 0>(ev);
  }
}

/* ----------------------------------------------------------------------
   Communication section
------------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_double_1d &buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  v_buf = buf.view<DeviceType>();
  if (communicationStage == 1) {
    k_cn.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3PackForwardCommFunctor<DeviceType>{d_cn, d_sendlist, v_buf});
  }
  if (communicationStage == 2) {
    k_dc6.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3PackForwardCommFunctor<DeviceType>{d_dc6, d_sendlist, v_buf});
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first, DAT::tdual_double_1d &buf)
{
  v_buf = buf.view<DeviceType>();
  if (communicationStage == 1) {
    k_cn.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3UnpackForwardCommFunctor<DeviceType>{d_cn, first, v_buf});
    k_cn.template modify<DeviceType>();
  }
  if (communicationStage == 2) {
    k_dc6.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3UnpackForwardCommFunctor<DeviceType>{d_dc6, first, v_buf});
    k_dc6.template modify<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first, DAT::tdual_double_1d &buf)
{
  v_buf = buf.view<DeviceType>();
  if (communicationStage == 1) {
    k_cn.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3PackReverseCommFunctor<DeviceType>{d_cn, first, v_buf});
  }
  if (communicationStage == 2) {
    k_dc6.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3PackReverseCommFunctor<DeviceType>{d_dc6, first, v_buf});
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_1d k_recvlist, DAT::tdual_double_1d &buf)
{
  d_recvlist = k_recvlist.view<DeviceType>();
  v_buf = buf.view<DeviceType>();
  if (communicationStage == 1) {
    k_cn.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3UnpackReverseCommFunctor<DeviceType>{d_cn, d_recvlist, v_buf});
    k_cn.template modify<DeviceType>();
  }
  if (communicationStage == 2) {
    k_dc6.template sync<DeviceType>();
    Kokkos::parallel_for(
      Kokkos::RangePolicy<DeviceType>(0, n),
      PairDispersionD3UnpackReverseCommFunctor<DeviceType>{d_dc6, d_recvlist, v_buf});
    k_dc6.template modify<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  if (communicationStage == 1) {

    k_cn.sync_host();

    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_cn(j);
    }
  }
  if (communicationStage == 2) {

    k_dc6.sync_host();

    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = h_dc6(j);
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (communicationStage == 1) {
    k_cn.sync_host();

    for (i = first; i < last; i++) { h_cn(i) = buf[m++]; }

    k_cn.modify_host();
  }
  if (communicationStage == 2) {

    k_dc6.sync_host();

    for (i = first; i < last; i++) { h_dc6(i) = buf[m++]; }

    k_dc6.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  if (communicationStage == 1) {
    k_cn.sync_host();
    for (i = first; i < last; i++) { buf[m++] = h_cn(i); }
  }
  if (communicationStage == 2) {
    k_dc6.sync_host();
    for (i = first; i < last; i++) { buf[m++] = h_dc6(i); }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  if (communicationStage == 1) {
    k_cn.sync_host();
    for (i = 0; i < n; i++) {
      j = list[i];
      h_cn(j) += buf[m++];
    }
    k_cn.modify_host();
  }
  if (communicationStage == 2) {
    k_dc6.sync_host();
    for (i = 0; i < n; i++) {
      j = list[i];
      h_dc6(j) += buf[m++];
    }
    k_dc6.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairDispersionD3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairDispersionD3Kokkos<LMPHostType>;
#endif
}
