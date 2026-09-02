/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_oxdna2_coaxstk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "fix_oxdna_npair_kokkos.h"
#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;
using MathConst::MY_PI;

namespace {
constexpr KK_FLOAT MY_PI_KK = static_cast<KK_FLOAT>(MY_PI);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna2CoaxstkKokkos<DeviceType>::PairOxdna2CoaxstkKokkos(LAMMPS *lmp) : PairOxdna2Coaxstk(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  oxdnaflag = EnabledOXDNAFlag::OXDNA2;
  fix_oxdna_lrfKK = nullptr;
  fix_oxdna_npairKK = nullptr;
  screened_pair_count = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna2CoaxstkKokkos<DeviceType>::~PairOxdna2CoaxstkKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2CoaxstkKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);

  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  torque = atomKK->k_torque.template view<DeviceType>();
  type = atomKK->k_type.template view<DeviceType>();
  id5p = atomKK->k_id5p.template view<DeviceType>();
  id3p = atomKK->k_id3p.template view<DeviceType>();

  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  // get the neighbor list and neighbors used in operator()

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_neighbors = k_list->d_neighbors;
  anum = list->inum;
  d_alist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;

  int need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterDuplicated>(f);
    dup_torque = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterDuplicated>(torque);
  } else {
    ndup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_torque = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterNonDuplicated>(torque);
  }

  copymode = 1;

  // d_n(x/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->k_nx.template view<DeviceType>();
  d_ny_xtrct = fix_oxdna_lrfKK->k_ny.template view<DeviceType>();
  d_nz_xtrct = fix_oxdna_lrfKK->k_nz.template view<DeviceType>();

  // If we're on a GPU, look up fix_oxdna_npairKK screened pair count and packed pair view.
  if (execution_space != HostKK) {
    screened_pair_count = fix_oxdna_npairKK->screened_pair_count;
    d_pairs_screened = fix_oxdna_npairKK->k_pairs_screened.template view<DeviceType>();
  }

  // loop over neighbors of my atoms for compute functors

  EV_FLOAT ev;

  auto run_compute_host = [&](auto host_tag, auto evflag_tag) {
    if constexpr (decltype(evflag_tag)::value) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, decltype(host_tag)>(0,anum), *this, ev);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, decltype(host_tag)>(0,anum), *this);
    }
  };

  auto run_compute_gpu = [&](auto gpu_tag, auto evflag_tag) {
    if constexpr (decltype(evflag_tag)::value) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, decltype(gpu_tag)>(0,screened_pair_count), *this, ev);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, decltype(gpu_tag)>(0,screened_pair_count), *this);
    }
  };

  auto run_compute_by_oxdnaflag = [&](auto oxdnaflag_tag, auto neighflag_tag, auto newtonpair_tag, auto evflag_tag) {
    if (execution_space == HostKK) {
      run_compute_host(TagPairOxdna2CoaxstkCompute<oxdnaflag_tag.value,neighflag_tag.value,newtonpair_tag.value,evflag_tag.value>{}, evflag_tag);
    } else {
      run_compute_gpu(TagPairOxdna2CoaxstkComputeGPUPair<oxdnaflag_tag.value,neighflag_tag.value,newtonpair_tag.value,evflag_tag.value>{}, evflag_tag);
    }
  };

  const int dispatch_neigh = (neighflag == HALF) ? 0 : (neighflag == HALFTHREAD) ? 1 : 2;
  const int dispatch_key = (evflag ? 8 : 0) | (newton_pair ? 4 : 0) | dispatch_neigh;

  switch (oxdnaflag) {
    case EnabledOXDNAFlag::OXDNA2:
      switch (dispatch_key) {
        case 0: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
        case 1: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
        case 2: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,FULL>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
        case 4: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALF>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
        case 5: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
        case 6: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,FULL>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
        case 8: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
        case 9: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
        case 10: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,FULL>{},      std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
        case 12: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALF>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
        case 13: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,HALFTHREAD>{},std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
        case 14: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA2>{}, std::integral_constant<int,FULL>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
        default: error->all(FLERR,"Illegal pair_style command");
      }
      break;
    case EnabledOXDNAFlag::OXDNA3:
      switch (dispatch_key) {
        case 0: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
        case 1: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
        case 2: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,FULL>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
        case 4: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALF>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
        case 5: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
        case 6: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,FULL>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
        case 8: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
        case 9: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
        case 10: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,FULL>{},      std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
        case 12: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALF>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
        case 13: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,HALFTHREAD>{},std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
        case 14: run_compute_by_oxdnaflag(std::integral_constant<int,OXDNA3>{}, std::integral_constant<int,FULL>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
        default: error->all(FLERR,"Illegal pair_style command");
      }
      break;
    default:
      error->all(FLERR,"Illegal pair_style command");
  }

  if (need_dup) {
    Kokkos::Experimental::contribute(f, dup_f);
    Kokkos::Experimental::contribute(torque, dup_torque);
  }

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f        = decltype(dup_f)();
    dup_torque   = decltype(dup_torque)();
    dup_eatom    = decltype(dup_eatom)();
    dup_vatom    = decltype(dup_vatom)();
  }
}

/* ----------------------------------------------------------------------
   Standard non-GPU Compute Functor(s)
-------------------------------------------------------------------------- */

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::operator()(TagPairOxdna2CoaxstkCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia, EV_FLOAT &ev) const
{
  // f and torque array are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int a = d_alist(ia);
  const int atype = type(a);
  // vectors COM-backbone site, COM-stacking site in lab frame
  KK_FLOAT ra_cbk[3], rb_cbk[3], ra_cstk[3], rb_cstk[3];

  KK_ACC_FLOAT delf[3],delta[3],deltb[3];    // force, torque increment
  KK_ACC_FLOAT evdwl, finc, tpair;           // energy, force, torque
  KK_FLOAT v1tmp[3];
  KK_FLOAT delr_bkbk[3],delr_bkbk_norm[3],rsq_bkbk,r_bkbk,rinv_bkbk;
  KK_FLOAT delr_stkstk[3],delr_stkstk_norm[3],rsq_stkstk,r_stkstk,rinv_stkstk;
  KK_FLOAT theta1,theta1p,t1dir[3],cost1;
  KK_FLOAT theta4,t4dir[3],cost4;
  KK_FLOAT theta5,theta5p,t5dir[3],cost5;
  KK_FLOAT theta6,theta6p,t6dir[3],cost6;
  KK_FLOAT cosphi3;

  KK_FLOAT f2,f4f6t1,f4t4,f4t5,f4t6;
  KK_FLOAT prime_cxst_ab;
  KK_FLOAT df2,df4f6t1,df4t4,df4t5,df4t6;

  // a has to be terminal nucleotide
  // NOTE: this is the same as 'continue' in vanilla code, but since within the functor this is effectively the outer a-loop,
  // we just return to skip the rest of the functor for this a.
  if (id3p(a) != -1 && id5p(a) != -1) return;

  // vector COM-backbone site a, COM-stacking site a
  if constexpr (OXDNAFLAG==OXDNA2) {
    constexpr KK_FLOAT dx_cbk_oxdna2 = -0.34;
    constexpr KK_FLOAT dy_cbk_oxdna2 = +0.3408;
    ra_cbk[0] = dx_cbk_oxdna2*d_nx_xtrct(a,0) + dy_cbk_oxdna2*d_ny_xtrct(a,0);
    ra_cbk[1] = dx_cbk_oxdna2*d_nx_xtrct(a,1) + dy_cbk_oxdna2*d_ny_xtrct(a,1);
    ra_cbk[2] = dx_cbk_oxdna2*d_nx_xtrct(a,2) + dy_cbk_oxdna2*d_ny_xtrct(a,2);
    constexpr KK_FLOAT dx_cstk_oxdna1 = +0.34;  // oxDNA2 uses same stacking site as oxDNA1
    ra_cstk[0] = dx_cstk_oxdna1*d_nx_xtrct(a,0);
    ra_cstk[1] = dx_cstk_oxdna1*d_nx_xtrct(a,1);
    ra_cstk[2] = dx_cstk_oxdna1*d_nx_xtrct(a,2);
  } else if constexpr (OXDNAFLAG==OXDNA3) {
    // oxDNA3 uses same backbone site as oxDNA2...
    constexpr KK_FLOAT dx_cbk_oxdna2 = -0.34;
    constexpr KK_FLOAT dy_cbk_oxdna2 = +0.3408;
    ra_cbk[0] = dx_cbk_oxdna2*d_nx_xtrct(a,0) + dy_cbk_oxdna2*d_ny_xtrct(a,0);
    ra_cbk[1] = dx_cbk_oxdna2*d_nx_xtrct(a,1) + dy_cbk_oxdna2*d_ny_xtrct(a,1);
    ra_cbk[2] = dx_cbk_oxdna2*d_nx_xtrct(a,2) + dy_cbk_oxdna2*d_ny_xtrct(a,2);
    // ...But the stacking site is different for oxDNA3.
    constexpr KK_FLOAT dx_cstk_oxdna3 = +0.37;
    ra_cstk[0] = dx_cstk_oxdna3*d_nx_xtrct(a,0);
    ra_cstk[1] = dx_cstk_oxdna3*d_nx_xtrct(a,1);
    ra_cstk[2] = dx_cstk_oxdna3*d_nx_xtrct(a,2);
  }

  const int bnum = d_numneigh(a);

  for (int ib = 0; ib < bnum; ib++) {

    int b = d_neighbors(a,ib);
    const KK_FLOAT factor_lj = special_lj[sbmask(b)];
    b &= NEIGHMASK;
    const int btype = type(b);

    // b has to be terminal nucleotide
    if(id3p(b)!=-1 && id5p(b)!=-1) continue;

    // vector COM-backbone site b, COM-stacking site b
    if constexpr (OXDNAFLAG==OXDNA2) {
      constexpr KK_FLOAT dx_cbk_oxdna2 = -0.34;
      constexpr KK_FLOAT dy_cbk_oxdna2 = +0.3408;
      rb_cbk[0] = dx_cbk_oxdna2*d_nx_xtrct(b,0) + dy_cbk_oxdna2*d_ny_xtrct(b,0);
      rb_cbk[1] = dx_cbk_oxdna2*d_nx_xtrct(b,1) + dy_cbk_oxdna2*d_ny_xtrct(b,1);
      rb_cbk[2] = dx_cbk_oxdna2*d_nx_xtrct(b,2) + dy_cbk_oxdna2*d_ny_xtrct(b,2);
      constexpr KK_FLOAT dx_cstk_oxdna1 = +0.34;  // oxDNA2 uses same stacking site as oxDNA1
      rb_cstk[0] = dx_cstk_oxdna1*d_nx_xtrct(b,0);
      rb_cstk[1] = dx_cstk_oxdna1*d_nx_xtrct(b,1);
      rb_cstk[2] = dx_cstk_oxdna1*d_nx_xtrct(b,2);
    } else if constexpr (OXDNAFLAG==OXDNA3) {
      // oxDNA3 uses same backbone site as oxDNA2...
      constexpr KK_FLOAT dx_cbk_oxdna2 = -0.34;
      constexpr KK_FLOAT dy_cbk_oxdna2 = +0.3408;
      rb_cbk[0] = dx_cbk_oxdna2*d_nx_xtrct(b,0) + dy_cbk_oxdna2*d_ny_xtrct(b,0);
      rb_cbk[1] = dx_cbk_oxdna2*d_nx_xtrct(b,1) + dy_cbk_oxdna2*d_ny_xtrct(b,1);
      rb_cbk[2] = dx_cbk_oxdna2*d_nx_xtrct(b,2) + dy_cbk_oxdna2*d_ny_xtrct(b,2);
      // ...But the stacking site is different for oxDNA3.
      constexpr KK_FLOAT dx_cstk_oxdna3 = +0.37;
      rb_cstk[0] = dx_cstk_oxdna3*d_nx_xtrct(b,0);
      rb_cstk[1] = dx_cstk_oxdna3*d_nx_xtrct(b,1);
      rb_cstk[2] = dx_cstk_oxdna3*d_nx_xtrct(b,2);
    }

    // vector stacking site b to a
    delr_stkstk[0] = x(a,0) + ra_cstk[0] - x(b,0) - rb_cstk[0];
    delr_stkstk[1] = x(a,1) + ra_cstk[1] - x(b,1) - rb_cstk[1];
    delr_stkstk[2] = x(a,2) + ra_cstk[2] - x(b,2) - rb_cstk[2];

    rsq_stkstk = delr_stkstk[0]*delr_stkstk[0] + delr_stkstk[1]*delr_stkstk[1] + delr_stkstk[2]*delr_stkstk[2];
    r_stkstk = sqrtf(rsq_stkstk);
    rinv_stkstk = 1.0 / r_stkstk;

    delr_stkstk_norm[0] = delr_stkstk[0] * rinv_stkstk;
    delr_stkstk_norm[1] = delr_stkstk[1] * rinv_stkstk;
    delr_stkstk_norm[2] = delr_stkstk[2] * rinv_stkstk;

    // vector backbone site b to a
    delr_bkbk[0] = x(a,0) + ra_cbk[0] - x(b,0) - rb_cbk[0];
    delr_bkbk[1] = x(a,1) + ra_cbk[1] - x(b,1) - rb_cbk[1];
    delr_bkbk[2] = x(a,2) + ra_cbk[2] - x(b,2) - rb_cbk[2];

    rsq_bkbk = delr_bkbk[0]*delr_bkbk[0] + delr_bkbk[1]*delr_bkbk[1] + delr_bkbk[2]*delr_bkbk[2];
    r_bkbk = sqrtf(rsq_bkbk);
    rinv_bkbk = 1.0 / r_bkbk;

    delr_bkbk_norm[0] = delr_bkbk[0] * rinv_bkbk;
    delr_bkbk_norm[1] = delr_bkbk[1] * rinv_bkbk;
    delr_bkbk_norm[2] = delr_bkbk[2] * rinv_bkbk;

    cost1 = -(d_nx_xtrct(a,0) * d_nx_xtrct(b,0) + d_nx_xtrct(a,1) * d_nx_xtrct(b,1) + d_nx_xtrct(a,2) * d_nx_xtrct(b,2));
    if (cost1 >  1.0) cost1 =  1.0;
    if (cost1 < -1.0) cost1 = -1.0;
    theta1 = acos(cost1);
    theta1p = 2 * MY_PI - theta1;

    // beginning of modulation factors

    // f4f6t1 = f4(theta1,..) + f6(theta1,..) modulation factors
    f4f6t1 = F4_KK(theta1, d_a_cxst1(atype,btype), d_theta_cxst1_0(atype,btype), d_dtheta_cxst1_ast(atype,btype),
                 d_b_cxst1(atype,btype), d_dtheta_cxst1_c(atype,btype)) + \
           F6_KK(theta1, d_AA_cxst1(atype,btype), d_BB_cxst1(atype,btype));

    // start early rejection criterium
    if (f4f6t1 != static_cast<KK_FLOAT>(0.0)) {
      // theta4 calculation
      cost4 = d_nz_xtrct(a,0)*d_nz_xtrct(b,0) + d_nz_xtrct(a,1)*d_nz_xtrct(b,1) + d_nz_xtrct(a,2)*d_nz_xtrct(b,2);
      if (cost4 > 1.0) cost4 = 1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);
      // f4t4 = f4 modulation factor
      f4t4 = F4_KK(theta4, d_a_cxst4(atype,btype), d_theta_cxst4_0(atype, btype), d_dtheta_cxst4_ast(atype, btype),
              d_b_cxst4(atype, btype), d_dtheta_cxst4_c(atype, btype)) +
             F4_KK(theta4, d_a_cxst4(atype,btype), MY_PI - d_theta_cxst4_0(atype, btype), d_dtheta_cxst4_ast(atype, btype),
              d_b_cxst4(atype, btype), d_dtheta_cxst4_c(atype, btype));
    // end of f4f6t1

    // f4t4 early rejection criterium
    if (f4t4 != static_cast<KK_FLOAT>(0.0)) {
      cost5 = (d_nz_xtrct(a,0)*delr_stkstk_norm[0] + d_nz_xtrct(a,1)*delr_stkstk_norm[1] + d_nz_xtrct(a,2)*delr_stkstk_norm[2]);
      if (cost5 > 1.0) cost5 = 1.0;
      if (cost5 < -1.0) cost5 = -1.0;
      theta5 = acos(cost5);
      theta5p = MY_PI - theta5;
      // f4t5 = f4(theta5,..) + f4(theta5p,..) modulation factors
      f4t5 = F4_KK(theta5, d_a_cxst5(atype,btype), d_theta_cxst5_0(atype,btype), d_dtheta_cxst5_ast(atype,btype),
              d_b_cxst5(atype,btype), d_dtheta_cxst5_c(atype,btype)) + \
             F4_KK(theta5p, d_a_cxst5(atype,btype), d_theta_cxst5_0(atype,btype), d_dtheta_cxst5_ast(atype,btype),
              d_b_cxst5(atype,btype), d_dtheta_cxst5_c(atype,btype));
    // end of f4t4

    // f4t5 early rejection criterium
    if (f4t5 != static_cast<KK_FLOAT>(0.0)) {
      cost6 = d_nz_xtrct(b,0)*delr_stkstk_norm[0] + d_nz_xtrct(b,1)*delr_stkstk_norm[1] + d_nz_xtrct(b,2)*delr_stkstk_norm[2];
      if (cost6 > 1.0) cost6 = 1.0;
      if (cost6 < -1.0) cost6 = -1.0;
      theta6 = acos(cost6);
      theta6p = MY_PI - theta6;
      // f4t6 = f4(theta6,..) + f4(theta6p,..) modulation factors
      f4t6 = F4_KK(theta6, d_a_cxst6(atype,btype), d_theta_cxst6_0(atype,btype), d_dtheta_cxst6_ast(atype,btype),
              d_b_cxst6(atype,btype), d_dtheta_cxst6_c(atype,btype)) + \
             F4_KK(theta6p, d_a_cxst6(atype,btype), d_theta_cxst6_0(atype,btype), d_dtheta_cxst6_ast(atype,btype),
              d_b_cxst6(atype,btype), d_dtheta_cxst6_c(atype,btype));

      v1tmp[0] = delr_bkbk_norm[1] * d_nx_xtrct(a,2) - delr_bkbk_norm[2] * d_nx_xtrct(a,1);
      v1tmp[1] = delr_bkbk_norm[2] * d_nx_xtrct(a,0) - delr_bkbk_norm[0] * d_nx_xtrct(a,2);
      v1tmp[2] = delr_bkbk_norm[0] * d_nx_xtrct(a,1) - delr_bkbk_norm[1] * d_nx_xtrct(a,0);
      cosphi3 = v1tmp[0] * delr_stkstk_norm[0] + v1tmp[1] * delr_stkstk_norm[1] + v1tmp[2] * delr_stkstk_norm[2];
      if (cosphi3 > 1.0) cosphi3 = 1.0;
      if (cosphi3 < -1.0) cosphi3 = -1.0;

      // Direction-dependent coaxial stacking strength.
      if (id5p(a) == -1 && id3p(b) == -1) {
        prime_cxst_ab = d_k_cxst(atype,btype);
      } else if (id3p(a) == -1 && id5p(b) == -1) {
        prime_cxst_ab = d_k_cxst(btype,atype);
      } else {
        prime_cxst_ab = 0.5 * (d_k_cxst(atype,btype) + d_k_cxst(btype,atype));
      }

      // f2 = f2 modulation factor
      f2 = F2_KK(r_stkstk, prime_cxst_ab, d_cut_cxst_0(atype,btype), d_cut_cxst_lc(atype,btype),
              d_cut_cxst_hc(atype,btype), d_cut_cxst_lo(atype,btype), d_cut_cxst_hi(atype,btype),
              d_b_cxst_lo(atype,btype), d_b_cxst_hi(atype,btype),
              d_cut_cxst_c(atype,btype));

      evdwl = f2 * f4f6t1 * f4t4 * f4t5 * f4t6 * factor_lj;
    // end of f4t5

    // evdwl early rejection criterium
    if (evdwl != static_cast<KK_FLOAT>(0.0)) {
      // df2 = DF2 modulation factor
      df2 = DF2_KK(r_stkstk, prime_cxst_ab, d_cut_cxst_0(atype,btype), d_cut_cxst_lc(atype,btype),
              d_cut_cxst_hc(atype,btype), d_cut_cxst_lo(atype,btype), d_cut_cxst_hi(atype,btype),
              d_b_cxst_lo(atype,btype), d_b_cxst_hi(atype,btype));
      // df4f6t1 = DF4(theta1,..)/sin(theta1) + DF6(theta1,..)/sin(theta1) modulation factors
      df4f6t1 = ( DF4_KK(theta1, d_a_cxst1(atype,btype), d_theta_cxst1_0(atype,btype), d_dtheta_cxst1_ast(atype,btype),
                     d_b_cxst1(atype,btype), d_dtheta_cxst1_c(atype,btype)) + \
              DF6_KK(theta1, d_AA_cxst1(atype,btype), d_BB_cxst1(atype,btype)) ) / sin(theta1);
      // df4t4 = DF4(theta4,..)/sin(theta4) + DF4(theta4, mirrored theta0)/sin(theta4)
      df4t4 = ( DF4_KK(theta4, d_a_cxst4(atype,btype), d_theta_cxst4_0(atype, btype), d_dtheta_cxst4_ast(atype, btype),
              d_b_cxst4(atype, btype), d_dtheta_cxst4_c(atype, btype)) +
            DF4_KK(theta4, d_a_cxst4(atype,btype), MY_PI - d_theta_cxst4_0(atype, btype), d_dtheta_cxst4_ast(atype, btype),
              d_b_cxst4(atype, btype), d_dtheta_cxst4_c(atype, btype)) ) / sin(theta4);
      // df4t5 = DF4(theta5,..)/sin(theta5) - DF4(theta5p,..)/sin(theta5) modulation factors
      df4t5 = ( DF4_KK(theta5, d_a_cxst5(atype,btype), d_theta_cxst5_0(atype,btype), d_dtheta_cxst5_ast(atype,btype),
                     d_b_cxst5(atype,btype), d_dtheta_cxst5_c(atype,btype)) - \
              DF4_KK(theta5p, d_a_cxst5(atype,btype), d_theta_cxst5_0(atype,btype), d_dtheta_cxst5_ast(atype,btype),
                     d_b_cxst5(atype,btype), d_dtheta_cxst5_c(atype,btype)) ) / sin(theta5);
      // df4t6 = DF4(theta6,..)/sin(theta6) - DF4(theta6p,..)/sin(theta6) modulation factors
      df4t6 = ( DF4_KK(theta6, d_a_cxst6(atype,btype), d_theta_cxst6_0(atype,btype), d_dtheta_cxst6_ast(atype,btype),
                     d_b_cxst6(atype,btype), d_dtheta_cxst6_c(atype,btype)) - \
              DF4_KK(theta6p, d_a_cxst6(atype,btype), d_theta_cxst6_0(atype,btype), d_dtheta_cxst6_ast(atype,btype),
                     d_b_cxst6(atype,btype), d_dtheta_cxst6_c(atype,btype)) ) / sin(theta6);

      // force, torque, and viral contributions for forces between h-bonding sites

      delf[0] = 0.0;
      delf[1] = 0.0;
      delf[2] = 0.0;

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;

      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // radial force
      finc  = -df2 * f4f6t1 * f4t4 * f4t5 * f4t6 * rinv_stkstk * factor_lj;

      delf[0] += delr_stkstk[0] * finc;
      delf[1] += delr_stkstk[1] * finc;
      delf[2] += delr_stkstk[2] * finc;

      // theta5 force
      if (theta5 != static_cast<KK_FLOAT>(0.0) && theta5p != static_cast<KK_FLOAT>(0.0)) {

        finc  = -f2 * f4f6t1 * f4t4 * df4t5 * f4t6 * rinv_stkstk * factor_lj;

        delf[0] += (delr_stkstk_norm[0]*cost5 - d_nz_xtrct(a,0)) * finc;
        delf[1] += (delr_stkstk_norm[1]*cost5 - d_nz_xtrct(a,1)) * finc;
        delf[2] += (delr_stkstk_norm[2]*cost5 - d_nz_xtrct(a,2)) * finc;
      }

      // theta6 force
      if (theta6 != static_cast<KK_FLOAT>(0.0) && theta6p != static_cast<KK_FLOAT>(0.0)) {

        finc  = -f2 * f4f6t1* f4t4 * f4t5 * df4t6 * rinv_stkstk * factor_lj;

        delf[0] += (delr_stkstk_norm[0]*cost6 - d_nz_xtrct(b,0)) * finc;
        delf[1] += (delr_stkstk_norm[1]*cost6 - d_nz_xtrct(b,1)) * finc;
        delf[2] += (delr_stkstk_norm[2]*cost6 - d_nz_xtrct(b,2)) * finc;
      }

      // increment forces and torques

      a_f(a,0) += delf[0];
      a_f(a,1) += delf[1];
      a_f(a,2) += delf[2];
      delta[0] = ra_cstk[1]*delf[2] - ra_cstk[2]*delf[1];
      delta[1] = ra_cstk[2]*delf[0] - ra_cstk[0]*delf[2];
      delta[2] = ra_cstk[0]*delf[1] - ra_cstk[1]*delf[0];
      a_torque(a,0) += delta[0];
      a_torque(a,1) += delta[1];
      a_torque(a,2) += delta[2];

      if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
        a_f(b,0) -= delf[0];
        a_f(b,1) -= delf[1];
        a_f(b,2) -= delf[2];
        deltb[0] = rb_cstk[1]*delf[2] - rb_cstk[2]*delf[1];
        deltb[1] = rb_cstk[2]*delf[0] - rb_cstk[0]*delf[2];
        deltb[2] = rb_cstk[0]*delf[1] - rb_cstk[1]*delf[0];
        a_torque(b,0) -= deltb[0];
        a_torque(b,1) -= deltb[1];
        a_torque(b,2) -= deltb[2];
      }

      // increment energy and virial
      // NOTE: The virial is calculated on the 'molecular' basis.
      // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

      if (EVFLAG) {
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;

        if (vflag_either || eflag_atom) {
          this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,evdwl,\
          delf[0],delf[1],delf[2],x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2));
        }
      }

      // pure torques not expressible as r x f

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;
      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // theta1 torque
      if (theta1 != static_cast<KK_FLOAT>(0.0) && theta1p != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f2 * df4f6t1 * f4t4 * f4t5 * f4t6 * factor_lj;

        t1dir[0] = d_nx_xtrct(a,1) * d_nx_xtrct(b,2) - d_nx_xtrct(a,2) * d_nx_xtrct(b,1);
        t1dir[1] = d_nx_xtrct(a,2) * d_nx_xtrct(b,0) - d_nx_xtrct(a,0) * d_nx_xtrct(b,2);
        t1dir[2] = d_nx_xtrct(a,0) * d_nx_xtrct(b,1) - d_nx_xtrct(a,1) * d_nx_xtrct(b,0);
        delta[0] += t1dir[0] * tpair;
        delta[1] += t1dir[1] * tpair;
        delta[2] += t1dir[2] * tpair;
        deltb[0] += t1dir[0] * tpair;
        deltb[1] += t1dir[1] * tpair;
        deltb[2] += t1dir[2] * tpair;
      }
      //theta4 torque
      if (theta4 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f2 * f4f6t1 * df4t4 * f4t5 * f4t6 * factor_lj;

        t4dir[0] = d_nz_xtrct(b,1) * d_nz_xtrct(a,2) - d_nz_xtrct(b,2) * d_nz_xtrct(a,1);
        t4dir[1] = d_nz_xtrct(b,2) * d_nz_xtrct(a,0) - d_nz_xtrct(b,0) * d_nz_xtrct(a,2);
        t4dir[2] = d_nz_xtrct(b,0) * d_nz_xtrct(a,1) - d_nz_xtrct(b,1) * d_nz_xtrct(a,0);
        delta[0] += t4dir[0] * tpair;
        delta[1] += t4dir[1] * tpair;
        delta[2] += t4dir[2] * tpair;
        deltb[0] += t4dir[0] * tpair;
        deltb[1] += t4dir[1] * tpair;
        deltb[2] += t4dir[2] * tpair;
      }
      //theta5 torque
      if (theta5 != static_cast<KK_FLOAT>(0.0) && theta5p != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f2 * f4f6t1 * f4t4 * df4t5 * f4t6 * factor_lj;

        t5dir[0] = delr_stkstk_norm[1] * d_nz_xtrct(a,2) - delr_stkstk_norm[2] * d_nz_xtrct(a,1);
        t5dir[1] = delr_stkstk_norm[2] * d_nz_xtrct(a,0) - delr_stkstk_norm[0] * d_nz_xtrct(a,2);
        t5dir[2] = delr_stkstk_norm[0] * d_nz_xtrct(a,1) - delr_stkstk_norm[1] * d_nz_xtrct(a,0);
        delta[0] += t5dir[0] * tpair;
        delta[1] += t5dir[1] * tpair;
        delta[2] += t5dir[2] * tpair;
      }
      // theta6 torque
      if (theta6 != static_cast<KK_FLOAT>(0.0) && theta6p != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f2 * f4f6t1 * f4t4 * f4t5 * df4t6 * factor_lj;

        t6dir[0] = delr_stkstk_norm[1] * d_nz_xtrct(b,2) - delr_stkstk_norm[2] * d_nz_xtrct(b,1);
        t6dir[1] = delr_stkstk_norm[2] * d_nz_xtrct(b,0) - delr_stkstk_norm[0] * d_nz_xtrct(b,2);
        t6dir[2] = delr_stkstk_norm[0] * d_nz_xtrct(b,1) - delr_stkstk_norm[1] * d_nz_xtrct(b,0);
        deltb[0] -= t6dir[0] * tpair;
        deltb[1] -= t6dir[1] * tpair;
        deltb[2] -= t6dir[2] * tpair;
      }

      // increment torques

      a_torque(a,0) += delta[0];
      a_torque(a,1) += delta[1];
      a_torque(a,2) += delta[2];

      if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
        a_torque(b,0) -= deltb[0];
        a_torque(b,1) -= deltb[1];
        a_torque(b,2) -= deltb[2];
      }
    // end of early rejection criterion
    } // evdwl
    } // f4t5
    } // f4t4
    } // f4f6t1
  }
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::operator()(TagPairOxdna2CoaxstkCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>\
  (TagPairOxdna2CoaxstkCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ia,ev);
}

/* ----------------------------------------------------------------------
   ComputeGPUPair Functor(s)
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_theta1_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  KK_FLOAT &theta1, KK_FLOAT &theta1p, KK_FLOAT &f4f6t1, KK_FLOAT &df4f6t1) const
{
  const KK_FLOAT a1 = d_a_cxst1(atype,btype);
  const KK_FLOAT t10 = d_theta_cxst1_0(atype,btype);
  const KK_FLOAT dt1a = d_dtheta_cxst1_ast(atype,btype);
  const KK_FLOAT b1 = d_b_cxst1(atype,btype);
  const KK_FLOAT dt1c = d_dtheta_cxst1_c(atype,btype);
  const KK_FLOAT aa1 = d_AA_cxst1(atype,btype);
  const KK_FLOAT bb1 = d_BB_cxst1(atype,btype);

  KK_FLOAT cost1 = -Kokkos::fma(a_nx[2], b_nx[2], Kokkos::fma(a_nx[1], b_nx[1], a_nx[0] * b_nx[0]));
  if (cost1 >  1.0) cost1 =  1.0;
  if (cost1 < -1.0) cost1 = -1.0;
  theta1 = acos(cost1);
  theta1p = 2 * MY_PI_KK - theta1;

  f4f6t1 = F4_KK(theta1, a1, t10, dt1a, b1, dt1c) + F6_KK(theta1, aa1, bb1);
  if (f4f6t1 == static_cast<KK_FLOAT>(0.0)) return false;

  // df4f6t1 = (DF4 + DF6) / sin(theta1)
  KK_FLOAT sin1_sq = Kokkos::fma(-cost1, cost1, static_cast<KK_FLOAT>(1.0));
  if (sin1_sq <= 0.0) return false;
  KK_FLOAT sin1 = sqrtf(sin1_sq);
  df4f6t1 = (DF4_KK(theta1, a1, t10, dt1a, b1, dt1c) + DF6_KK(theta1, aa1, bb1)) / sin1;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_theta4_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  KK_FLOAT &theta4, KK_FLOAT &f4t4, KK_FLOAT &df4t4) const
{
  const KK_FLOAT a4 = d_a_cxst4(atype,btype);
  const KK_FLOAT t40 = d_theta_cxst4_0(atype,btype);
  const KK_FLOAT dt4a = d_dtheta_cxst4_ast(atype,btype);
  const KK_FLOAT b4 = d_b_cxst4(atype,btype);
  const KK_FLOAT dt4c = d_dtheta_cxst4_c(atype,btype);

  KK_FLOAT cost4 = Kokkos::fma(a_nz[2], b_nz[2], Kokkos::fma(a_nz[1], b_nz[1], a_nz[0] * b_nz[0]));
  if (cost4 > 1.0) cost4 = 1.0;
  if (cost4 < -1.0) cost4 = -1.0;
  theta4 = acos(cost4);
  f4t4 = F4_KK(theta4, a4, t40, dt4a, b4, dt4c) +
         F4_KK(theta4, a4, MY_PI - t40, dt4a, b4, dt4c);
  if (f4t4 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin4_sq = Kokkos::fma(-cost4, cost4, static_cast<KK_FLOAT>(1.0));
  if (sin4_sq <= 0.0) return false;
  df4t4 = ( DF4_KK(theta4, a4, t40, dt4a, b4, dt4c) +
            DF4_KK(theta4, a4, MY_PI - t40, dt4a, b4, dt4c) ) / sqrtf(sin4_sq);
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_theta5_terms(const int &atype, const int &btype, const KK_FLOAT (&a_nz)[3],
  const KK_FLOAT (&delr_stkstk_norm)[3], KK_FLOAT &theta5, KK_FLOAT &theta5p, KK_FLOAT &f4t5, KK_FLOAT &df4t5, KK_FLOAT &cost5) const
{
  const KK_FLOAT a5 = d_a_cxst5(atype,btype);
  const KK_FLOAT t50 = d_theta_cxst5_0(atype,btype);
  const KK_FLOAT dt5a = d_dtheta_cxst5_ast(atype,btype);
  const KK_FLOAT b5 = d_b_cxst5(atype,btype);
  const KK_FLOAT dt5c = d_dtheta_cxst5_c(atype,btype);

  cost5 = Kokkos::fma(a_nz[2], delr_stkstk_norm[2], Kokkos::fma(a_nz[1], delr_stkstk_norm[1], a_nz[0] * delr_stkstk_norm[0]));
  if (cost5 > 1.0) cost5 = 1.0;
  if (cost5 < -1.0) cost5 = -1.0;
  theta5 = acos(cost5);
  theta5p = MY_PI_KK - theta5;
  f4t5 = F4_KK(theta5, a5, t50, dt5a, b5, dt5c) +
         F4_KK(theta5p, a5, t50, dt5a, b5, dt5c);
  if (f4t5 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin5_sq = Kokkos::fma(-cost5, cost5, static_cast<KK_FLOAT>(1.0));
  if (sin5_sq <= 0.0) return false;
  df4t5 = ( DF4_KK(theta5, a5, t50, dt5a, b5, dt5c) -
            DF4_KK(theta5p, a5, t50, dt5a, b5, dt5c) ) / sqrtf(sin5_sq);
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_theta6_terms(const int &atype, const int &btype, const KK_FLOAT (&b_nz)[3],
  const KK_FLOAT (&delr_stkstk_norm)[3], KK_FLOAT &theta6, KK_FLOAT &theta6p, KK_FLOAT &f4t6, KK_FLOAT &df4t6, KK_FLOAT &cost6) const
{
  const KK_FLOAT a6 = d_a_cxst6(atype,btype);
  const KK_FLOAT t60 = d_theta_cxst6_0(atype,btype);
  const KK_FLOAT dt6a = d_dtheta_cxst6_ast(atype,btype);
  const KK_FLOAT b6 = d_b_cxst6(atype,btype);
  const KK_FLOAT dt6c = d_dtheta_cxst6_c(atype,btype);

  cost6 = Kokkos::fma(b_nz[2], delr_stkstk_norm[2], Kokkos::fma(b_nz[1], delr_stkstk_norm[1], b_nz[0] * delr_stkstk_norm[0]));
  if (cost6 > 1.0) cost6 = 1.0;
  if (cost6 < -1.0) cost6 = -1.0;
  theta6 = acos(cost6);
  theta6p = MY_PI_KK - theta6;
  f4t6 = F4_KK(theta6, a6, t60, dt6a, b6, dt6c) +
         F4_KK(theta6p, a6, t60, dt6a, b6, dt6c);
  if (f4t6 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin6_sq = Kokkos::fma(-cost6, cost6, static_cast<KK_FLOAT>(1.0));
  if (sin6_sq <= 0.0) return false;
  df4t6 = ( DF4_KK(theta6, a6, t60, dt6a, b6, dt6c) -
            DF4_KK(theta6p, a6, t60, dt6a, b6, dt6c) ) / sqrtf(sin6_sq);
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_cosphi3_terms(const int &a, const int &b,
  const KK_FLOAT (&ra_cbk)[3], const KK_FLOAT (&rb_cbk)[3], const KK_FLOAT (&a_nx)[3],
  const KK_FLOAT (&delr_stkstk_norm)[3], KK_FLOAT &cosphi3) const
{
  KK_FLOAT delr_bkbk[3], delr_bkbk_norm[3];
  delr_bkbk[0] = x(a,0) + ra_cbk[0] - x(b,0) - rb_cbk[0];
  delr_bkbk[1] = x(a,1) + ra_cbk[1] - x(b,1) - rb_cbk[1];
  delr_bkbk[2] = x(a,2) + ra_cbk[2] - x(b,2) - rb_cbk[2];

  KK_FLOAT rsq_bkbk = Kokkos::fma(delr_bkbk[0], delr_bkbk[0],
                    Kokkos::fma(delr_bkbk[1], delr_bkbk[1], delr_bkbk[2] * delr_bkbk[2]));
  KK_FLOAT rinv_bkbk = 1.0 / sqrtf(rsq_bkbk);
  delr_bkbk_norm[0] = delr_bkbk[0] * rinv_bkbk;
  delr_bkbk_norm[1] = delr_bkbk[1] * rinv_bkbk;
  delr_bkbk_norm[2] = delr_bkbk[2] * rinv_bkbk;

  KK_FLOAT v1tmp0 = delr_bkbk_norm[1] * a_nx[2] - delr_bkbk_norm[2] * a_nx[1];
  KK_FLOAT v1tmp1 = delr_bkbk_norm[2] * a_nx[0] - delr_bkbk_norm[0] * a_nx[2];
  KK_FLOAT v1tmp2 = delr_bkbk_norm[0] * a_nx[1] - delr_bkbk_norm[1] * a_nx[0];
  cosphi3 = Kokkos::fma(v1tmp2, delr_stkstk_norm[2], Kokkos::fma(v1tmp1, delr_stkstk_norm[1], v1tmp0 * delr_stkstk_norm[0]));
  if (cosphi3 > 1.0) cosphi3 = 1.0;
  if (cosphi3 < -1.0) cosphi3 = -1.0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_radial_terms(const int &atype, const int &btype, const KK_FLOAT &r_st,
  const KK_FLOAT &prime_cxst_ab,
  KK_FLOAT &f2, KK_FLOAT &df2) const
{
  f2 = F2_KK(r_st, prime_cxst_ab, d_cut_cxst_0(atype,btype), d_cut_cxst_lc(atype,btype),
             d_cut_cxst_hc(atype,btype), d_cut_cxst_lo(atype,btype), d_cut_cxst_hi(atype,btype),
             d_b_cxst_lo(atype,btype), d_b_cxst_hi(atype,btype), d_cut_cxst_c(atype,btype));
  if (f2 == static_cast<KK_FLOAT>(0.0)) return false;
  df2 = DF2_KK(r_st, prime_cxst_ab, d_cut_cxst_0(atype,btype), d_cut_cxst_lc(atype,btype),
               d_cut_cxst_hc(atype,btype), d_cut_cxst_lo(atype,btype), d_cut_cxst_hi(atype,btype),
               d_b_cxst_lo(atype,btype), d_b_cxst_hi(atype,btype));
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_force_contrib(const KK_FLOAT &df2, const KK_FLOAT &f2, const KK_FLOAT &f4f6t1,
  const KK_FLOAT &f4t4, const KK_FLOAT &f4t5, const KK_FLOAT &f4t6, const KK_FLOAT &df4t5, const KK_FLOAT &df4t6, const KK_FLOAT &rinv_st,
  const KK_FLOAT &factor_lj, const KK_FLOAT &cost5, const KK_FLOAT &cost6,
  const KK_FLOAT &theta5, const KK_FLOAT &theta5p, const KK_FLOAT &theta6, const KK_FLOAT &theta6p,
  const KK_FLOAT (&delr_stkstk)[3], const KK_FLOAT (&delr_stkstk_norm)[3], const KK_FLOAT (&a_nz)[3],
  const KK_FLOAT (&b_nz)[3], const KK_FLOAT (&ra_cstk)[3], const KK_FLOAT (&rb_cstk)[3],
  KK_ACC_FLOAT (&delf)[3], KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const
{
  delf[0] = delf[1] = delf[2] = 0.0;

  // radial contribution
  KK_FLOAT product = f4f6t1;
  product = product * f4t4;
  product = product * f4t5;
  product = product * f4t6;
  KK_FLOAT finc = -df2 * product * rinv_st * factor_lj;
  delf[0] = Kokkos::fma(delr_stkstk[0], finc, delf[0]);
  delf[1] = Kokkos::fma(delr_stkstk[1], finc, delf[1]);
  delf[2] = Kokkos::fma(delr_stkstk[2], finc, delf[2]);

  // theta5
  if (theta5 != static_cast<KK_FLOAT>(0.0) && theta5p != static_cast<KK_FLOAT>(0.0)) {
    KK_FLOAT prod5 = f4f6t1 * f4t4 * df4t5 * f4t6;
    KK_FLOAT finc5 = -f2 * prod5 * rinv_st * factor_lj;
    delf[0] = Kokkos::fma(Kokkos::fma(delr_stkstk_norm[0], cost5, -a_nz[0]), finc5, delf[0]);
    delf[1] = Kokkos::fma(Kokkos::fma(delr_stkstk_norm[1], cost5, -a_nz[1]), finc5, delf[1]);
    delf[2] = Kokkos::fma(Kokkos::fma(delr_stkstk_norm[2], cost5, -a_nz[2]), finc5, delf[2]);
  }

  // theta6
  if (theta6 != static_cast<KK_FLOAT>(0.0) && theta6p != static_cast<KK_FLOAT>(0.0)) {
    KK_FLOAT prod6 = f4f6t1 * f4t4 * f4t5 * df4t6;
    KK_FLOAT finc6 = -f2 * prod6 * rinv_st * factor_lj;
    delf[0] = Kokkos::fma(Kokkos::fma(delr_stkstk_norm[0], cost6, -b_nz[0]), finc6, delf[0]);
    delf[1] = Kokkos::fma(Kokkos::fma(delr_stkstk_norm[1], cost6, -b_nz[1]), finc6, delf[1]);
    delf[2] = Kokkos::fma(Kokkos::fma(delr_stkstk_norm[2], cost6, -b_nz[2]), finc6, delf[2]);
  }

  // r x f torques
  delta[0] = Kokkos::fma(ra_cstk[1], delf[2], - ra_cstk[2] * delf[1]);
  delta[1] = Kokkos::fma(ra_cstk[2], delf[0], - ra_cstk[0] * delf[2]);
  delta[2] = Kokkos::fma(ra_cstk[0], delf[1], - ra_cstk[1] * delf[0]);
  deltb[0] = Kokkos::fma(rb_cstk[1], delf[2], - rb_cstk[2] * delf[1]);
  deltb[1] = Kokkos::fma(rb_cstk[2], delf[0], - rb_cstk[0] * delf[2]);
  deltb[2] = Kokkos::fma(rb_cstk[0], delf[1], - rb_cstk[1] * delf[0]);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::coaxstk_torque_contrib(const KK_FLOAT &f2, const KK_FLOAT &df4f6t1, const KK_FLOAT &f4f6t1,
  const KK_FLOAT &f4t4, const KK_FLOAT &f4t5, const KK_FLOAT &f4t6, const KK_FLOAT &df4t4, const KK_FLOAT &df4t5, const KK_FLOAT &df4t6,
  const KK_FLOAT &factor_lj, const KK_FLOAT &theta1, const KK_FLOAT &theta1p, const KK_FLOAT &theta4, const KK_FLOAT &theta5,
  const KK_FLOAT &theta5p, const KK_FLOAT &theta6, const KK_FLOAT &theta6p,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3], const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  const KK_FLOAT (&delr_stkstk_norm)[3], KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const
{
  // theta1
  if (theta1 != static_cast<KK_FLOAT>(0.0) && theta1p != static_cast<KK_FLOAT>(0.0)) {
    KK_FLOAT tpair = -f2 * df4f6t1 * f4t4 * f4t5 * f4t6 * factor_lj;
    KK_FLOAT t1x = Kokkos::fma(a_nx[1], b_nx[2], - a_nx[2] * b_nx[1]);
    KK_FLOAT t1y = Kokkos::fma(a_nx[2], b_nx[0], - a_nx[0] * b_nx[2]);
    KK_FLOAT t1z = Kokkos::fma(a_nx[0], b_nx[1], - a_nx[1] * b_nx[0]);
    delta[0] = Kokkos::fma(t1x, tpair, delta[0]);
    delta[1] = Kokkos::fma(t1y, tpair, delta[1]);
    delta[2] = Kokkos::fma(t1z, tpair, delta[2]);
    deltb[0] = Kokkos::fma(t1x, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(t1y, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(t1z, tpair, deltb[2]);
  }

  // theta4
  if (theta4 != static_cast<KK_FLOAT>(0.0)) {
    KK_FLOAT tpair = -f2 * f4f6t1 * df4t4 * f4t5 * f4t6 * factor_lj;
    KK_FLOAT t4x = Kokkos::fma(b_nz[1], a_nz[2], - b_nz[2] * a_nz[1]);
    KK_FLOAT t4y = Kokkos::fma(b_nz[2], a_nz[0], - b_nz[0] * a_nz[2]);
    KK_FLOAT t4z = Kokkos::fma(b_nz[0], a_nz[1], - b_nz[1] * a_nz[0]);
    delta[0] = Kokkos::fma(t4x, tpair, delta[0]);
    delta[1] = Kokkos::fma(t4y, tpair, delta[1]);
    delta[2] = Kokkos::fma(t4z, tpair, delta[2]);
    deltb[0] = Kokkos::fma(t4x, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(t4y, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(t4z, tpair, deltb[2]);
  }

  // theta5
  if (theta5 != static_cast<KK_FLOAT>(0.0) && theta5p != static_cast<KK_FLOAT>(0.0)) {
    KK_FLOAT tpair = -f2 * f4f6t1 * f4t4 * df4t5 * f4t6 * factor_lj;
    KK_FLOAT t5x = Kokkos::fma(delr_stkstk_norm[1], a_nz[2], - delr_stkstk_norm[2] * a_nz[1]);
    KK_FLOAT t5y = Kokkos::fma(delr_stkstk_norm[2], a_nz[0], - delr_stkstk_norm[0] * a_nz[2]);
    KK_FLOAT t5z = Kokkos::fma(delr_stkstk_norm[0], a_nz[1], - delr_stkstk_norm[1] * a_nz[0]);
    delta[0] = Kokkos::fma(t5x, tpair, delta[0]);
    delta[1] = Kokkos::fma(t5y, tpair, delta[1]);
    delta[2] = Kokkos::fma(t5z, tpair, delta[2]);
  }

  // theta6
  if (theta6 != static_cast<KK_FLOAT>(0.0) && theta6p != static_cast<KK_FLOAT>(0.0)) {
    KK_FLOAT tpair = -f2 * f4f6t1 * f4t4 * f4t5 * df4t6 * factor_lj;
    KK_FLOAT t6x = Kokkos::fma(delr_stkstk_norm[1], b_nz[2], - delr_stkstk_norm[2] * b_nz[1]);
    KK_FLOAT t6y = Kokkos::fma(delr_stkstk_norm[2], b_nz[0], - delr_stkstk_norm[0] * b_nz[2]);
    KK_FLOAT t6z = Kokkos::fma(delr_stkstk_norm[0], b_nz[1], - delr_stkstk_norm[1] * b_nz[0]);
    deltb[0] = Kokkos::fma(-t6x, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(-t6y, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(-t6z, tpair, deltb[2]);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::operator()(TagPairOxdna2CoaxstkComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ipair, EV_FLOAT &ev) const
{
  // f and torque array are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  // Direct packed pair lookup: high 32 bits = a, low 32 bits = b.
  const uint64_t pair = d_pairs_screened(ipair);
  // "pair >> 32" shifts the pair to the right by 32 bits, so the upper 32 bits
  // becomes the lower 32 bits to recover the atom-a index.
  const int a = static_cast<int>(pair >> 32);
  // "pair & 0xffffffffu" keeps only the lower 32 bits to recover the atom-b index.
  int b = static_cast<int>(pair & 0xffffffffu);
  const KK_FLOAT factor_lj = special_lj[sbmask(b)];
  if (factor_lj == static_cast<KK_FLOAT>(0.0)) return;
  b &= NEIGHMASK;

  // a has to be terminal nucleotide
  if(id3p[a]!=-1 && id5p[a]!=-1) return;
  // b has to be terminal nucleotide
  if(id3p[b]!=-1 && id5p[b]!=-1) return;

  const int atype = type(a);
  const int btype = type(b);

  // vectors COM-backbone site, COM-stacking site in lab frame
  KK_FLOAT ra_cbk[3], rb_cbk[3], ra_cstk[3], rb_cstk[3];

  KK_ACC_FLOAT delf[3],delta[3],deltb[3];    // force, torque increment
  KK_ACC_FLOAT evdwl;                        // energy
  KK_FLOAT delr_stkstk[3],delr_stkstk_norm[3],rsq_stkstk,r_stkstk,rinv_stkstk;
  // NOTE: delr_bkbk[]3, etc is scoped out into coaxstk_cosphi3_terms to reduce register pressure
  KK_FLOAT theta1,theta1p;
  KK_FLOAT theta4;
  KK_FLOAT theta5,theta5p,cost5;
  KK_FLOAT theta6,theta6p,cost6;
  KK_FLOAT cosphi3;
  KK_FLOAT prime_cxst_ab;

  KK_FLOAT f2,f4f6t1,f4t4,f4t5,f4t6;
  KK_FLOAT df2,df4f6t1,df4t4,df4t5,df4t6;

  // single loads for local axes to reduce repeated global reads
  const KK_FLOAT a_nx_loc[3] = { d_nx_xtrct(a,0), d_nx_xtrct(a,1), d_nx_xtrct(a,2) };
  const KK_FLOAT b_nx_loc[3] = { d_nx_xtrct(b,0), d_nx_xtrct(b,1), d_nx_xtrct(b,2) };
  const KK_FLOAT a_ny_loc[3] = { d_ny_xtrct(a,0), d_ny_xtrct(a,1), d_ny_xtrct(a,2) };
  const KK_FLOAT b_ny_loc[3] = { d_ny_xtrct(b,0), d_ny_xtrct(b,1), d_ny_xtrct(b,2) };

  // vector COM-backbone site [a/b], COM-stacking site [a/b]
  if constexpr (OXDNAFLAG==OXDNA2) {
    constexpr KK_FLOAT dx_cbk_oxdna2 = -0.34;
    constexpr KK_FLOAT dy_cbk_oxdna2 = +0.3408;
    ra_cbk[0] = dx_cbk_oxdna2*a_nx_loc[0] + dy_cbk_oxdna2*a_ny_loc[0];
    ra_cbk[1] = dx_cbk_oxdna2*a_nx_loc[1] + dy_cbk_oxdna2*a_ny_loc[1];
    ra_cbk[2] = dx_cbk_oxdna2*a_nx_loc[2] + dy_cbk_oxdna2*a_ny_loc[2];
    rb_cbk[0] = dx_cbk_oxdna2*b_nx_loc[0] + dy_cbk_oxdna2*b_ny_loc[0];
    rb_cbk[1] = dx_cbk_oxdna2*b_nx_loc[1] + dy_cbk_oxdna2*b_ny_loc[1];
    rb_cbk[2] = dx_cbk_oxdna2*b_nx_loc[2] + dy_cbk_oxdna2*b_ny_loc[2];
    constexpr KK_FLOAT dx_cstk_oxdna1 = +0.34;  // oxDNA2 uses same stacking site as oxDNA1
    ra_cstk[0] = dx_cstk_oxdna1*a_nx_loc[0];
    ra_cstk[1] = dx_cstk_oxdna1*a_nx_loc[1];
    ra_cstk[2] = dx_cstk_oxdna1*a_nx_loc[2];
    rb_cstk[0] = dx_cstk_oxdna1*b_nx_loc[0];
    rb_cstk[1] = dx_cstk_oxdna1*b_nx_loc[1];
    rb_cstk[2] = dx_cstk_oxdna1*b_nx_loc[2];
  } else if constexpr (OXDNAFLAG==OXDNA3) {
    // oxDNA3 uses same backbone site as oxDNA2...
    constexpr KK_FLOAT dx_cbk_oxdna2 = -0.34;
    constexpr KK_FLOAT dy_cbk_oxdna2 = +0.3408;
    ra_cbk[0] = dx_cbk_oxdna2*a_nx_loc[0] + dy_cbk_oxdna2*a_ny_loc[0];
    ra_cbk[1] = dx_cbk_oxdna2*a_nx_loc[1] + dy_cbk_oxdna2*a_ny_loc[1];
    ra_cbk[2] = dx_cbk_oxdna2*a_nx_loc[2] + dy_cbk_oxdna2*a_ny_loc[2];
    rb_cbk[0] = dx_cbk_oxdna2*b_nx_loc[0] + dy_cbk_oxdna2*b_ny_loc[0];
    rb_cbk[1] = dx_cbk_oxdna2*b_nx_loc[1] + dy_cbk_oxdna2*b_ny_loc[1];
    rb_cbk[2] = dx_cbk_oxdna2*b_nx_loc[2] + dy_cbk_oxdna2*b_ny_loc[2];
    // ...But the stacking site is different for oxDNA3.
    constexpr KK_FLOAT dx_cstk_oxdna3 = +0.37;
    ra_cstk[0] = dx_cstk_oxdna3*a_nx_loc[0];
    ra_cstk[1] = dx_cstk_oxdna3*a_nx_loc[1];
    ra_cstk[2] = dx_cstk_oxdna3*a_nx_loc[2];
    rb_cstk[0] = dx_cstk_oxdna3*b_nx_loc[0];
    rb_cstk[1] = dx_cstk_oxdna3*b_nx_loc[1];
    rb_cstk[2] = dx_cstk_oxdna3*b_nx_loc[2];
  }

  // vector stacking site b to a
  // stkstk is needed for theta5/6 and radial terms, so we do not scope out....
  delr_stkstk[0] = x(a,0) + ra_cstk[0] - x(b,0) - rb_cstk[0];
  delr_stkstk[1] = x(a,1) + ra_cstk[1] - x(b,1) - rb_cstk[1];
  delr_stkstk[2] = x(a,2) + ra_cstk[2] - x(b,2) - rb_cstk[2];
  rsq_stkstk = Kokkos::fma(delr_stkstk[0], delr_stkstk[0], Kokkos::fma(delr_stkstk[1], delr_stkstk[1], delr_stkstk[2]*delr_stkstk[2]));
  r_stkstk = sqrtf(rsq_stkstk);
  rinv_stkstk = 1.0 / r_stkstk;
  delr_stkstk_norm[0] = delr_stkstk[0] * rinv_stkstk;
  delr_stkstk_norm[1] = delr_stkstk[1] * rinv_stkstk;
  delr_stkstk_norm[2] = delr_stkstk[2] * rinv_stkstk;
  // .... but bkbk (vector backbone site b to a) is only needed for cosphi3, so we scope out to reduce register pressure

  const KK_FLOAT a_nz_loc[3] = { d_nz_xtrct(a,0), d_nz_xtrct(a,1), d_nz_xtrct(a,2) };
  const KK_FLOAT b_nz_loc[3] = { d_nz_xtrct(b,0), d_nz_xtrct(b,1), d_nz_xtrct(b,2) };

  // beginning of modulation factors
  if (!coaxstk_theta1_terms(atype,btype,a_nx_loc,b_nx_loc,theta1,theta1p,f4f6t1,df4f6t1)) return;
  if (!coaxstk_theta4_terms(atype,btype,a_nz_loc,b_nz_loc,theta4,f4t4,df4t4)) return;
  if (!coaxstk_theta5_terms(atype,btype,a_nz_loc,delr_stkstk_norm,theta5,theta5p,f4t5,df4t5,cost5)) return;
  if (!coaxstk_theta6_terms(atype,btype,b_nz_loc,delr_stkstk_norm,theta6,theta6p,f4t6,df4t6,cost6)) return;
  // cosphi3 is just scoped out for sake of register pressure, but does not feature any early exit criteria
  coaxstk_cosphi3_terms(a, b, ra_cbk, rb_cbk, a_nx_loc, delr_stkstk_norm, cosphi3);

  // Direction-dependent coaxial stacking strength.
  if (id5p(a) == -1 && id3p(b) == -1) {
    prime_cxst_ab = d_k_cxst(atype,btype);
  } else if (id3p(a) == -1 && id5p(b) == -1) {
    prime_cxst_ab = d_k_cxst(btype,atype);
  } else {
    prime_cxst_ab = 0.5 * (d_k_cxst(atype,btype) + d_k_cxst(btype,atype));
  }

  if (!coaxstk_radial_terms(atype,btype,r_stkstk,prime_cxst_ab,f2,df2)) return;

  evdwl = f2 * f4f6t1 * f4t4 * f4t5 * f4t6 * factor_lj;

  // force, torque, and viral contributions for forces between h-bonding sites

  delf[0] = delf[1] = delf[2] = 0.0;
  delta[0] = delta[1] = delta[2] = 0.0;
  deltb[0] = deltb[1] = deltb[2] = 0.0;

  coaxstk_force_contrib(df2,f2,f4f6t1,f4t4,f4t5,f4t6,df4t5,df4t6,
      rinv_stkstk,factor_lj,cost5,cost6,theta5,theta5p,theta6,theta6p,
      delr_stkstk,delr_stkstk_norm,a_nz_loc,b_nz_loc,ra_cstk,rb_cstk,delf,delta,deltb);

  // increment forces and torques
  a_f(a,0) += delf[0];
  a_f(a,1) += delf[1];
  a_f(a,2) += delf[2];
  a_torque(a,0) += delta[0];
  a_torque(a,1) += delta[1];
  a_torque(a,2) += delta[2];

  if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
    a_f(b,0) -= delf[0];
    a_f(b,1) -= delf[1];
    a_f(b,2) -= delf[2];
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }

  // increment energy and virial
  // NOTE: The virial is calculated on the 'molecular' basis.
  // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

  if (EVFLAG) {
    ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;

    if (vflag_either || eflag_atom) {
      this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,evdwl,\
      delf[0],delf[1],delf[2],x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2));
    }
  }

  // pure torques not expressible as r x f

  delta[0] = delta[1] = delta[2] = 0.0;
  deltb[0] = deltb[1] = deltb[2] = 0.0;

  // compute pure torques via staged torque helper
  coaxstk_torque_contrib(f2,df4f6t1,f4f6t1,f4t4,f4t5,f4t6,df4t4,df4t5,df4t6,factor_lj,
    theta1,theta1p,theta4,theta5,theta5p,theta6,theta6p,
    a_nx_loc,b_nx_loc,a_nz_loc,b_nz_loc,delr_stkstk_norm,delta,deltb);

  // increment torques
  a_torque(a,0) += delta[0];
  a_torque(a,1) += delta[1];
  a_torque(a,2) += delta[2];

  if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }
// end of early rejection criterion
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::operator()(TagPairOxdna2CoaxstkComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ipair) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>\
  (TagPairOxdna2CoaxstkComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ipair,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2CoaxstkKokkos<DeviceType>::allocate()
{
  PairOxdna2Coaxstk::allocate();

  int n = atom->ntypes;

  memoryKK->create_kokkos(k_k_cxst,n+1,n+1,"PairOxdna2Coaxstk:k_cxst");
  memoryKK->create_kokkos(k_cut_cxst_0,n+1,n+1,"PairOxdna2Coaxstk:cut_cxst_0");
  memoryKK->create_kokkos(k_cut_cxst_c,n+1,n+1,"PairOxdna2Coaxstk:cut_cxst_c");
  memoryKK->create_kokkos(k_cut_cxst_lo,n+1,n+1,"PairOxdna2Coaxstk:cut_cxst_lo");
  memoryKK->create_kokkos(k_cut_cxst_hi,n+1,n+1,"PairOxdna2Coaxstk:cut_cxst_hi");
  memoryKK->create_kokkos(k_cut_cxst_lc,n+1,n+1,"PairOxdna2Coaxstk:cut_cxst_lc");
  memoryKK->create_kokkos(k_cut_cxst_hc,n+1,n+1,"PairOxdna2Coaxstk:cut_cxst_hc");
  memoryKK->create_kokkos(k_b_cxst_lo,n+1,n+1,"PairOxdna2Coaxstk:b_cxst_lo");
  memoryKK->create_kokkos(k_b_cxst_hi,n+1,n+1,"PairOxdna2Coaxstk:b_cxst_hi");
  memoryKK->create_kokkos(k_cutsq_cxst_hc,n+1,n+1,"PairOxdna2Coaxstk:cutsq_cxst_hc");

  memoryKK->create_kokkos(k_a_cxst1,n+1,n+1,"PairOxdna2Coaxstk:a_cxst1");
  memoryKK->create_kokkos(k_theta_cxst1_0,n+1,n+1,"PairOxdna2Coaxstk:theta_cxst1_0");
  memoryKK->create_kokkos(k_dtheta_cxst1_ast,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst1_ast");
  memoryKK->create_kokkos(k_b_cxst1,n+1,n+1,"PairOxdna2Coaxstk:b_cxst1");
  memoryKK->create_kokkos(k_dtheta_cxst1_c,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst1_c");

  memoryKK->create_kokkos(k_a_cxst4,n+1,n+1,"PairOxdna2Coaxstk:a_cxst4");
  memoryKK->create_kokkos(k_theta_cxst4_0,n+1,n+1,"PairOxdna2Coaxstk:theta_cxst4_0");
  memoryKK->create_kokkos(k_dtheta_cxst4_ast,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst4_ast");
  memoryKK->create_kokkos(k_b_cxst4,n+1,n+1,"PairOxdna2Coaxstk:b_cxst4");
  memoryKK->create_kokkos(k_dtheta_cxst4_c,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst4_c");

  memoryKK->create_kokkos(k_a_cxst5,n+1,n+1,"PairOxdna2Coaxstk:a_cxst5");
  memoryKK->create_kokkos(k_theta_cxst5_0,n+1,n+1,"PairOxdna2Coaxstk:theta_cxst5_0");
  memoryKK->create_kokkos(k_dtheta_cxst5_ast,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst5_ast");
  memoryKK->create_kokkos(k_b_cxst5,n+1,n+1,"PairOxdna2Coaxstk:b_cxst5");
  memoryKK->create_kokkos(k_dtheta_cxst5_c,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst5_c");

  memoryKK->create_kokkos(k_a_cxst6,n+1,n+1,"PairOxdna2Coaxstk:a_cxst6");
  memoryKK->create_kokkos(k_theta_cxst6_0,n+1,n+1,"PairOxdna2Coaxstk:theta_cxst6_0");
  memoryKK->create_kokkos(k_dtheta_cxst6_ast,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst6_ast");
  memoryKK->create_kokkos(k_b_cxst6,n+1,n+1,"PairOxdna2Coaxstk:b_cxst6");
  memoryKK->create_kokkos(k_dtheta_cxst6_c,n+1,n+1,"PairOxdna2Coaxstk:dtheta_cxst6_c");

  memoryKK->create_kokkos(k_AA_cxst1,n+1,n+1,"PairOxdna2Coaxstk:AA_cxst1");
  memoryKK->create_kokkos(k_BB_cxst1,n+1,n+1,"PairOxdna2Coaxstk:BB_cxst1");

  d_k_cxst = k_k_cxst.template view<DeviceType>();
  d_cut_cxst_0 = k_cut_cxst_0.template view<DeviceType>();
  d_cut_cxst_c = k_cut_cxst_c.template view<DeviceType>();
  d_cut_cxst_lo = k_cut_cxst_lo.template view<DeviceType>();
  d_cut_cxst_hi = k_cut_cxst_hi.template view<DeviceType>();
  d_cut_cxst_lc = k_cut_cxst_lc.template view<DeviceType>();
  d_cut_cxst_hc = k_cut_cxst_hc.template view<DeviceType>();
  d_b_cxst_lo = k_b_cxst_lo.template view<DeviceType>();
  d_b_cxst_hi = k_b_cxst_hi.template view<DeviceType>();
  d_cutsq_cxst_hc = k_cutsq_cxst_hc.template view<DeviceType>();

  d_a_cxst1 = k_a_cxst1.template view<DeviceType>();
  d_theta_cxst1_0 = k_theta_cxst1_0.template view<DeviceType>();
  d_dtheta_cxst1_ast = k_dtheta_cxst1_ast.template view<DeviceType>();
  d_b_cxst1 = k_b_cxst1.template view<DeviceType>();
  d_dtheta_cxst1_c = k_dtheta_cxst1_c.template view<DeviceType>();

  d_a_cxst4 = k_a_cxst4.template view<DeviceType>();
  d_theta_cxst4_0 = k_theta_cxst4_0.template view<DeviceType>();
  d_dtheta_cxst4_ast = k_dtheta_cxst4_ast.template view<DeviceType>();
  d_b_cxst4 = k_b_cxst4.template view<DeviceType>();
  d_dtheta_cxst4_c = k_dtheta_cxst4_c.template view<DeviceType>();

  d_a_cxst5 = k_a_cxst5.template view<DeviceType>();
  d_theta_cxst5_0 = k_theta_cxst5_0.template view<DeviceType>();
  d_dtheta_cxst5_ast = k_dtheta_cxst5_ast.template view<DeviceType>();
  d_b_cxst5 = k_b_cxst5.template view<DeviceType>();
  d_dtheta_cxst5_c = k_dtheta_cxst5_c.template view<DeviceType>();

  d_a_cxst6 = k_a_cxst6.template view<DeviceType>();
  d_theta_cxst6_0 = k_theta_cxst6_0.template view<DeviceType>();
  d_dtheta_cxst6_ast = k_dtheta_cxst6_ast.template view<DeviceType>();
  d_b_cxst6 = k_b_cxst6.template view<DeviceType>();
  d_dtheta_cxst6_c = k_dtheta_cxst6_c.template view<DeviceType>();

  d_AA_cxst1 = k_AA_cxst1.template view<DeviceType>();
  d_BB_cxst1 = k_BB_cxst1.template view<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2CoaxstkKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2CoaxstkKokkos<DeviceType>::init_style()
{
  neighbor->add_request(this);
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  fix_oxdna_lrfKK = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);

  fix_oxdna_npairKK = nullptr;
  auto npair_fixes = modify->get_fix_by_style("^OXDNA/NPAIR/kk");
  if (npair_fixes.size() == 0) {
    fix_oxdna_npairKK = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(modify->add_fix("npair_kk all OXDNA/NPAIR/kk"));
  } else {
    fix_oxdna_npairKK = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(npair_fixes[0]);
  }
  if (!fix_oxdna_npairKK) error->all(FLERR, "Fix OXDNA/NPAIR/kk lookup failed");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxdna2CoaxstkKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxdna2Coaxstk::init_one(i,j);

  // Assign directionally: [i][j] gets [i][j], [j][i] gets [j][i]
  k_k_cxst.view_host()(i,j) = k_cxst[i][j]; k_k_cxst.view_host()(j,i) = k_cxst[j][i];
  k_cut_cxst_0.view_host()(i,j) = cut_cxst_0[i][j]; k_cut_cxst_0.view_host()(j,i) = cut_cxst_0[j][i];
  k_cut_cxst_c.view_host()(i,j) = cut_cxst_c[i][j]; k_cut_cxst_c.view_host()(j,i) = cut_cxst_c[j][i];
  k_cut_cxst_lo.view_host()(i,j) = cut_cxst_lo[i][j]; k_cut_cxst_lo.view_host()(j,i) = cut_cxst_lo[j][i];
  k_cut_cxst_hi.view_host()(i,j) = cut_cxst_hi[i][j]; k_cut_cxst_hi.view_host()(j,i) = cut_cxst_hi[j][i];
  k_cut_cxst_lc.view_host()(i,j) = cut_cxst_lc[i][j]; k_cut_cxst_lc.view_host()(j,i) = cut_cxst_lc[j][i];
  k_cut_cxst_hc.view_host()(i,j) = cut_cxst_hc[i][j]; k_cut_cxst_hc.view_host()(j,i) = cut_cxst_hc[j][i];
  k_b_cxst_lo.view_host()(i,j) = b_cxst_lo[i][j]; k_b_cxst_lo.view_host()(j,i) = b_cxst_lo[j][i];
  k_b_cxst_hi.view_host()(i,j) = b_cxst_hi[i][j]; k_b_cxst_hi.view_host()(j,i) = b_cxst_hi[j][i];
  k_cutsq_cxst_hc.view_host()(i,j) = cutsq_cxst_hc[i][j]; k_cutsq_cxst_hc.view_host()(j,i) = cutsq_cxst_hc[j][i];

  k_a_cxst1.view_host()(i,j) = a_cxst1[i][j]; k_a_cxst1.view_host()(j,i) = a_cxst1[j][i];
  k_theta_cxst1_0.view_host()(i,j) = theta_cxst1_0[i][j]; k_theta_cxst1_0.view_host()(j,i) = theta_cxst1_0[j][i];
  k_dtheta_cxst1_ast.view_host()(i,j) = dtheta_cxst1_ast[i][j]; k_dtheta_cxst1_ast.view_host()(j,i) = dtheta_cxst1_ast[j][i];
  k_b_cxst1.view_host()(i,j) = b_cxst1[i][j]; k_b_cxst1.view_host()(j,i) = b_cxst1[j][i];
  k_dtheta_cxst1_c.view_host()(i,j) = dtheta_cxst1_c[i][j]; k_dtheta_cxst1_c.view_host()(j,i) = dtheta_cxst1_c[j][i];

  k_a_cxst4.view_host()(i,j) = a_cxst4[i][j]; k_a_cxst4.view_host()(j,i) = a_cxst4[j][i];
  k_theta_cxst4_0.view_host()(i,j) = theta_cxst4_0[i][j]; k_theta_cxst4_0.view_host()(j,i) = theta_cxst4_0[j][i];
  k_dtheta_cxst4_ast.view_host()(i,j) = dtheta_cxst4_ast[i][j]; k_dtheta_cxst4_ast.view_host()(j,i) = dtheta_cxst4_ast[j][i];
  k_b_cxst4.view_host()(i,j) = b_cxst4[i][j]; k_b_cxst4.view_host()(j,i) = b_cxst4[j][i];
  k_dtheta_cxst4_c.view_host()(i,j) = dtheta_cxst4_c[i][j]; k_dtheta_cxst4_c.view_host()(j,i) = dtheta_cxst4_c[j][i];

  k_a_cxst5.view_host()(i,j) = a_cxst5[i][j]; k_a_cxst5.view_host()(j,i) = a_cxst5[j][i];
  k_theta_cxst5_0.view_host()(i,j) = theta_cxst5_0[i][j]; k_theta_cxst5_0.view_host()(j,i) = theta_cxst5_0[j][i];
  k_dtheta_cxst5_ast.view_host()(i,j) = dtheta_cxst5_ast[i][j]; k_dtheta_cxst5_ast.view_host()(j,i) = dtheta_cxst5_ast[j][i];
  k_b_cxst5.view_host()(i,j) = b_cxst5[i][j]; k_b_cxst5.view_host()(j,i) = b_cxst5[j][i];
  k_dtheta_cxst5_c.view_host()(i,j) = dtheta_cxst5_c[i][j]; k_dtheta_cxst5_c.view_host()(j,i) = dtheta_cxst5_c[j][i];

  k_a_cxst6.view_host()(i,j) = a_cxst6[i][j]; k_a_cxst6.view_host()(j,i) = a_cxst6[j][i];
  k_theta_cxst6_0.view_host()(i,j) = theta_cxst6_0[i][j]; k_theta_cxst6_0.view_host()(j,i) = theta_cxst6_0[j][i];
  k_dtheta_cxst6_ast.view_host()(i,j) = dtheta_cxst6_ast[i][j]; k_dtheta_cxst6_ast.view_host()(j,i) = dtheta_cxst6_ast[j][i];
  k_b_cxst6.view_host()(i,j) = b_cxst6[i][j]; k_b_cxst6.view_host()(j,i) = b_cxst6[j][i];
  k_dtheta_cxst6_c.view_host()(i,j) = dtheta_cxst6_c[i][j]; k_dtheta_cxst6_c.view_host()(j,i) = dtheta_cxst6_c[j][i];

  k_AA_cxst1.view_host()(i,j) = AA_cxst1[i][j]; k_AA_cxst1.view_host()(j,i) = AA_cxst1[j][i];
  k_BB_cxst1.view_host()(i,j) = BB_cxst1[i][j]; k_BB_cxst1.view_host()(j,i) = BB_cxst1[j][i];

  k_k_cxst.template modify<LMPHostType>();
  k_cut_cxst_0.template modify<LMPHostType>();
  k_cut_cxst_c.template modify<LMPHostType>();
  k_cut_cxst_lo.template modify<LMPHostType>();
  k_cut_cxst_hi.template modify<LMPHostType>();
  k_cut_cxst_lc.template modify<LMPHostType>();
  k_cut_cxst_hc.template modify<LMPHostType>();
  k_b_cxst_lo.template modify<LMPHostType>();
  k_b_cxst_hi.template modify<LMPHostType>();
  k_cutsq_cxst_hc.template modify<LMPHostType>();

  k_a_cxst1.template modify<LMPHostType>();
  k_theta_cxst1_0.template modify<LMPHostType>();
  k_dtheta_cxst1_ast.template modify<LMPHostType>();
  k_b_cxst1.template modify<LMPHostType>();
  k_dtheta_cxst1_c.template modify<LMPHostType>();

  k_a_cxst4.template modify<LMPHostType>();
  k_theta_cxst4_0.template modify<LMPHostType>();
  k_dtheta_cxst4_ast.template modify<LMPHostType>();
  k_b_cxst4.template modify<LMPHostType>();
  k_dtheta_cxst4_c.template modify<LMPHostType>();

  k_a_cxst5.template modify<LMPHostType>();
  k_theta_cxst5_0.template modify<LMPHostType>();
  k_dtheta_cxst5_ast.template modify<LMPHostType>();
  k_b_cxst5.template modify<LMPHostType>();
  k_dtheta_cxst5_c.template modify<LMPHostType>();

  k_a_cxst6.template modify<LMPHostType>();
  k_theta_cxst6_0.template modify<LMPHostType>();
  k_dtheta_cxst6_ast.template modify<LMPHostType>();
  k_b_cxst6.template modify<LMPHostType>();
  k_dtheta_cxst6_c.template modify<LMPHostType>();

  k_AA_cxst1.template modify<LMPHostType>();
  k_BB_cxst1.template modify<LMPHostType>();

  // Sync to device
  k_k_cxst.template sync<DeviceType>();
  k_cut_cxst_0.template sync<DeviceType>();
  k_cut_cxst_c.template sync<DeviceType>();
  k_cut_cxst_lo.template sync<DeviceType>();
  k_cut_cxst_hi.template sync<DeviceType>();
  k_cut_cxst_lc.template sync<DeviceType>();
  k_cut_cxst_hc.template sync<DeviceType>();
  k_b_cxst_lo.template sync<DeviceType>();
  k_b_cxst_hi.template sync<DeviceType>();
  k_cutsq_cxst_hc.template sync<DeviceType>();

  k_a_cxst1.template sync<DeviceType>();
  k_theta_cxst1_0.template sync<DeviceType>();
  k_dtheta_cxst1_ast.template sync<DeviceType>();
  k_b_cxst1.template sync<DeviceType>();
  k_dtheta_cxst1_c.template sync<DeviceType>();

  k_a_cxst4.template sync<DeviceType>();
  k_theta_cxst4_0.template sync<DeviceType>();
  k_dtheta_cxst4_ast.template sync<DeviceType>();
  k_b_cxst4.template sync<DeviceType>();
  k_dtheta_cxst4_c.template sync<DeviceType>();

  k_a_cxst5.template sync<DeviceType>();
  k_theta_cxst5_0.template sync<DeviceType>();
  k_dtheta_cxst5_ast.template sync<DeviceType>();
  k_b_cxst5.template sync<DeviceType>();
  k_dtheta_cxst5_c.template sync<DeviceType>();

  k_a_cxst6.template sync<DeviceType>();
  k_theta_cxst6_0.template sync<DeviceType>();
  k_dtheta_cxst6_ast.template sync<DeviceType>();
  k_b_cxst6.template sync<DeviceType>();
  k_dtheta_cxst6_c.template sync<DeviceType>();

  k_AA_cxst1.template sync<DeviceType>();
  k_BB_cxst1.template sync<DeviceType>();

  // Register the COM screen cutoff for this pair: coaxial stacking acts at the
  // stacking site (COM +/- 0.34*nx on each atom); use the same conservative
  // 2*0.4 margin (>= 2*0.34) so a COM-distance screen never drops an interacting
  // pair. The npair fix takes the max over all consuming styles and type-pairs.
  if (fix_oxdna_npairKK) fix_oxdna_npairKK->request_screen_cutoff(cutone + 0.8);

  // "cutone" is "cut_cxst_hc[i][j]", sets the master list distance cutoff
  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairOxdna2CoaxstkKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (EFLAG) {
    if (eflag_atom) {
      const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * epair);
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
      } else {
        a_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(delx*fx);
    const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(dely*fy);
    const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(delz*fz);
    const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(delx*fy);
    const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(delx*fz);
    const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(dely*fz);

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
        if (NEWTON_PAIR || j < nlocal) {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
      } else {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
      }
    }

    if (vflag_atom) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
          a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
          a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
          a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
          a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
          a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
        if (NEWTON_PAIR || j < nlocal) {
        a_vatom(j,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        a_vatom(j,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        a_vatom(j,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        a_vatom(j,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        a_vatom(j,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        a_vatom(j,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
      } else {
        a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairOxdna2CoaxstkKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}


namespace LAMMPS_NS {
template class PairOxdna2CoaxstkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdna2CoaxstkKokkos<LMPHostType>;
#endif
}
